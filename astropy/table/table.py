from copy import deepcopy
import numpy as np
import collections
from itertools import count

from astropy.utils import OrderedDict
from structhelper import _drop_fields


class ArgumentError(Exception):
    pass


def insert_odict(d_old, position, key, value):
    '''Convenience function to insert values into an OrderedDict'''
    items = d_old.items()
    items.insert(position, (key, value))
    return OrderedDict(items)


def rename_odict(d_old, before, after):
    '''Convenience function to rename keys in an OrderedDict'''
    keys = [(after if key == before else key) for key in d_old.keys()]
    values = d_old.values()
    return OrderedDict(zip(keys, values))


class Column(np.ndarray):
    """Define a data column for use in a Table object.

    Parameters
    ----------
    name : str
        Column name and key for reference within Table
    data : list, ndarray or None
        Column data values
    dtype : numpy.dtype compatible value
        Data type for column
    shape : tuple or ()
        Dimensions of a single row element in the column data
    length : int or 0
        Number of row elements in column data
    description : str or None
        Full description of column
    units : str or None
        Physical units
    format : str or None
        Sprintf-style format string for outputting column values
    meta : dict or None
        Meta-data associated with the column

    Examples
    --------
    A Column can be created in two different ways:

    - Provide a ``data`` value and optionally a ``dtype`` value

      Examples::

        col = Column('name', data=[1, 2, 3])         # shape=(3,)
        col = Column('name', data=[[1, 2], [3, 4]])  # shape=(2, 2)
        col = Column('name', data=[1, 2, 3], dtype=float)
        col = Column('name', np.array([1, 2, 3]))
        col = Column('name', ['hello', 'world'])

      The ``dtype`` argument can be any value which is an acceptable
      fixed-size data-type initializer for the numpy.dtype() method.  See
      `http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html`_.
      Examples include:

      - Python non-string type (float, int, bool)
      - Numpy non-string type (e.g. np.float32, np.int64, np.bool)
      - Numpy.dtype array-protocol type strings (e.g. 'i4', 'f8', 'S15')

      If no ``dtype`` value is provide then the type is inferred using
      ``np.array(data)``.  When ``data`` is provided then the ``shape``
      and ``length`` arguments are ignored.

    - Provide zero or more of ``dtype``, ``shape``, ``length``

      Examples::

        col = Column('name')
        col = Column('name', dtype=int, length=10, shape=(3,4))

      The default ``dtype`` is ``np.float64`` and the default ``length`` is zero.
      The ``shape`` argument is the array shape of a single cell in the column.
      The default ``shape`` is () which means a single value in each element.
    """

    def __new__(cls, name, data=None,
                 dtype=None, shape=(), length=0,
                 description=None, units=None, format=None, meta=None):

        if data is None:
            self_data = np.zeros(length,
                                 dtype=(np.dtype(dtype).str, shape))
        else:
            self_data = np.asarray(data, dtype=dtype)

        self = self_data.view(cls)
        self.name = name
        self.units = units
        self.format = format
        self.description = description
        self.parent_table = None

        self.meta = OrderedDict()
        if meta is not None:
            self.meta.update(meta)

        return self

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        for attr in ('name', 'units', 'format', 'description'):
            val = getattr(obj, attr, None)
            setattr(self, attr, val)
        self.meta = deepcopy(getattr(obj, 'meta', {}))

    @property
    def data(self):
        return self.view(np.ndarray)

    def copy(self):
        """Return a copy of the current Column instance.
        """
        newdata = self.view(np.ndarray).copy()
        return Column(self.name, newdata, units=self.units, format=self.format,
                      description=self.description, meta=deepcopy(self.meta))

    @property
    def descr(self):
        """Array-interface compliant full description of the column.

        This returns a 3-tuple (name, type, shape) that can always be
        used in a structured array dtype definition.
        """
        return (self.name, self.dtype.str, self.shape[1:])

    # @property
    # def data(self):
    #     """Column data attribute.  Only available if parent data exist"""
    #     if self.parent_table is not None:
    #         return self.parent_table[self.name]
    #     else:
    #         return self

    def __repr__(self):
        s = "<Column name='{0} units='{1}' " \
            "format='{2}' description='{3}'>".format(
            self.name, self.units, self.format, self.description)
        return s

    def __eq__(self, c):
        attrs = ('name', 'units', 'dtype', 'format', 'description', 'meta')
        equal = all(getattr(self, attr) == getattr(c, attr) for attr in attrs)
        return equal

    def __ne__(self, other):
        return not self.__eq__(other)


class Table(object):
    """A class to represent tables of data.

    Parameters
    ----------
    data : numpy ndarray, dict, list, or Table
        Data to initialize table.
    names : list
        Specify column names (behavior depends on data type)
    dtypes : list
        Specify column data types (behavior depends on data type)
    meta : dict
        Metadata associated with the table

    Initialization
    --------------

    The Table object can be initialized with several different forms
    for the ``data`` argument.

    numpy ndarray (structured array)

        The base column names are the field names of the ``data`` structured
        array.  The ``names`` list (optional) can be used to select
        particular fields and/or reorder the base names.  The ``dtypes`` list
        (optional) must match the length of ``names`` and is used to
        override the existing ``data`` types.

    numpy ndarray (homogeneous)

        The ``data`` ndarray must be at least 2-dimensional, with the first
        (left-most) index corresponding to row number (table length) and the
        second index corresponding to column number (table width).  Higher
        dimensions get absorbed in the shape of each table cell.

        If provided the ``names`` list must match the "width" of the ``data``
        argument.  The default for ``names`` is to auto-generate column names
        in the form "col<N>".  If provided the ``dtypes`` list overrides the
        base column types and must match the length of ``names``.

    dict-like

        The keys of the ``data`` object define the base column names.  The
        corresponding values can be Column objects, numpy arrays, or list-like
        objects.  The ``names`` list (optional) can be used to select
        particular fields and/or reorder the base names.  The ``dtypes`` list
        (optional) must match the length of ``names`` and is used to override
        the existing or default data types.

    list-like

        Each item in the ``data`` list provides a column of data values and can
        can be a Column object, numpy array, or list-like object.  The
        ``names`` list defines the name of each column.  The names will be
        auto-generated if not provided (either from the ``names`` argument or
        by Column objects).  If provided the ``names`` argument must match the
        number of items in the ``data`` list.  The optional ``dtypes`` list
        will override the existing or default data types and must match
        ``names`` in length.

    None

        Initialize a zero-length table.  If ``names`` and optionally ``dtypes``
        are provided then the corresponding columns are created.
    """

    def __init__(self, data=None, names=None, dtypes=None, meta={}):
        self._data = None
        self.columns = OrderedDict()
        self.meta = deepcopy(meta)

        if isinstance(data, (list, tuple)):
            init_func = self._init_from_list

        elif isinstance(data, np.ndarray):
            if data.dtype.names:
                init_func = self._init_from_ndarray_struct
            else:
                init_func = self._init_from_ndarray_homog

        elif isinstance(data, dict):
            init_func = self._init_from_dict

        elif isinstance(data, Table):
            init_func = self._init_from_table

        elif data is None:
            if names is None:
                return
            else:
                init_func = self._init_from_list
                data = [] * len(names)
        else:
            raise ValueError('Data type {0} not allowed to init Table'
                             .format(type(data)))

        init_func(data, names, dtypes)

    def _check_names_dtypes(self, names, dtypes):
        """Make sure that names and dtypes are boths iterable and have
        the same length.
        """
        for inp_list, inp_str in ((dtypes, 'dtypes'), (names, 'names')):
            if not isinstance(inp_list, collections.Iterable):
                raise ValueError('{0} must be a list or None'.format(inp_str))

        if len(names) != len(dtypes):
            raise ValueError(
                'Arguments "names" and "dtypes" must match in length'
                .format(inp_str))

    def _init_from_list(self, data, names, dtypes):
        """Initialize table from a list of columns"""
        n_cols = len(data)
        if names is None:
            names = [None] * n_cols
        if dtypes is None:
            dtypes = [None] * n_cols
        self._check_names_dtypes(names, dtypes)

        cols = []
        for i_col, col, name, dtype in zip(count(), data, names, dtypes):
            def_name = 'col{0}'.format(i_col)
            if isinstance(col, Column):
                col = Column((name or col.name), col.data, dtype=dtype,
                             meta=deepcopy(col.meta))
            elif isinstance(col, (np.ndarray, collections.Iterable)):
                col = Column((name or def_name), col, dtype=dtype)
            else:
                raise ValueError('Elements in list initialization must be either '
                                 'Column or list-like')
            cols.append(col)

        self._init_from_cols(cols)

    def _init_from_ndarray_struct(self, data, names, dtypes):
        """Initialize table from an ndarray structured array"""
        if names is None:
            names = data.dtype.names
        if dtypes is None:
            dtypes = [None] * len(names)

        cols = [data[name] for name in names]
        self._init_from_list(cols, names, dtypes)

    def _init_from_ndarray_homog(self, data, names, dtypes):
        """Initialize table from an ndarray homogeneous array"""
        n_cols = data.shape[1]
        if names is None:
            names = [None] * n_cols
        if dtypes is None:
            dtypes = [None] * n_cols

        cols = [data[:, i] for i in range(n_cols)]
        self._init_from_list(cols, names, dtypes)

    def _init_from_dict(self, data, names, dtypes):
        if names is None:
            names = data.keys()
        if dtypes is None:
            dtypes = [None] * len(names)

        data_list = [data[name] for name in names]
        self._init_from_list(data_list, names, dtypes)

    def _init_from_table(self, data, names, dtypes):
        if names is None:
            names = data.colnames
        if dtypes is None:
            dtypes = [None] * len(names)

        cols = [data.columns[name] for name in names]
        self._init_from_list(cols, names, dtypes)
        self.meta = deepcopy(data.meta)

    def _init_from_cols(self, cols):
        dtypes = [col.descr for col in cols]

        lengths = set(len(col.data) for col in cols)
        if len(lengths) != 1:
            raise ValueError(
                'Inconsistent data column lengths: {0}'.format(lengths))

        self._data = np.ndarray(lengths.pop(), dtype=dtypes)
        self.columns = OrderedDict()
        for col in cols:
            self._data[col.name] = col.data
            newcol = Column(col.name, self._data[col.name], units=col.units,
                            format=col.format, description=col.description,
                            meta=deepcopy(col.meta))
            self.columns[col.name] = newcol

    def __repr__(self):
        s = "<Table "
        s += "name='{0}' ".format(self.name)
        s += "rows='{0}' ".format(self.__len__())
        s += "columns='{0}'>".format(len(self.columns))
        return s

    def __getitem__(self, item):
        try:
            return self._data[item]
        except (ValueError, KeyError, TypeError):
            raise KeyError("Column {0} does not exist".format(item))
        except:
            # Bad index raises "IndexError: index out of bounds", but also
            # re-raise any other exception
            raise

    def __setitem__(self, item, value):
        try:
            self._data[item] = value
        except (ValueError, KeyError, TypeError):
            raise KeyError("Column {0} does not exist".format(item))
        except:
            raise

    @property
    def colnames(self):
        return self.columns.keys()

    def keys(self):
        return self.columns.keys()

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data)

    def index_column(self, name):
        """
        Return the index of column ``name``.
        """
        try:
            return self.colnames.index(name)
        except ValueError:
            raise ValueError("Column {0} does not exist".format(name))

    def add_column(self, col, index=None):
        """
        Add a new Column object ``col`` to the table.  If ``index``
        is supplied then insert column before ``index`` position
        in the list of columns, otherwise append column to the end
        of the list.

        Parameters
        ----------
        col : Column
            Column object to add.
        index : int or None
            Insert column before this position or at end (default)
        """
        if index is None:
            index = len(self.columns)
        self.add_columns([col], [index])

    def add_columns(self, cols, indexes=None):
        """
        Add a list of new Column objects ``cols`` to the table.  If a
        corresponding list of ``indexes`` is supplied then insert column before
        each ``index`` position in the *original* list of columns, otherwise
        append columns to the end of the list.

        Parameters
        ----------
        cols : list of Columns
            Column objects to add.
        indexes : list of ints or None
            Insert column before this position or at end (default)
        """
        if indexes is None:
            indexes = [len(self.columns)] * len(cols)
        elif len(indexes) != len(cols):
            raise ValueError('Number of indexes must match number of cols')

        if self._data is None:
            # No existing table data, init from cols
            newcols = cols
        else:
            new_indexes = range(len(newcols) + 1)
            newcols = self.columns.values()
            for col, index in zip(cols, indexes):
                i = new_indexes.index(index)
                new_indexes.insert(i, None)
                newcols.insert(i, col)

        self._init_from_cols(newcols)

    def remove_column(self, name):
        """
        Remove a column from the table

        Parameters
        ----------
        name : str
            Name of column to remove
        """

        self.remove_columns([name])

    def remove_columns(self, names):
        '''
        Remove several columns from the table

        Parameters
        ----------
        names : list
            A list containing the names of the columns to remove
        '''

        for name in names:
            if name not in self.columns:
                raise KeyError("Column {0} does not exist".format(name))

        for name in names:
            self.columns.pop(name)

        self._data = _drop_fields(self._data, names)

    def keep_columns(self, names):
        '''
        Keep only the columns specified (remove the others)

        Parameters
        ----------
        names : list
            A list containing the names of the columns to keep. All other
            columns will be removed.
        '''

        if isinstance(names, basestring):
            names = [names]

        for name in names:
            if name not in self.columns:
                raise KeyError("Column {0} does not exist".format(name))

        remove = list(set(self.keys()) - set(names))

        self.remove_columns(remove)

    def rename_column(self, name, new_name):
        '''
        Rename a column

        Parameters
        ----------
        name : str
            The current name of the column.
        new_name : str
            The new name for the column
        '''

        if name not in self.keys():
            raise KeyError("Column {0} does not exist".format(name))

        if new_name in self.keys():
            raise KeyError("Column {0} already exists".format(new_name))

        pos = self.columns.keys().index(name)
        self._data.dtype.names = (self.keys()[:pos] + [new_name, ]
                                  + self.keys()[pos + 1:])

        self.columns = rename_odict(self.columns, name, new_name)

    def add_row(self, vals=None):
        """Add a new row to the end of the table.

        The ``vals`` argument can be:

        sequence (e.g. tuple or list)
            Column values in the same order as table columns.
        mapping (e.g. dict)
            Keys corresponding to column names.  Missing values will be
            filled with np.zeros for the column dtype.
        None
            All values filled with np.zeros for the column dtype.

        Parameters
        ----------
        vals : tuple, list, dict or None
            Use the specified values in the new row
        """
        newlen = len(self._data) + 1
        self._data.resize((newlen,), refcheck=False)

        if isinstance(vals, collections.Mapping):
            row = self._data[-1]
            for name, val in vals.items():
                try:
                    row[name] = val
                except IndexError:
                    raise ValueError("No column {0} in table".format(name))

        elif isinstance(vals, collections.Iterable):
            if len(self.columns) != len(vals):
                raise ValueError('Mismatch between number of vals and columns')
            if not isinstance(vals, tuple):
                vals = tuple(vals)
            self._data[-1] = vals

        elif vals is not None:
            raise TypeError('Vals must be an iterable or mapping or None')

    def sort(self, keys):
        '''
        Sort the table according to one or more keys. This operates
        on the existing table (and does not return a new table).

        Parameters
        ----------
        keys : str or list of str
            The key(s) to order the table by
        '''
        if type(keys) is not list:
            keys = [keys]
        self._data.sort(order=keys)

from __future__ import print_function

from copy import deepcopy
import numpy as np
import collections
from itertools import count

from astropy.utils import OrderedDict
from .structhelper import _drop_fields

try:
    unicode
except NameError:
    unicode = basestring = str


class TableColumns(OrderedDict):
    def __init__(self, table=None, cols={}):
        self.table = table  # the parent table containing these columns
        if isinstance(cols, (list, tuple)):
            cols = [(col.name, col) for col in cols]
        super(TableColumns, self).__init__(cols)

    def __getitem__(self, item):
        """Get times from a TableColumns object.
        ::

          tc = TableColumns(cols=[Column('a'), Column('b'), Column('c')])
          tc['a']  # Column('a')
          tc[1] # Column('b')
          tc['a', 'b'] # <TableColumns names=('a', 'b')>
          tc[1:3] # <TableColumns names=('b', 'c')>
        """
        if isinstance(item, basestring):
            return OrderedDict.__getitem__(self, item)
        elif isinstance(item, int):
            return self.values()[item]
        elif isinstance(item, tuple):
            return TableColumns(self.table, [self[x] for x in item])
        elif isinstance(item, slice):
            return TableColumns(self.table,
                                [self[x] for x in self.keys()[item]])
        else:
            raise IndexError('Illegal key or index value for TableColumns '
                             'object')

    def __repr__(self):
        names = ("'{0}'".format(x) for x in self.keys())
        return "<TableColumns names=({0})>".format(
            ",".join(names))

    def _rename_column(self, name, new_name):
        if new_name in self:
            raise KeyError("Column {0} already exists".format(new_name))

        mapper = {name: new_name}
        new_names = [mapper.get(name, name) for name in self]
        cols = self.values()
        self.clear()
        self.update(zip(new_names, cols))

    def keys(self):
        return list(OrderedDict.keys(self))

    def values(self):
        return list(OrderedDict.values(self))


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

      The default ``dtype`` is ``np.float64`` and the default ``length`` is
      zero.  The ``shape`` argument is the array shape of a single cell in the
      column.  The default ``shape`` is () which means a single value in each
      element.
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
        self._dtype = self_data.dtype
        self._name = name
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
        for attr in ('units', 'format', 'description'):
            val = getattr(obj, attr, None)
            setattr(self, attr, val)
        self.meta = deepcopy(getattr(obj, 'meta', {}))

    def _get_name(self):
        return self._name

    def _set_name(self, val):
        if self.parent_table is not None:
            table = self.parent_table
            table.columns._rename_column(self.name, val)
            table._data.dtype.names = table.columns.keys()

        self._name = val

    name = property(_get_name, _set_name)

    @property
    def dtype(self):
        return self._dtype

    @property
    def data(self):
        return self.view(np.ndarray)

    def copy(self, data=None):
        """Return a copy of the current Column instance.
        """
        if data is None:
            data = self.view(np.ndarray).copy()

        return Column(self.name, data, units=self.units, format=self.format,
                      description=self.description, meta=deepcopy(self.meta))

    @property
    def descr(self):
        """Array-interface compliant full description of the column.

        This returns a 3-tuple (name, type, shape) that can always be
        used in a structured array dtype definition.
        """
        return (self.name, self.dtype.str, self.shape[1:])

    def __repr__(self):
        s = "<Column name='{0} units='{1}' format='{2}' " \
            "description='{3}'>\n{4}".format(self.name, self.units,
                 self.format, self.description, repr(self.data))
        return s

    def __eq__(self, col):
        if isinstance(col, Column):
            attrs = ('name', 'units', 'dtype', 'format', 'description', 'meta')
            equal = all(getattr(self, x) == getattr(col, x) for x in attrs)
        else:
            equal = self.data == col
        return equal

    def __ne__(self, other):
        return not self.__eq__(other)


class Row(object):
    """A class to represent one row of a Table object.

    A Row object is returned when a Table object is indexed with an integer
    or when iterating over a table::

      >>> table = Table([(1, 2), (3, 4)], names=('a', 'b'))
      >>> row = table[1]
      >>> row
      <Row 1 of table
       values=(2, 4)
       dtype=[('a', '<i8'), ('b', '<i8')]>
      >>> row['a']
      2
      >>> row[1]
      4
    """

    def __init__(self, table, index):
        self.table = table
        self.index = index
        self.data = table._data[index]

    def __getitem__(self, item):
        return self.data[item]

    def __setitem__(self, item, val):
        self.data[item] = val

    @property
    def meta(self):
        return self.table.meta

    @property
    def columns(self):
        return self.table.columns

    @property
    def colnames(self):
        return self.table.colnames

    @property
    def dtype(self):
        return self.data.dtype

    def __repr__(self):
        return "<Row {0} of table\n values={1}\n dtype={2}>".format(
            self.index, self.data, self.dtype)


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
        self.columns = TableColumns(self)
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
                raise ValueError('Elements in list initialization must be '
                                 'either Column or list-like')
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
        columns = TableColumns(self)
        dtypes = [col.descr for col in cols]

        lengths = set(len(col.data) for col in cols)
        if len(lengths) != 1:
            raise ValueError(
                'Inconsistent data column lengths: {0}'.format(lengths))

        data = np.ndarray(lengths.pop(), dtype=dtypes)
        for col in cols:
            data[col.name] = col.data
            newcol = Column(col.name, data[col.name], units=col.units,
                            format=col.format, description=col.description,
                            meta=deepcopy(col.meta))
            columns[col.name] = newcol
            newcol.parent_table = self

        self.columns = columns
        self._data = data

    def __repr__(self):
        names = ("'{0}'".format(x) for x in self.colnames)
        s = "<Table rows={0} names=({1})>\n{2}".format(
            self.__len__(),  ','.join(names), repr(self._data))
        return s

    def __getitem__(self, item):
        if isinstance(item, basestring):
            return self.columns[item]
        elif isinstance(item, int):
            # XXX should return Row instance
            return Row(self, item)
        elif isinstance(item, tuple):
            return Table([self[x] for x in item], meta=deepcopy(self.meta))
        else:
            # XXX Losing other column attrs like description format etc.
            # XXX Error checking?
            return Table(self._data[item], names=self.colnames,
                         meta=deepcopy(self.meta))

    def __setitem__(self, item, value):
        try:
            self._data[item] = value
        except (ValueError, KeyError, TypeError):
            raise KeyError("Column {0} does not exist".format(item))
        except:
            raise

    @property
    def colnames(self):
        return list(self.columns.keys())

    def keys(self):
        return list(self.columns.keys())

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
            newcols = list(self.columns.values())
            new_indexes = list(range(len(newcols) + 1))
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
        Rename a column.

        This can also be done directly with by setting the ``name`` attribute
        for a column::

          table[name].name = new_name

        Parameters
        ----------
        name : str
            The current name of the column.
        new_name : str
            The new name for the column
        '''

        if name not in self.keys():
            raise KeyError("Column {0} does not exist".format(name))

        self.columns[name].name = new_name

#         pos = self.columns.keys().index(name)
#         self._data.dtype.names = (self.keys()[:pos] + [new_name, ]
#                                   + self.keys()[pos + 1:])

#         keys = [(after if key == before else key) for key in d_old.keys()]
#         values = d_old.values()
#         return OrderedDict(zip(keys, values))

#         self.columns = rename_odict(self.columns, name, new_name)

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

        # Add_row() probably corrupted the Column views of self._data.  Rebuild
        # self.columns.  Col.copy() takes an optional data reference that it
        # uses in the copy.
        cols = [c.copy(self._data[c.name]) for c in self.columns.values()]
        self.columns = TableColumns(self, cols)

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

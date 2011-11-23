import copy
import numpy as np
import collections

from astropy.utils import OrderedDict
from structhelper import _append_field, _drop_fields


class ArgumentError(Exception):
    pass


def insert_odict(d_old, position, key, value):
    '''Convenience function to insert values into an OrderedDict'''
    d_new = OrderedDict()
    for i, k in enumerate(d_old):
        if position == i:
            d_new[key] = value
        d_new[k] = d_old[k]
    if position == len(d_old):
        d_new[key] = value
    return d_new


def rename_odict(d_old, before, after):
    '''Convenience function to rename keys in an OrderedDict'''
    d_new = OrderedDict()
    for k in d_old:
        if k == before:
            d_new[after] = d_old[k]
        else:
            d_new[k] = d_old[k]
    return d_new


class Column(object):
    """Define a data column for use in a Table object.

    Parameters
    ----------
    name : str
        Column name and key for reference within Table
    data : list, ndarray or None
        Column data values
    datatype : see examples for type
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

    *Provide a `data` value and optionally a `datatype` value*
    ::

      col = Column('name', data=[1, 2, 3])         # shape=(3,)
      col = Column('name', data=[[1, 2], [3, 4]])  # shape=(2, 2)
      col = Column('name', data=[1, 2, 3], datatype=float)  # float type
      col = Column('name', np.array([1, 2, 3]))
      col = Column('name', ['hello', 'world'])

    The `datatype` value can be one of the following (see
    `http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html`_):

    - Python non-string type (float, int, bool)
    - Numpy non-string type (e.g. np.float32, np.int64, np.bool)
    - Numpy.dtype array-protocol type strings (e.g. 'i4', 'f8', 'S15')

    If no `datatype` value is provide then the type is inferred using
    `np.array(data)`.  When `data` is provided then the `shape` and `length`
    args are ignored.

    *Provide zero or more of `datatype`, `shape`, `length`*
    ::

      col = Column('name')
      col = Column('name', datatype=int, length=10, shape=(3,4))

    The default `datatype` is `np.float` and the default `length` is zero.
    The `shape` argument is the array shape of a single cell in the column.
    The default `shape` is () which means a single value in each element.
    """

    def __init__(self, name, data=None,
                 datatype=None, shape=(), length=0,
                 description=None, units=None, format=None, meta=None):

        self.name = name
        self.units = units
        self.format = format
        self.description = description
        self.datatype = datatype
        self.parent_table = None

        self.meta = OrderedDict()
        if meta is not None:
            self.meta.update(meta)

        if data is None:
            self._data = np.zeros(length,
                                  dtype=(datatype or np.float, shape))
        else:
            try:
                dtype = (datatype or data.dtype, data.shape[1:])
            except AttributeError:
                data = np.array(data)
                dtype = (datatype or data.dtype, data.shape[1:])
            self._data = np.ndarray(len(data), dtype=dtype)
            self._data[:] = data

    def clear_data(self):
        """Set the internal column data attribute to None in order
        to release the reference.  This should typically be done
        in combination with setting the parent_table attribute so
        that self.data returns the desired column data.
        """
        self._data = None

    def copy(self):
        """Return a copy of the current Column instance.
        """
        # Use a minimal constructor then manually copy attributes
        newcol = Column(self.name)
        for attr in ('units', 'format', 'description', 'datatype',
                     'parent_table', '_data'):
            val = getattr(self, attr)
            setattr(newcol, attr, val)
        newcol.meta = copy.deepcopy(self.meta)

        return newcol

    @property
    def dtype(self):
        """Data type attribute.  This checks that column data reference has
        been made."""
        if self.data is None:
            raise ValueError('No column data reference available '
                             'so dtype is undefined')
        return self.data.dtype

    @property
    def data(self):
        """Column data attribute.  Only available if parent data exist"""
        if self.parent_table is not None:
            return self.parent_table[self.name]
        else:
            return self._data

    def __repr__(self):
        s = "<Column name='{0} units='{1}' " \
            "format='{2}' description='{3}'>".format(
            self.name, self.units, self.format, self.description)
        return s

    def __eq__(self, c):
        attrs = ('name', 'units', 'datatype', 'format', 'description', 'meta')
        equal = all(getattr(self, attr) == getattr(c, attr) for attr in attrs)
        return equal

    def __ne__(self, other):
        return not self.__eq__(other)


class Table(object):
    '''A class to represent tables of data'''

    def __init__(self, cols=None, name=None):
        self.name = name
        self._data = None
        self.columns = OrderedDict()
        if cols is not None:
            self._new_from_cols(cols)

    def _new_from_cols(self, cols):
        dtype = [(col.name, col.data.dtype, col.data.shape[1:])
                 for col in cols]

        lengths = set(len(col.data) for col in cols)
        if len(lengths) != 1:
            raise ValueError(
                'Inconsistent data column lengths: {0}'.format(lengths))

        self._data = np.ndarray(lengths.pop(), dtype=dtype)
        for col in cols:
            self._data[col.name] = col.data
            if col.parent_table is not None:
                col = col.copy()
            col.parent_table = self  # see same in insert_column()
            col.clear_data()
            self.columns[col.name] = col

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

        dtype = (col.name, col.data.dtype, col.data.shape[1:])

        if self._data is None:
            # Table has no existing data so make a new structured array
            # with a single column and then copy the column data.
            self._data = np.ndarray(len(col.data), dtype=[dtype])
            self._data[col.name] = col.data
        else:
            if len(col.data) != len(self._data):
                raise ValueError(
                    "Column data length does not match table length")

            self._data = _append_field(self._data, col.data, dtype=dtype,
                                       position=index)

        # If the user-supplied col is already part of a table then copy.
        if col.parent_table is not None:
            col = col.copy()

        # Now that column data are copied into the Table _data table set
        # the column parent_table to self and clear the data reference.
        col.parent_table = self
        col.clear_data()

        self.columns = insert_odict(self.columns, index, col.name, col)

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

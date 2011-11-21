import numpy as np

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
    """A class to contain information about columns"""

    def __init__(self, name=None, datatype=None, shape=tuple(),
                 units=None, format=None, description=None, length=0,
                 meta=None, data=None):

        self.name = name
        self.units = units
        self.format = format
        self.description = description
        self.datatype = datatype
        self._parent_table = None

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
                dtype = (data.dtype, data.shape[1:])
            self._data = np.ndarray(len(data), dtype=dtype)
            self._data[:] = data

    def _get_parent_table(self):
        return self._parent_table

    def _set_parent_table(self, val):
        self._parent_table = val
        self._data = None

    parent_table = property(_get_parent_table, _set_parent_table)

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
        if self._parent_table is not None:
            return self._parent_table[self.name]
        else:
            return self._data

    def __repr__(self):
        s = "<Column name='{0} units='{1}' format='{2}' description='{3}'>".format(
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

    def __init__(self, name=None, length=None):
        self.name = name
        # self._length = length
        self._data = None
        self.columns = OrderedDict()

        # ??? Maybe make masking read-only, need to understand scenarios when a
        # table can convert from normal to masked etc.  This is TBD.
        self.masked = False

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

    def add_column(self, *args, **kwargs):
        self.append_column(*args, **kwargs)

    def append_column(self, name=None, data=None, datatype=None, shape=tuple(),
                      units=None, format=None, description=None, length=0,
                      meta=None):
        """
        Add a new Column object ``column`` after the last existing column.
        """
        self.insert_column(len(self.columns),
                           name=name, data=data, datatype=datatype,
                           shape=shape, units=units, format=format,
                           description=description,
                           length=length, meta=meta)

    def insert_column(self, index, name=None, data=None, datatype=None,
                      shape=tuple(), units=None, format=None, description=None,
                      length=0, meta=None):
        """
        Insert a new Column object ``column`` at given ``index`` position.
        """
        # Once self._data table is defined then length is set by table length
        if self._data is not None:
            length = len(self._data)

        col = Column(name=name, datatype=datatype, shape=shape,
                     units=units, format=format, description=description,
                     length=length, meta=meta, data=data)

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

        # Now that column data are copied into the Table _data table set
        # the column parent_table to self.  This causes the column to
        # drop its reference to data and makes col.data return
        # col.parent_table[col.name].
        col.parent_table = self

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

        self._data = _drop_fields(self._data, names)  # XXX Doesn't set mask kwarg 

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

    def rename_column(self, before, after):
        '''
        Rename a column

        Parameters
        ----------
        before : str
            The current name of the column.
        after : str
            The new name for the column
        '''

        if before not in self.keys():
            raise KeyError("Column {0} does not exist".format(before))

        if after in self.keys():
            raise KeyError("Column {0} already exists".format(after))

        pos = self.columns.keys().index(before)
        self._data.dtype.names = self.keys()[:pos] + [after, ] + self.keys()[pos + 1:]

        self.columns = rename_odict(self.columns, before, after)

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

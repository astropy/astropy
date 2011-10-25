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
    '''A class to contain information about columns'''

    def __init__(self, name=None, dtype=None, units=None, format=None,
                 description=None, meta=None):

        self.name = name
        self.units = units
        self.format = format
        self.description = description
        self.meta = OrderedDict()
        self._dtype = dtype

        if meta is not None:
            self.meta.update(meta)

    @property
    def dtype(self):
        """Read-only data type attribute"""
        return self._dtype

    def __repr__(self):
        s = "<Column name='{0} units='{1}' format='{2}' description='{3}'>".format(
            self.name, self.units, self.format, self.description)
        return s

    def __eq__(self, c):
        attrs = ('name', 'units', 'dtype', 'format', 'description', 'meta')
        equal = all(getattr(self, attr) == getattr(c, attr) for attr in attrs)
        return equal

    def __ne__(self, other):
        return not self.__eq__(other)


class Table(object):
    '''A class to represent tables of data'''

    def __init__(self, name=None, length=None):
        self.name = name
        self._length = length
        self._data = None
        self.columns = OrderedDict()

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

    def insert_column(self, index, column):
        """
        Insert a new Column object ``column`` at given ``index`` position.
        """
        pass

    def append_column(self, column):
        """
        Add a new Column object ``column`` after the last existing column.
        """
        self.insert_column(len(self.columns), column)

    def add_column(self, name, data=None, dtype=None, shape=None,
                   units=None, format=None, description=None, meta=None,
                   before=None, after=None, position=None):
        '''
        Add a column to the table

        Parameters
        ----------
        name : string
            The name of the column to add.
        data : ~numpy.array
            The column data.
        dtype : type or np.dtype
            Type to convert the data to. This is the equivalent to the `dtype`
            argument in numpy.array.
        shape : tuple, optional
            Shape for the cells
        units : str, optional
            The units of the values in the column.
        format : str, optional
            The format to use for ASCII printing.
        description : str, optional
            A description of the content of the column.
        meta : dict or OrderedDict, optional
            A dictionary of metadata for the column.
        before : str, optional
            Column before which the new column should be inserted.
        after : str, optional
            Column after which the new column should be inserted.
        position : int, optional
            Position at which the new column should be inserted.
        '''

        # Convert data to requested type
        if data is None:
            if dtype is None:
                raise ArgumentError("dtype is required if data is not specified")
            else:
                if self.__len__() > 0:
                    length = self.__len__()
                else:
                    length = self._length

                if length is not None:
                    if shape:
                        data = np.zeros((length, ) + shape, dtype=dtype)
                    else:
                        data = np.zeros(length, dtype=dtype)
                else:
                    raise ArgumentError("If you want to add an empty column to an empty table, "
                                        "you need to specify the length of the table when "
                                        "initializing the Table instance.")
        else:
            if self._length is not None and len(data) != self._length:
                raise ValueError("data length does not match length specified when initializing Table")
            if self._data is not None and len(data) != self.__len__():
                raise ValueError("data length does not match table length")
            data = np.array(data, dtype=dtype)

        # Create Column instance to describe the column
        column = Column(name=name, dtype=data.dtype, units=units,
                        format=format, description=description, meta=meta)

        if (before, after, position).count(None) < 2:
            raise ArgumentError("Only one of before/after/position can be specified")

        if before is not None:
            if before in self.keys():
                position = self.keys().index(before)
            else:
                raise KeyError("Column {0} does not exist".format(before))

        if after is not None:
            if after in self.keys():
                position = self.keys().index(after) + 1
            else:
                raise KeyError("Column {0} does not exist".format(after))

        if position is not None:
            if position > len(self.columns):
                raise ValueError("position cannot be larger than the number of columns")
            if position < 0:
                raise ValueError("position cannot be negative")

        if data.ndim > 1:
            dtype_new = (name, data.dtype, data.shape[1:])
        else:
            dtype_new = (name, data.dtype)

        # Add the column to the structured array
        if self._data is None:
            self._data = np.array(zip(data), dtype=[dtype_new])
        else:
            self._data = _append_field(self._data, data, dtype=dtype_new,
                                       position=position)

        # Add the metadata to the columns attribute
        if not np.equal(position, None):
            self.columns = insert_odict(self.columns, position, name, column)
        else:
            self.columns[name] = column

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

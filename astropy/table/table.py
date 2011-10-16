import numpy as np

from astropy.utils import OrderedDict
from structhelper import _append_field


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


class Column(object):
    '''A class to contain information about columns'''

    def __init__(self, name=None, dtype=None, units=None, format=None,
                 description=None):

        self.name = name
        self.units = units
        self.format = format
        self.description = description
        self.meta = OrderedDict()

        object.__setattr__(self, 'dtype', dtype)

    def __setattr__(self, attribute, value):
        if attribute == 'dtype':
            raise Exception("Cannot change dtype through Column")
        else:
            object.__setattr__(self, attribute, value)

    def __repr__(self):
        s = "<Column "
        s += "name='{0}' ".format(self.name)
        s += "units='{0}' ".format(self.units)
        s += "format='{0}' ".format(self.format)
        s += "description='{0}'>".format(self.description)
        return s

    def __eq__(self, c):
        if self.name != c.name:
            return False
        if self.units != c.units:
            return False
        if self.dtype != c.dtype:
            return False
        if self.format != c.format:
            return False
        if self.description != c.description:
            return False
        if self.meta.keys != c.meta.keys:
            return False
        if self.meta.values != c.meta.values:
            return False
        return True

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
        if self._data is None or item not in self.columns:
            raise KeyError("Column {0} does not exist".format(item))
        else:
            return self._data[item]

    def __setitem__(self, item, value):
        if self._data is None or item not in self.columns:
            raise KeyError("Column {0} does not exist".format(item))
        else:
            self._data[item] = value

    def keys(self):
        return self.columns.keys

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data)

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
                    raise ArgumentError("If you want to add an empty column to an empty table, you need to specify the length of the table when initializing the Table instance.")
        else:
            if self._length is not None and len(data) != self._length:
                raise ValueError("data length does not match length specified when initializing Table")
            if self._data is not None and len(data) != self.__len__():
                raise ValueError("data length does not match table length")
            data = np.array(data, dtype=dtype)

        # Create Column instance to describe the column
        column = Column(name=name, dtype=data.dtype, units=units,
                        format=format, description=description)

        if (before, after, position).count(None) < 2:
            raise ArgumentError("Only one of before/after/position can be specified")

        if before is not None:
            if before in self.columns.keys():
                position = self.columns.keys().index(before)
            else:
                raise KeyError("Column {0} does not exist".format(before))

        if after is not None:
            if after in self.columns.keys():
                position = self.columns.keys().index(after) + 1
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

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import abc
import collections
import sys
from copy import deepcopy
import functools
import warnings

import numpy as np
from numpy import ma

from ..units import Unit
from .. import log
from ..utils import OrderedDict, isiterable, deprecated, deprecated_attribute
from .pprint import _pformat_table, _pformat_col, _pformat_col_iter, _more_tabcol
from ..utils.console import color_print
from ..config import ConfigurationItem
from ..io import registry as io_registry
from . import operations

# Python 2 and 3 source compatibility
try:
    unicode
except NameError:
    unicode = basestring = str

NUMPY_LT_1P5 = [int(x) for x in np.__version__.split('.')[:2]] < [1, 5]

AUTO_COLNAME = ConfigurationItem(
    'auto_colname', 'col{0}',
    'The template that determines the name of a column if it cannot be '
    'determined. Uses new-style (format method) string formatting')

ERROR_COLUMN_ARGS_MESSAGE = """
The first argument to {class_name} is the string {first_arg}, which was probably intended
as the column name.  Starting in Astropy 0.3 the argument order for initializing
a {class_name} object is {class_name}(data=None, name=None, ...)."""


def _check_column_new_args(func):
    """
    Decorator for transition from 0.2 arg order (name, data, ..) to 0.3 order (data,
    name, ...).  Check if user provided a string as the first arg (note that a string
    cannot be valid as ``data``).  Raise an error with a useful message.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if len(args) > 1 and isinstance(args[1], basestring):
            cls = args[0]  # Column or MaskedColumn class from __new__(cls, ..)
            raise ValueError(ERROR_COLUMN_ARGS_MESSAGE.format(class_name=cls.__name__,
                                                              first_arg=repr(args[1])))
        return func(*args, **kwargs)
    return wrapper


def _auto_names(n_cols):
    return [AUTO_COLNAME().format(i) for i in range(n_cols)]


class TableColumns(OrderedDict):
    """OrderedDict subclass for a set of columns.

    This class enhances item access to provide convenient access to columns
    by name or index, including slice access.  It also handles renaming
    of columns.

    The initialization argument ``cols`` can be any structure that is valid
    for initializing a Python dict.  This includes a dict, list of
    (key, val) tuple pairs, list of [key, val] lists, etc.

    Parameters
    ----------
    cols : dict, list, tuple; optional
        Column objects as data structure that can init dict (see above)
    """

    def __init__(self, cols={}):
        if isinstance(cols, (list, tuple)):
            cols = [(col.name, col) for col in cols]
        super(TableColumns, self).__init__(cols)

    def __getitem__(self, item):
        """Get items from a TableColumns object.
        ::

          tc = TableColumns(cols=[Column(name='a'), Column(name='b'), Column(name='c')])
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
            return TableColumns([self[x] for x in item])
        elif isinstance(item, slice):
            return TableColumns([self[x] for x in self.keys()[item]])
        else:
            raise IndexError('Illegal key or index value for TableColumns '
                             'object')

    def __repr__(self):
        names = ("'{0}'".format(x) for x in self.keys())
        return "<TableColumns names=({0})>".format(",".join(names))

    def _rename_column(self, name, new_name):
        if new_name in self:
            raise KeyError("Column {0} already exists".format(new_name))

        mapper = {name: new_name}
        new_names = [mapper.get(name, name) for name in self]
        cols = self.values()
        self.clear()
        self.update(zip(new_names, cols))

    # Define keys and values for Python 2 and 3 source compatibility
    def keys(self):
        return list(OrderedDict.keys(self))

    def values(self):
        return list(OrderedDict.values(self))


class BaseColumn(object):

    __metaclass__ = abc.ABCMeta

    def __array_finalize__(self, obj):
        # Obj will be none for direct call to Column() creator
        if obj is None:
            return

        # Self was created from template (e.g. obj[slice] or (obj * 2))
        # or viewcast e.g. obj.view(Column).  In either case we want to
        # init Column attributes for self from obj if possible.
        self.parent_table = None
        for attr in ('name', 'unit', 'format', 'description'):
            val = getattr(obj, attr, None)
            setattr(self, attr, val)
        self.meta = deepcopy(getattr(obj, 'meta', {}))

    def __array_wrap__(self, out_arr, context=None):
        """
        __array_wrap__ is called at the end of every ufunc.

        If the output is the same shape as the column then call the standard
        ndarray __array_wrap__ which will return a Column object.

        If instead the output shape is different (e.g. for reduction ufuncs
        like sum() or mean()) then return the output viewed as a standard
        np.ndarray.  The "[()]" selects everything, but also converts a zero
        rank array to a scalar.  For some reason np.sum() returns a zero rank
        scalar array while np.mean() returns a scalar.  So the [()] is needed
        for this case.
        """
        if self.shape == out_arr.shape:
            return np.ndarray.__array_wrap__(self, out_arr, context)
        else:
            return out_arr.view(np.ndarray)[()]

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, val):
        if self.parent_table is not None:
            table = self.parent_table
            table.columns._rename_column(self.name, val)
            table._data.dtype.names = table.columns.keys()
            if table.masked:
                table._data.mask.dtype.names = table.columns.keys()

        self._name = val

    @property
    def descr(self):
        """Array-interface compliant full description of the column.

        This returns a 3-tuple (name, type, shape) that can always be
        used in a structured array dtype definition.
        """
        return (self.name, self.dtype.str, self.shape[1:])

    def __repr__(self):
        unit = None if self.unit is None else str(self.unit)
        out = "<{0} name={1} unit={2} format={3} " \
            "description={4}>\n{5}".format(
            self.__class__.__name__,
            repr(self.name), repr(unit),
            repr(self.format), repr(self.description), repr(self.data))

        return out

    def iter_str_vals(self):
        """
        Return an iterator that yields the string-formatted values of this
        column.

        Returns
        -------
        str_vals : iterator
            Column values formatted as strings
        """
        # pprint._pformat_col_iter(col, max_lines, show_name, show_unit, outs)
        # Iterate over formatted values with no max number of lines, no column
        # name, no unit, and ignoring the returned header info in outs.
        for str_val in _pformat_col_iter(self, -1, False, False, {}):
            yield str_val

    def attrs_equal(self, col):
        """Compare the column attributes of ``col`` to this object.

        The comparison attributes are: name, unit, dtype, format, description,
        and meta.

        Parameters
        ----------
        col: Column
            Comparison column

        Returns
        -------
        equal: boolean
            True if all attributes are equal
        """
        if not isinstance(col, BaseColumn):
            raise ValueError('Comparison `col` must be a Column or MaskedColumn object')

        attrs = ('name', 'unit', 'dtype', 'format', 'description', 'meta')
        equal = all(getattr(self, x) == getattr(col, x) for x in attrs)

        return equal

    def pformat(self, max_lines=None, show_name=True, show_unit=False):
        """Return a list of formatted string representation of column values.

        If no value of `max_lines` is supplied then the height of the screen
        terminal is used to set `max_lines`.  If the terminal height cannot be
        determined then the default will be determined using the
        `astropy.table.pprint.MAX_LINES` configuration item. If a negative
        value of `max_lines` is supplied then there is no line limit applied.

        Parameters
        ----------
        max_lines : int
            Maximum lines of output (header + data rows)

        show_name : bool
            Include column name (default=True)

        show_unit : bool
            Include a header row for unit (default=False)

        Returns
        -------
        lines : list
            List of lines with header and formatted column values

        """
        lines, n_header = _pformat_col(self, max_lines, show_name, show_unit)
        return lines

    def pprint(self, max_lines=None, show_name=True, show_unit=False):
        """Print a formatted string representation of column values.

        If no value of `max_lines` is supplied then the height of the screen
        terminal is used to set `max_lines`.  If the terminal height cannot be
        determined then the default will be determined using the
        `astropy.table.pprint.MAX_LINES` configuration item. If a negative
        value of `max_lines` is supplied then there is no line limit applied.

        Parameters
        ----------
        max_lines : int
            Maximum number of values in output

        show_name : bool
            Include column name (default=True)

        show_unit : bool
            Include a header row for unit (default=False)
        """
        lines, n_header = _pformat_col(self, max_lines, show_name, show_unit)
        for i, line in enumerate(lines):
            if i < n_header:
                color_print(line, 'red')
            else:
                print line

    def more(self, max_lines=None, show_name=True, show_unit=False):
        """Interactively browse column with a paging interface.

        Supported keys::

          f, <space> : forward one page
          b : back one page
          r : refresh same page
          n : next row
          p : previous row
          < : go to beginning
          > : go to end
          q : quit browsing
          h : print this help

        Parameters
        ----------
        max_lines : int
            Maximum number of lines in table output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit (default=False)

        """
        _more_tabcol(self, max_lines=max_lines, show_name=show_name,
                     show_unit=show_unit)

    @property
    def unit(self):
        """
        The unit associated with this column.  May be a string or a
        `astropy.units.UnitBase` instance.

        Setting the `unit` property does not change the values of the
        data.  To perform a unit conversion, use `convert_unit_to`.
        """
        return self._unit

    @unit.setter
    def unit(self, unit):
        if unit is None:
            self._unit = None
        else:
            self._unit = Unit(unit, parse_strict='silent')

    @unit.deleter
    def unit(self):
        self._unit = None
    
    @property
    @deprecated('0.3', alternative=':attr:`Column.unit`')
    def units(self):
        return self.unit
        
    @units.setter
    @deprecated('0.3', alternative=':attr:`Column.unit`')
    def units(self, unit):
        self.unit = unit
    
    @units.deleter
    @deprecated('0.3', alternative=':attr:`Column.unit`')
    def units(self):
        del self.unit
    
    def convert_unit_to(self, new_unit, equivalencies=[]):
        """
        Converts the values of the column in-place from the current
        unit to the given unit.

        To change the unit associated with this column without
        actually changing the data values, simply set the `unit`
        property.

        Parameters
        ----------
        new_unit : str or `astropy.units.UnitBase` instance
            The unit to convert to.

        equivalencies : list of equivalence pairs, optional
           A list of equivalence pairs to try if the unit are not
           directly convertible.  See :ref:`unit_equivalencies`.

        Raises
        ------
        astropy.units.UnitException
            If units are inconsistent
        """
        if self.unit is None:
            raise ValueError("No unit set on column")
        self.data[:] = self.unit.to(
            new_unit, self.data, equivalencies=equivalencies)
        self.unit = new_unit

    def __str__(self):
        lines, n_header = _pformat_col(self)
        return '\n'.join(lines)


class Column(BaseColumn, np.ndarray):
    """Define a data column for use in a Table object.

    Parameters
    ----------
    data : list, ndarray or None
        Column data values
    name : str
        Column name and key for reference within Table
    dtype : numpy.dtype compatible value
        Data type for column
    shape : tuple or ()
        Dimensions of a single row element in the column data
    length : int or 0
        Number of row elements in column data
    description : str or None
        Full description of column
    unit : str or None
        Physical unit
    format : str or None or function
        Format string for outputting column values.  This can be an
        "old-style" (``format % value``) or "new-style" (`str.format`)
        format specification string or a function that accepts a single
        value and returns a string.
    meta : dict-like or None
        Meta-data associated with the column

    Examples
    --------
    A Column can be created in two different ways:

    - Provide a ``data`` value but not ``shape`` or ``length`` (which are
      inferred from the data).

      Examples::

        col = Column(data=[1, 2], name='name')  # shape=(2,)
        col = Column(data=[[1, 2], [3, 4]], name='name')  # shape=(2, 2)
        col = Column(data=[1, 2], name='name', dtype=float)
        col = Column(data=np.array([1, 2]), name='name')
        col = Column(data=['hello', 'world'], name='name')

      The ``dtype`` argument can be any value which is an acceptable
      fixed-size data-type initializer for the numpy.dtype() method.  See
      `<http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html>`_.
      Examples include:

      - Python non-string type (float, int, bool)
      - Numpy non-string type (e.g. np.float32, np.int64, np.bool)
      - Numpy.dtype array-protocol type strings (e.g. 'i4', 'f8', 'S15')

      If no ``dtype`` value is provide then the type is inferred using
      ``np.array(data)``.

    - Provide ``length`` and optionally ``shape``, but not ``data``

      Examples::

        col = Column(name='name', length=5)
        col = Column(name='name', dtype=int, length=10, shape=(3,4))

      The default ``dtype`` is ``np.float64``.  The ``shape`` argument is the
      array shape of a single cell in the column.

    .. warning::

       In the next major release of `astropy` (0.3), the order of function
       arguments for creating a |Column| will change.  Currently the order is
       ``Column(name, data, ...)``, but in 0.3 and later it will be
       ``Column(data, name, ...)``.  This improves consistency with |Table| and
       `numpy`.

       In order to use the same code for Astropy 0.2 and 0.3, column objects
       should always be created using named keyword arguments for ``data`` and
       ``name``, for instance ``c = Column(data=[1, 2], name='col')``.  When
       Astropy 0.3 is released then the the keyword identifiers can be dropped,
       allowing for ``c = Column([1, 2], 'c')``.
    """

    @_check_column_new_args
    def __new__(cls, data=None, name=None,
                dtype=None, shape=(), length=0,
                description=None, unit=None, format=None, meta=None,
                dtypes=None, units=None):
        
        if dtypes is not None:
            dtype = dtypes
            warnings.warn("'dtypes' has been renamed to the singular 'dtype'.",
                          DeprecationWarning)
        
        if units is not None:
            unit = units
            warnings.warn("'units' has been renamed to the singular 'unit'.",
                          DeprecationWarning)
            
        if data is None:
            dtype = (np.dtype(dtype).str, shape)
            self_data = np.zeros(length, dtype=dtype)
        elif isinstance(data, Column):
            self_data = np.asarray(data.data, dtype=dtype)
            if description is None:
                description = data.description
            if unit is None:
                unit = unit or data.unit
            if format is None:
                format = data.format
            if meta is None:
                meta = deepcopy(data.meta)
        elif isinstance(data, MaskedColumn):
            raise TypeError("Cannot convert a MaskedColumn to a Column")
        else:
            self_data = np.asarray(data, dtype=dtype)

        self = self_data.view(cls)
        self._name = name
        self.unit = unit
        self.format = format
        self.description = description
        self.parent_table = None

        self.meta = OrderedDict()
        if meta is not None:
            self.meta.update(meta)

        return self

    @property
    def data(self):
        return self.view(np.ndarray)

    def copy(self, data=None, copy_data=True):
        """Return a copy of the current Column instance.  If ``data`` is supplied
        then a view (reference) of ``data`` is used, and ``copy_data`` is ignored.
        """
        if data is None:
            data = self.view(np.ndarray)
            if copy_data:
                data = data.copy()

        return Column(name=self.name, data=data, unit=self.unit, format=self.format,
                      description=self.description, meta=deepcopy(self.meta))


class MaskedColumn(BaseColumn, ma.MaskedArray):
    """Define a masked data column for use in a Table object.

    Parameters
    ----------
    data : list, ndarray or None
        Column data values
    name : str
        Column name and key for reference within Table
    mask : list, ndarray or None
        Boolean mask for which True indicates missing or invalid data
    fill_value : float, int, str or None
        Value used when filling masked column elements
    dtype : numpy.dtype compatible value
        Data type for column
    shape : tuple or ()
        Dimensions of a single row element in the column data
    length : int or 0
        Number of row elements in column data
    description : str or None
        Full description of column
    unit : str or None
        Physical unit
    format : str or None or function
        Format string for outputting column values.  This can be an
        "old-style" (``format % value``) or "new-style" (`str.format`)
        format specification string or a function that accepts a single
        value and returns a string.
    meta : dict-like or None
        Meta-data associated with the column

    Examples
    --------
    A MaskedColumn is similar to a Column except that it includes ``mask`` and
    ``fill_value`` attributes.  It can be created in two different ways:

    - Provide a ``data`` value but not ``shape`` or ``length`` (which are
      inferred from the data).

      Examples::

        col = MaskedColumn(data=[1, 2], name='name')
        col = MaskedColumn(data=[1, 2], name='name', mask=[True, False])
        col = MaskedColumn(data=[1, 2], name='name', dtype=float, fill_value=99)

      The ``mask`` argument will be cast as a boolean array and specifies
      which elements are considered to be missing or invalid.

      The ``dtype`` argument can be any value which is an acceptable
      fixed-size data-type initializer for the numpy.dtype() method.  See
      `<http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html>`_.
      Examples include:

      - Python non-string type (float, int, bool)
      - Numpy non-string type (e.g. np.float32, np.int64, np.bool)
      - Numpy.dtype array-protocol type strings (e.g. 'i4', 'f8', 'S15')

      If no ``dtype`` value is provide then the type is inferred using
      ``np.array(data)``.  When ``data`` is provided then the ``shape``
      and ``length`` arguments are ignored.

    - Provide ``length`` and optionally ``shape``, but not ``data``

      Examples::

        col = MaskedColumn(name='name', length=5)
        col = MaskedColumn(name='name', dtype=int, length=10, shape=(3,4))

      The default ``dtype`` is ``np.float64``.  The ``shape`` argument is the
      array shape of a single cell in the column.

    .. warning::

       In the next major release of `astropy` (0.3), the order of function
       arguments for creating a |MaskedColumn| will change.  Currently the order is
       ``MaskedColumn(name, data, ...)``, but in 0.3 and later it will be
       ``MaskedColumn(data, name, ...)``.  This improves consistency with |Table|
       and `numpy`.

       In order to use the same code for Astropy 0.2 and 0.3, column objects
       should always be created using named keyword arguments for ``data`` and
       ``name``, for instance ``c = MaskedColumn(data=[1, 2], name='col')``.  When
       Astropy 0.3 is released then the the keyword identifiers can be dropped,
       allowing for ``c = MaskedColumn([1, 2], 'c')``.
    """

    @_check_column_new_args
    def __new__(cls, data=None, name=None, mask=None, fill_value=None,
                dtype=None, shape=(), length=0,
                description=None, unit=None, format=None, meta=None,
                units=None, dtypes=None):
        
        if dtypes is not None:
            dtype = dtypes
            warnings.warn("'dtypes' has been renamed to the singular 'dtype'.",
                          DeprecationWarning)
        
        if units is not None:
            unit = units
            warnings.warn("'units' has been renamed to the singular 'unit'.",
                          DeprecationWarning)
        
        if NUMPY_LT_1P5:
            raise ValueError('MaskedColumn requires NumPy version 1.5 or later')

        if data is None:
            dtype = (np.dtype(dtype).str, shape)
            self_data = ma.zeros(length, dtype=dtype)
        elif isinstance(data, (Column, MaskedColumn)):
            self_data = ma.asarray(data.data, dtype=dtype)
            if description is None:
                description = data.description
            if unit is None:
                unit = unit or data.unit
            if format is None:
                format = data.format
            if meta is None:
                meta = deepcopy(data.meta)
        else:
            self_data = ma.asarray(data, dtype=dtype)

        self = self_data.view(MaskedColumn)
        if mask is None and hasattr(data, 'mask'):
            mask = data.mask
        if fill_value is None and hasattr(data, 'fill_value'):
            fill_value = data.fill_value
        self.mask = mask
        self.fill_value = fill_value
        self._name = name
        self.unit = unit
        self.format = format
        self.description = description
        self.parent_table = None

        self.meta = OrderedDict()
        if meta is not None:
            self.meta.update(meta)

        return self

    def __array_finalize__(self, obj):
        BaseColumn.__array_finalize__(self, obj)
        ma.MaskedArray.__array_finalize__(self, obj)

    def _fix_fill_value(self, val):
        """Fix a fill value (if needed) to work around a bug with setting the fill
        value of a string array in MaskedArray with Python 3.x.  See
        https://github.com/numpy/numpy/pull/2733.  This mimics the check in
        numpy.ma.core._check_fill_value() (version < 1.7) which incorrectly sets
        fill_value to a default if self.dtype.char is 'U' (which is the case for Python
        3).  Here we change the string to a byte string so that in Python 3 the
        isinstance(val, basestring) part fails.
        """
        if isinstance(val, basestring) and (self.dtype.char not in 'SV'):
            val = val.encode()
        return val

    @property
    def fill_value(self):
        return self.get_fill_value()  # defer to native ma.MaskedArray method

    @fill_value.setter
    def fill_value(self, val):
        """Set fill value both in the masked column view and in the parent table
        if it exists.  Setting one or the other alone doesn't work."""
        val = self._fix_fill_value(val)

        if self.parent_table:
            self.parent_table._data[self._name].fill_value = val

        # Yet another ma bug workaround: If the value of fill_value for a string array is
        # requested but not yet set then it gets created as 'N/A'.  From this point onward
        # any new fill_values are truncated to 3 characters.  Note that this does not
        # occur if the masked array is a structured array (as in the previous block that
        # deals with the parent table).
        #
        # >>> x = ma.array(['xxxx'])
        # >>> x.fill_value  # fill_value now gets represented as an 'S3' array
        # 'N/A'
        # >>> x.fill_value='yyyy'
        # >>> x.fill_value
        # 'yyy'
        #
        # To handle this we are forced to reset a private variable first:
        self._fill_value = None

        self.set_fill_value(val)  # defer to native ma.MaskedArray method

    @property
    def data(self):
        out = self.view(ma.MaskedArray)
        # The following is necessary because of a bug in Numpy, which was
        # fixed in numpy/numpy#2703. The fix should be included in Numpy 1.8.0.
        out.fill_value = self.fill_value
        return out

    def filled(self, fill_value=None):
        """Return a copy of self, with masked values filled with a given value.

        Parameters
        ----------
        fill_value : scalar; optional
            The value to use for invalid entries (None by default).
            If None, the `fill_value` attribute of the array is used instead.

        Returns
        -------
        filled_column : Column
            A copy of ``self`` with masked entries replaced by `fill_value`
            (be it the function argument or the attribute of ``self``).
        """
        if fill_value is None:
            fill_value = self.fill_value
        fill_value = self._fix_fill_value(fill_value)

        data = super(MaskedColumn, self).filled(fill_value)
        out = Column(name=self.name, data=data, unit=self.unit, format=self.format,
                     description=self.description, meta=deepcopy(self.meta))
        return out

    def copy(self, data=None, copy_data=True):
        """
        Return a copy of the current MaskedColumn instance.  If ``data`` is supplied
        then a view (reference) of ``data`` is used, and ``copy_data`` is ignored.

        Parameters
        ----------
        data : array; optional
            Data to use when creating MaskedColumn copy.  If not supplied the
            column data array is used.
        copy_data : boolean; optional
            Make a copy of input data instead of using a reference (default=True)

        Returns
        -------
        column : MaskedColumn
            A copy of ``self``
        """
        if data is None:
            data = self.view(ma.MaskedArray)
            if copy_data:
                data = data.copy()

        return MaskedColumn(name=self.name, data=data, unit=self.unit, format=self.format,
                            # Do not include mask=self.mask since `data` has the mask
                            fill_value=self.fill_value,
                            description=self.description, meta=deepcopy(self.meta))


class Row(object):
    """A class to represent one row of a Table object.

    A Row object is returned when a Table object is indexed with an integer
    or when iterating over a table::

      >>> table = Table([(1, 2), (3, 4)], names=('a', 'b'),
      ...               dtypes=('int32', 'int32'))
      >>> row = table[1]
      >>> row
      <Row 1 of table
       values=(2, 4)
       dtype=[('a', '<i4'), ('b', '<i4')]>
      >>> row['a']
      2
      >>> row[1]
      4
    """

    def __init__(self, table, index):
        self._table = table
        self._index = index
        self._data = table._data[index]

        # MaskedArray __getitem__ has a strange behavior where if a
        # row mask is all False then it returns a np.void which
        # has no mask attribute. This makes it impossible to then set
        # the mask. Here we recast back to mvoid. This was fixed in
        # Numpy following issue numpy/numpy#483, and the fix should be
        # included in Numpy 1.8.0.
        if self._table.masked and isinstance(self._data, np.void):
            self._data = ma.core.mvoid(self._data,
                                       mask=self._table._mask[index])

    def __getitem__(self, item):
        return self.data[item]

    def __setitem__(self, item, val):
        self.data[item] = val

    def __eq__(self, other):
        if self._table.masked:
            # Sent bug report to numpy-discussion group on 2012-Oct-21, subject:
            # "Comparing rows in a structured masked array raises exception"
            # No response, so this is still unresolved.
            raise ValueError('Unable to compare rows for masked table due to numpy.ma bug')
        return self.data == other

    def __ne__(self, other):
        if self._table.masked:
            raise ValueError('Unable to compare rows for masked table due to numpy.ma bug')
        return self.data != other

    @property
    def _mask(self):
        return self._data.mask

    def __array__(self, dtype=None):
        """Support converting Row to np.array via np.array(table).

        Coercion to a different dtype via np.array(table, dtype) is not
        supported and will raise a ValueError.
        """
        if dtype is not None:
            raise ValueError('Datatype coercion is not allowed')

        return np.array(self._data)

    def __len__(self):
        return len(self._data.dtype)

    @property
    def table(self):
        return self._table

    @property
    def index(self):
        return self._index

    @property
    def data(self):
        return self._data

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
    
    @property
    @deprecated('0.3', alternative=':attr:`Row.dtype`', pending=False)
    def dtypes(self):
        return self.dtype

    def __repr__(self):
        return "<Row {0} of table\n values={1!r}\n dtype={2}>".format(
            self.index, self.data, self.dtype)


collections.Sequence.register(Row)


class Table(object):
    """A class to represent tables of heterogeneous data.

    `Table` provides a class for heterogeneous tabular data, making use of a
    `numpy` structured array internally to store the data values.  A key
    enhancement provided by the `Table` class is the ability to easily modify
    the structure of the table by adding or removing columns, or adding new
    rows of data.  In addition table and column metadata are fully supported.

    `Table` differs from `NDData` by the assumption that the input data
    consists of columns of homogeneous data, where each column has a unique
    identifier and may contain additional metadata such as the data
    unit, format, and description.

    Parameters
    ----------
    data : numpy ndarray, dict, list, or Table, optional
        Data to initialize table.
    masked : boolean, optional
        Specify whether the table is masked.
    names : list, optional
        Specify column names
    dtype : list, optional
        Specify column data types
    meta : dict, optional
        Metadata associated with the table.
    copy : boolean, optional
        Copy the input data (default=True).

    """

    def __init__(self, data=None, masked=None, names=None, dtype=None,
                 meta=None, copy=True, dtypes=None):
        
        if dtypes is not None:
            dtype = dtypes
            warnings.warn("'dtypes' has been renamed to the singular 'dtype'.",
                          DeprecationWarning)
        
        # Set up a placeholder empty table
        self._data = None
        self._set_masked(masked)
        self.columns = TableColumns()
        self._meta = OrderedDict() if meta is None else deepcopy(meta)

        # Must copy if dtype are changing
        if not copy and dtype is not None:
            raise ValueError('Cannot specify dtype when copy=False')

        # Infer the type of the input data and set up the initialization
        # function, number of columns, and potentially the default col names

        default_names = None

        if isinstance(data, (list, tuple)):
            init_func = self._init_from_list
            if data and all(isinstance(row, dict) for row in data):
                n_cols = len(data[0])
            else:
                n_cols = len(data)

        elif isinstance(data, np.ndarray):
            if data.dtype.names:
                init_func = self._init_from_ndarray  # _struct
                n_cols = len(data.dtype.names)
                default_names = data.dtype.names
            else:
                init_func = self._init_from_ndarray  # _homog
                n_cols = data.shape[1]

        elif isinstance(data, dict):
            init_func = self._init_from_dict
            n_cols = len(data.keys())
            default_names = data.keys()

        elif isinstance(data, Table):
            init_func = self._init_from_table
            n_cols = len(data.colnames)
            default_names = data.colnames

        elif data is None:
            if names is None:
                return  # Empty table
            else:
                init_func = self._init_from_list
                n_cols = len(names)
                data = [[]] * n_cols
        else:
            raise ValueError('Data type {0} not allowed to init Table'
                             .format(type(data)))

        # Set up defaults if names and/or dtype are not specified.
        # A value of None means the actual value will be inferred
        # within the appropriate initialization routine, either from
        # existing specification or auto-generated.

        if names is None:
            names = default_names or [None] * n_cols
        if dtype is None:
            dtype = [None] * n_cols
        self._check_names_dtype(names, dtype, n_cols)

        # Finally do the real initialization
        init_func(data, names, dtype, n_cols, copy)

        # Whatever happens above, the masked property should be set to a boolean
        if type(self.masked) != bool:
            raise TypeError("masked property has not been set to True or False")

        if NUMPY_LT_1P5 and self.masked:
            raise ValueError('Masked table requires NumPy version 1.5 or later')

    @property
    def mask(self):
        return self._data.mask if self.masked else None

    @property
    def _mask(self):
        """This is needed due to intricacies in numpy.ma, don't remove it."""
        return self._data.mask

    def filled(self, fill_value=None):
        """Return a copy of self, with masked values filled.

        If input ``fill_value`` supplied then that value is used for all masked entries
        in the table.  Otherwise the individual ``fill_value`` defined for each
        table column is used.

        Returns
        -------
        filled_table : Table
            New table with masked values filled
        """
        if self.masked:
            data = [col.filled(fill_value) for col in self.columns.values()]
        else:
            data = self
        return self.__class__(data, meta=deepcopy(self.meta))

    def __array__(self, dtype=None):
        """Support converting Table to np.array via np.array(table).

        Coercion to a different dtype via np.array(table, dtype) is not
        supported and will raise a ValueError.
        """
        if dtype is not None:
            raise ValueError('Datatype coercion is not allowed')

        # This limitation is because of the following unexpected result that
        # should have made a table copy while changing the column names.
        #
        # >>> d = astropy.table.Table([[1,2],[3,4]])
        # >>> np.array(d, dtype=[('a', 'i8'), ('b', 'i8')])
        # array([(0, 0), (0, 0)],
        #       dtype=[('a', '<i8'), ('b', '<i8')])

        return self._data.data if self.masked else self._data

    def _rebuild_table_column_views(self):
        """
        Some table manipulations can corrupt the Column views of self._data.  This
        function will cleanly rebuild the columns and self.columns.  This is a slightly
        subtle operation, see comments.
        """
        cols = []
        for col in self.columns.values():
            # First make a new column based on the name and the original column.  This
            # step is needed because the table manipulation may have changed the table
            # masking so that the original data columns no longer correspond to
            # self.ColumnClass.  This uses data refs, not copies.
            newcol = self.ColumnClass(name=col.name, data=col)

            # Now use the copy() method to copy the column and its metadata, but at
            # the same time set the column data to a view of self._data[col.name].
            # Somewhat confusingly in this case copy() refers to copying the
            # column attributes, but the data are used by reference.
            newcol = newcol.copy(data=self._data[col.name])
            cols.append(newcol)

        self.columns = TableColumns(cols)

    def _check_names_dtype(self, names, dtype, n_cols):
        """Make sure that names and dtype are boths iterable and have
        the same length as data.
        """
        for inp_list, inp_str in ((dtype, 'dtype'), (names, 'names')):
            if not isiterable(inp_list):
                raise ValueError('{0} must be a list or None'.format(inp_str))

        if len(names) != n_cols or len(dtype) != n_cols:
            raise ValueError(
                'Arguments "names" and "dtype" must match number of columns'
                .format(inp_str))

    def _set_masked_from_cols(self, cols):
        if self.masked is None:
            if any(isinstance(col, (MaskedColumn, ma.MaskedArray)) for col in cols):
                self._set_masked(True)
            else:
                self._set_masked(False)
        elif not self.masked:
            if any(isinstance(col, (MaskedColumn, ma.MaskedArray)) for col in cols):
                self._set_masked(True)

    def _init_from_list(self, data, names, dtype, n_cols, copy):
        """Initialize table from a list of columns.  A column can be a
        Column object, np.ndarray, or any other iterable object.
        """
        if not copy:
            raise ValueError('Cannot use copy=False with a list data input')

        # Set self.masked appropriately, then get class to create column instances.
        self._set_masked_from_cols(data)

        cols = []
        def_names = _auto_names(n_cols)

        if data and all(isinstance(row, dict) for row in data):
            names_from_data = set()
            for row in data:
                names_from_data.update(row)

            cols = {}
            for name in names_from_data:
                cols[name] = []
                for i, row in enumerate(data):
                    try:
                        cols[name].append(row[name])
                    except KeyError:
                        raise ValueError('Row {0} has no value for column {1}'.format(i, name))
            if all(name is None for name in names):
                names = sorted(names_from_data)
            self._init_from_dict(cols, names, dtype, n_cols, copy)
            return

        for col, name, def_name, dtype in zip(data, names, def_names, dtype):
            if isinstance(col, (Column, MaskedColumn)):
                col = self.ColumnClass(name=(name or col.name), data=col, dtype=dtype)
            elif isinstance(col, np.ndarray) or isiterable(col):
                col = self.ColumnClass(name=(name or def_name), data=col, dtype=dtype)
            else:
                raise ValueError('Elements in list initialization must be '
                                 'either Column or list-like')
            cols.append(col)

        self._init_from_cols(cols)

    def _init_from_ndarray(self, data, names, dtype, n_cols, copy):
        """Initialize table from an ndarray structured array"""

        data_names = data.dtype.names or _auto_names(n_cols)
        struct = data.dtype.names is not None
        names = [name or data_names[i] for i, name in enumerate(names)]

        cols = ([data[name] for name in data_names] if struct else
                [data[:, i] for i in range(n_cols)])

        # Set self.masked appropriately, then get class to create column instances.
        self._set_masked_from_cols(cols)

        if copy:
            self._init_from_list(cols, names, dtype, n_cols, copy)
        else:
            dtype = [(name, col.dtype) for name, col in zip(names, cols)]
            self._data = data.view(dtype).ravel()
            columns = TableColumns()

            for name in names:
                columns[name] = self.ColumnClass(name=name, data=self._data[name])
                columns[name].parent_table = self
            self.columns = columns

    def _init_from_dict(self, data, names, dtype, n_cols, copy):
        """Initialize table from a dictionary of columns"""

        if not copy:
            raise ValueError('Cannot use copy=False with a dict data input')

        data_list = [data[name] for name in names]
        self._init_from_list(data_list, names, dtype, n_cols, copy)

    def _init_from_table(self, data, names, dtype, n_cols, copy):
        """Initialize table from an existing Table object """

        table = data  # data is really a Table, rename for clarity
        data_names = table.colnames
        self.meta.clear()
        self.meta.update(deepcopy(table.meta))
        cols = table.columns.values()

        # Set self.masked appropriately from cols
        self._set_masked_from_cols(cols)

        if copy:
            self._init_from_list(cols, names, dtype, n_cols, copy)
        else:
            names = [vals[0] or vals[1] for vals in zip(names, data_names)]
            dtype = [(name, col.dtype) for name, col in zip(names, cols)]
            data = table._data.view(dtype)

            self._update_table_from_cols(self, data, cols, names)

    def _init_from_cols(self, cols):
        """Initialize table from a list of Column objects"""

        lengths = set(len(col.data) for col in cols)
        if len(lengths) != 1:
            raise ValueError('Inconsistent data column lengths: {0}'
                             .format(lengths))

        self._set_masked_from_cols(cols)
        cols = [self.ColumnClass(name=col.name, data=col) for col in cols]

        names = [col.name for col in cols]
        dtype = [col.descr for col in cols]
        empty_init = ma.empty if self.masked else np.empty
        data = empty_init(lengths.pop(), dtype=dtype)
        for col in cols:
            data[col.name] = col.data

        self._update_table_from_cols(self, data, cols, names)

    def _new_from_slice(self, slice_):
        """Create a new table as a referenced slice from self."""

        table = self.__class__(masked=self.masked)
        table.meta.clear()
        table.meta.update(deepcopy(self.meta))
        cols = self.columns.values()
        names = [col.name for col in cols]
        data = self._data[slice_]

        self._update_table_from_cols(table, data, cols, names)

        return table

    @staticmethod
    def _update_table_from_cols(table, data, cols, names):
        """Update the existing ``table`` so that it represents the given
        ``data`` (a structured ndarray) with ``cols`` and ``names``."""

        columns = TableColumns()
        table._data = data

        for name, col in zip(names, cols):
            newcol = col.copy(data=data[name], copy_data=False)
            newcol.name = name
            newcol.parent_table = table
            columns[name] = newcol
        table.columns = columns

    def __repr__(self):
        names = ("'{0}'".format(x) for x in self.colnames)
        s = "<Table rows={0} names=({1})>\n{2}".format(
            self.__len__(), ','.join(names), repr(self._data))
        return s

    def __str__(self):
        lines, n_header = _pformat_table(self)
        return '\n'.join(lines)

    def pprint(self, max_lines=None, max_width=None, show_name=True,
               show_unit=False):
        """Print a formatted string representation of the table.

        If no value of `max_lines` is supplied then the height of the screen
        terminal is used to set `max_lines`.  If the terminal height cannot
        be determined then the default is taken from the configuration item
        `astropy.table.pprint.MAX_LINES`.  If a negative value of `max_lines`
        is supplied then there is no line limit applied.

        The same applies for max_width except the configuration item is
        `astropy.table.pprint.MAX_WIDTH`.

        Parameters
        ----------
        max_lines : int
            Maximum number of lines in table output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit (default=False)
        """

        lines, n_header = _pformat_table(self, max_lines, max_width, show_name,
                                         show_unit)
        for i, line in enumerate(lines):
            if i < n_header:
                color_print(line, 'red')
            else:
                print line

    def pformat(self, max_lines=None, max_width=None, show_name=True,
                show_unit=False, html=False):
        """Return a list of lines for the formatted string representation of
        the table.

        If no value of `max_lines` is supplied then the height of the screen
        terminal is used to set `max_lines`.  If the terminal height cannot
        be determined then the default is taken from the configuration item
        `astropy.table.pprint.MAX_LINES`.  If a negative value of `max_lines`
        is supplied then there is no line limit applied.

        The same applies for max_width except the configuration item  is
        `astropy.table.pprint.MAX_WIDTH`.

        Parameters
        ----------
        max_lines : int or None
            Maximum number of rows to output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit (default=False)

        html : bool
            Format the output as an HTML table (default=False)

        Returns
        -------
        lines : list
            Formatted table as a list of strings
        """
        lines, n_header = _pformat_table(self, max_lines, max_width,
                                         show_name, show_unit, html)
        return lines

    def more(self, max_lines=None, max_width=None, show_name=True,
             show_unit=False):
        """Interactively browse table with a paging interface.

        Supported keys::

          f, <space> : forward one page
          b : back one page
          r : refresh same page
          n : next row
          p : previous row
          < : go to beginning
          > : go to end
          q : quit browsing
          h : print this help

        Parameters
        ----------
        max_lines : int
            Maximum number of lines in table output

        max_width : int or None
            Maximum character width of output

        show_name : bool
            Include a header row for column names (default=True)

        show_unit : bool
            Include a header row for unit (default=False)
        """
        _more_tabcol(self, max_lines, max_width, show_name,
                     show_unit)

    def _repr_html_(self):
        lines = self.pformat(html=True)
        return ''.join(lines)

    def __getitem__(self, item):
        if isinstance(item, basestring):
            return self.columns[item]
        elif isinstance(item, int):
            return Row(self, item)
        elif isinstance(item, tuple):
            if all(isinstance(x, np.ndarray) for x in item):
                # Item is a tuple of ndarrays as in the output of np.where, e.g.
                # t[np.where(t['a'] > 2)]
                return self._new_from_slice(item)
            elif (all(x in self.colnames for x in item)):
                # Item is a tuple of strings that are valid column names
                return self.__class__([self[x] for x in item], meta=deepcopy(self.meta))
            else:
                raise ValueError('Illegal item for table item access')

        elif (isinstance(item, slice) or isinstance(item, np.ndarray)
              or isinstance(item, list)):
            return self._new_from_slice(item)
        else:
            raise ValueError('Illegal type {0} for table item access'
                             .format(type(item)))

    def __setitem__(self, item, value):
        # If the item is a string then it must be the name of a column.
        # If that column doesn't already exist then create it now.
        if isinstance(item, basestring) and item not in self.colnames:
            NewColumn = MaskedColumn if self.masked else Column

            # Make sure value is an ndarray so we can get the dtype
            if not isinstance(value, np.ndarray):
                value = np.asarray(value)

            # Make new column and assign the value.  If the table currently has no rows
            # (len=0) of the value is already a Column then define new column directly
            # from value.  In the latter case this allows for propagation of Column
            # metadata.  Otherwise define a new column with the right length and shape and
            # then set it from value.  This allows for broadcasting, e.g. t['a'] = 1.
            if isinstance(value, BaseColumn):
                new_column = value.copy(copy_data=False)
                new_column.name = item
            elif len(self) == 0:
                new_column = NewColumn(name=item, data=value)
            else:
                new_column = NewColumn(name=item, length=len(self), dtype=value.dtype,
                                       shape=value.shape[1:])
                new_column[:] = value

            # Now add new column to the table
            self.add_column(new_column)
        elif isinstance(value, Row):
            # Value is another row
            self._data[item] = value.data
        else:
            # Otherwise just delegate to the numpy item setter.
            self._data[item] = value

    def __delitem__(self, item):
        if isinstance(item, basestring):
            self.remove_column(item)
        elif isinstance(item, tuple):
            self.remove_columns(item)

    def __iter__(self):
        self._iter_index = 0
        return self

    def __next__(self):
        """Python 3 iterator"""
        if self._iter_index < len(self._data):
            val = self[self._iter_index]
            self._iter_index += 1
            return val
        else:
            raise StopIteration

    if sys.version_info[0] < 3:  # pragma: py2
        next = __next__

    def field(self, item):
        """Return column[item] for recarray compatibility."""
        return self.columns[item]

    @property
    def masked(self):
        return self._masked

    @masked.setter
    def masked(self, masked):
        raise Exception('Masked attribute is read-only (use t = Table(t, masked=True)'
                        ' to convert to a masked table)')

    def _set_masked(self, masked):
        """
        Set the table masked property.

        Parameters
        ----------
        masked : bool
            State of table masking (True or False)
        """
        if hasattr(self, '_masked'):
            # The only allowed change is from None to False or True, or False to True
            if self._masked is None and masked in [False, True]:
                self._masked = masked
            elif self._masked is False and masked is True:
                log.info("Upgrading Table to masked Table")
                self._masked = masked
            elif self._masked is masked:
                raise Exception("Masked attribute is already set to {0}".format(masked))
            else:
                raise Exception("Cannot change masked attribute to {0} once it is set to {1}"
                                .format(masked, self._masked))
        else:
            if masked in [True, False, None]:
                self._masked = masked
            else:
                raise ValueError("masked should be one of True, False, None")
        if self._masked:
            self._column_class = MaskedColumn
        else:
            self._column_class = Column

    @property
    def ColumnClass(self):
        if self._column_class is None:
            return Column
        else:
            return self._column_class

    @property
    def dtype(self):
        return self._data.dtype

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

    def create_mask(self):
        if isinstance(self._data, ma.MaskedArray):
            raise Exception("data array is already masked")
        else:
            self._data = ma.array(self._data)

    def index_column(self, name):
        """
        Return the positional index of column ``name``.

        Parameters
        ----------
        name : str
            column name

        Returns
        -------
        index : int
            Positional index of column ``name``.

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> print t
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Get index of column 'b' of the table::

            >>> t.index_column('b')
            1
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

        Examples
        --------
        Create a table with two columns 'a' and 'b'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> print t
             a   b
            --- ---
              1 0.1
              2 0.2
              3 0.3

        Create a third column 'c' and append it to the end of the table::

            >>> col_c = Column(name='c', data=['x', 'y', 'z'])
            >>> t.add_column(col_c)
            >>> print t
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Add column 'd' at position 1. Note that the column is inserted
        before the given index::

            >>> col_d = Column(name='d', data=['a', 'b', 'c'])
            >>> t.add_column(col_d, 1)
            >>> print t
             a   d   b   c
            --- --- --- ---
              1   a 0.1   x
              2   b 0.2   y
              3   c 0.3   z

        To add several columns use add_columns.
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

        Examples
        --------
        Create a table with two columns 'a' and 'b'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> print t
             a   b
            --- ---
              1 0.1
              2 0.2
              3 0.3

        Create column 'c' and 'd' and append them to the end of the table::

            >>> col_c = Column(name='c', data=['x', 'y', 'z'])
            >>> col_d = Column(name='d', data=['u', 'v', 'w'])
            >>> t.add_columns([col_c, col_d])
            >>> print t
             a   b   c   d
            --- --- --- ---
              1 0.1   x   u
              2 0.2   y   v
              3 0.3   z   w

        Add column 'c' at position 0 and column 'd' at position 1. Note that
        the columns are inserted before the given position::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3]], names=('a', 'b'))
            >>> col_c = Column(name='c', data=['x', 'y', 'z'])
            >>> col_d = Column(name='d', data=['u', 'v', 'w'])
            >>> t.add_columns([col_c, col_d], [0, 1])
            >>> print t
             c   a   d   b
            --- --- --- ---
              x   1   u 0.1
              y   2   v 0.2
              z   3   w 0.3
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
        Remove a column from the table.

        This can also be done with::

          del table[name]

        Parameters
        ----------
        name : str
            Name of column to remove

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> print t
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Remove column 'b' from the table::

            >>> t.remove_column('b')
            >>> print t
             a   c
            --- ---
              1   x
              2   y
              3   z

        To remove several columns at the same time use remove_columns.
        """

        self.remove_columns([name])

    def remove_columns(self, names):
        '''
        Remove several columns from the table

        Parameters
        ----------
        names : list
            A list containing the names of the columns to remove

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...     names=('a', 'b', 'c'))
            >>> print t
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Remove columns 'b' and 'c' from the table::

            >>> t.remove_columns(['b', 'c'])
            >>> print t
             a
            ---
              1
              2
              3

        Specifying only a single column also works. Remove column 'b' from the table:

            >>> t = Table([[1, 2, 3], [0.1, 0.2, 0.3], ['x', 'y', 'z']],
            ...     names=('a', 'b', 'c'))
            >>> t.remove_columns('b')
            >>> print t
             a   c
            --- ---
              1   x
              2   y
              3   z

        This gives the same as using remove_column.
        '''

        for name in names:
            if name not in self.columns:
                raise KeyError("Column {0} does not exist".format(name))

        for name in names:
            self.columns.pop(name)

        newdtype = [(name, self._data.dtype[name]) for name in self._data.dtype.names
                if name not in names]
        newdtype = np.dtype(newdtype)

        if newdtype:
            if self.masked:
                table = np.ma.empty(self._data.shape, dtype=newdtype)
            else:
                table = np.empty(self._data.shape, dtype=newdtype)

            for field in newdtype.fields:
                table[field] = self._data[field]
                if self.masked:
                    table[field].fill_value = self._data[field].fill_value
        else:
            table = None

        self._data = table


    def keep_columns(self, names):
        '''
        Keep only the columns specified (remove the others).

        Parameters
        ----------
        names : list
            A list containing the names of the columns to keep. All other
            columns will be removed.

        Examples
        --------
        Create a table with three columns 'a', 'b' and 'c'::

            >>> t = Table([[1, 2, 3],[0.1, 0.2, 0.3],['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> print t
             a   b   c
            --- --- ---
              1 0.1   x
              2 0.2   y
              3 0.3   z

        Specifying only a single column name keeps only this column.
        Keep only column 'a' of the table::

            >>> t.keep_columns('a')
            >>> print t
             a
            ---
              1
              2
              3

        Specifying a list of column names is keeps is also possible.
        Keep columns 'a' and 'c' of the table::

            >>> t = Table([[1, 2, 3],[0.1, 0.2, 0.3],['x', 'y', 'z']],
            ...           names=('a', 'b', 'c'))
            >>> t.keep_columns(['a', 'c'])
            >>> print t
             a   c
            --- ---
              1   x
              2   y
              3   z
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

    def add_row(self, vals=None, mask=None):
        """Add a new row to the end of the table.

        The ``vals`` argument can be:

        sequence (e.g. tuple or list)
            Column values in the same order as table columns.
        mapping (e.g. dict)
            Keys corresponding to column names.  Missing values will be
            filled with np.zeros for the column dtype.
        None
            All values filled with np.zeros for the column dtype.

        This method requires that the Table object "owns" the underlying array
        data.  In particular one cannot add a row to a Table that was
        initialized with copy=False from an existing array.

        The ``mask`` attribute should give (if desired) the mask for the
        values. The type of the mask should match that of the values, i.e. if
        ``vals`` is an iterable, then ``mask`` should also be an iterable
        with the same length, and if ``vals`` is a mapping, then ``mask``
        should be a dictionary.

        Parameters
        ----------
        vals : tuple, list, dict or None
            Use the specified values in the new row
        """

        def _is_mapping(obj):
            """Minimal checker for mapping (dict-like) interface for obj"""
            attrs = ('__getitem__', '__len__', '__iter__', 'keys', 'values', 'items')
            return all(hasattr(obj, attr) for attr in attrs)

        newlen = len(self._data) + 1

        if vals is None:
            vals = np.zeros(1, dtype=self._data.dtype)[0]

        if mask is not None and not self.masked:
            self._set_masked(True)

        if self.masked:
            if newlen == 1:
                self._data = ma.empty(1, dtype=self._data.dtype)
            else:
                self._data = ma.resize(self._data, (newlen,))
        else:
            self._data.resize((newlen,), refcheck=False)

        if _is_mapping(vals):

            if mask is not None and not _is_mapping(mask):
                raise TypeError("Mismatch between type of vals and mask")

            # Now check that the mask is specified for the same keys as the
            # values, otherwise things get really confusing.
            if mask is not None and set(vals.keys()) != set(mask.keys()):
                raise ValueError('keys in mask should match keys in vals')

            if self.masked:
                # We set the mask to True regardless of whether a mask value
                # is specified or not - that is, any cell where a new row
                # value is not specified should be treated as missing.
                self._data.mask[-1] = (True,) * len(self._data.dtype)

            # First we copy the values
            for name, val in vals.items():
                try:
                    self._data[name][-1] = val
                except IndexError:
                    raise ValueError("No column {0} in table".format(name))
                if mask:
                    self._data[name].mask[-1] = mask[name]

        elif isiterable(vals):

            if mask is not None and (not isiterable(mask) or _is_mapping(mask)):
                raise TypeError("Mismatch between type of vals and mask")

            if len(self.columns) != len(vals):
                raise ValueError('Mismatch between number of vals and columns')

            if not isinstance(vals, tuple):
                vals = tuple(vals)

            self._data[-1] = vals

            if mask is not None:

                if len(self.columns) != len(mask):
                    raise ValueError('Mismatch between number of masks and columns')

                if not isinstance(mask, tuple):
                    mask = tuple(mask)

                self._data.mask[-1] = mask

        else:
            raise TypeError('Vals must be an iterable or mapping or None')

        self._rebuild_table_column_views()

    def sort(self, keys):
        '''
        Sort the table according to one or more keys. This operates
        on the existing table and does not return a new table.

        Parameters
        ----------
        keys : str or list of str
            The key(s) to order the table by
        '''
        if type(keys) is not list:
            keys = [keys]
        self._data.sort(order=keys)
        self._rebuild_table_column_views()

    def reverse(self):
        '''
        Reverse the row order of table rows.  The table is reversed
        in place and there are no function arguments.
        '''
        self._data[:] = self._data[::-1].copy()
        self._rebuild_table_column_views()

    read = classmethod(io_registry.read)
    write = io_registry.write

    @property
    def meta(self):
        return self._meta

    def copy(self, copy_data=True):
        '''
        Return a copy of the table


        Parameters
        ----------
        copy_data : bool
            If True (the default), copy the underlying data array.
            Otherwise, use the same data array
        '''
        return self.__class__(self, copy=copy_data)

    def __deepcopy__(self, memo=None):
        return self.copy(True)

    def __copy__(self):
        return self.copy(False)

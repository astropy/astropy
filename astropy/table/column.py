# Licensed under a 3-clause BSD style license - see LICENSE.rst

import abc
import functools
import operator

import warnings

from copy import deepcopy

import numpy as np
from numpy import ma

from ..units import Unit, Quantity
from ..utils import deprecated
from ..utils.console import color_print
from ..utils.exceptions import AstropyDeprecationWarning
from ..utils.metadata import MetaData
from . import groups
from .pprint import (_pformat_col, _pformat_col_iter, _more_tabcol)

from ..config import ConfigurationItem

# Python 2 and 3 source compatibility
try:
    unicode
except NameError:
    unicode = basestring = str

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


def _column_compare(op):
    """
    Convenience function to return a function that properly does a
    comparison between a column object and something else.
    """
    def compare(self, other):
        # We have to define this to ensure that we always return boolean arrays
        # (otherwise in some cases, Column objects are returned).
        if isinstance(other, BaseColumn):
            other = other.data
        return op(self.data, other)
    return compare


class BaseColumn(object):

    __metaclass__ = abc.ABCMeta

    meta = MetaData()

    # Define comparison operators
    __eq__ = _column_compare(operator.eq)
    __ne__ = _column_compare(operator.ne)
    __lt__ = _column_compare(operator.lt)
    __le__ = _column_compare(operator.le)
    __gt__ = _column_compare(operator.gt)
    __ge__ = _column_compare(operator.ge)

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
        col : Column
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
        astropy.units.UnitsError
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

    @property
    def groups(self):
        if not hasattr(self, '_groups'):
            self._groups = groups.ColumnGroups(self)
        return self._groups

    def group_by(self, keys):
        """
        Group this column by the specified ``keys``

        This effectively splits the column into groups which correspond to unique values of
        the ``keys`` grouping object.  The output is a new `Column` or `MaskedColumn` which
        contains a copy of this column but sorted by row according to ``keys``.

        The ``keys`` input to `group_by` must be a numpy array with the same length as
        this column.

        Parameters
        ----------
        keys : numpy array
            Key grouping object

        Returns
        -------
        out : Column
            New column with groups attribute set accordingly
        """
        return groups.column_group_by(self, keys)

    def _copy_groups(self, out):
        """
        Copy current groups into a copy of self ``out``
        """
        if self.parent_table:
            if hasattr(self.parent_table, '_groups'):
                out._groups = groups.ColumnGroups(out, indices=self.parent_table._groups._indices)
        elif hasattr(self, '_groups'):
            out._groups = groups.ColumnGroups(out, indices=self._groups._indices)


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
                          AstropyDeprecationWarning)

        if units is not None:
            unit = units
            warnings.warn("'units' has been renamed to the singular 'unit'.",
                          AstropyDeprecationWarning)

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
            if name is None:
                name = data.name
        elif isinstance(data, MaskedColumn):
            raise TypeError("Cannot convert a MaskedColumn to a Column")
        elif isinstance(data, Quantity):
            if unit is None:
                self_data = np.asarray(data, dtype=dtype)
                unit = data.unit
            else:
                self_data = np.asarray(data.to(unit), dtype=dtype)
        else:
            self_data = np.asarray(data, dtype=dtype)

        self = self_data.view(cls)
        self._name = name
        self.unit = unit
        self.format = format
        self.description = description
        self.parent_table = None
        self.meta = meta

        return self

    @property
    def data(self):
        return self.view(np.ndarray)

    def copy(self, order='C', data=None, copy_data=True):
        """Return a copy of the current Column instance.  If ``data`` is supplied
        then a view (reference) of ``data`` is used, and ``copy_data`` is ignored.
        """
        if data is None:
            data = self.view(np.ndarray)
            if copy_data:
                data = data.copy(order)

        out = Column(name=self.name, data=data, unit=self.unit, format=self.format,
                     description=self.description, meta=deepcopy(self.meta))
        self._copy_groups(out)

        return out


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
                          AstropyDeprecationWarning)

        if units is not None:
            unit = units
            warnings.warn("'units' has been renamed to the singular 'unit'.",
                          AstropyDeprecationWarning)

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
            if name is None:
                name = data.name
        elif isinstance(data, Quantity):
            if unit is None:
                self_data = ma.asarray(data, dtype=dtype)
                unit = data.unit
            else:
                self_data = ma.asarray(data.to(unit), dtype=dtype)
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
        self.meta = meta

        return self

    def __array_finalize__(self, obj):
        BaseColumn.__array_finalize__(self, obj)
        ma.MaskedArray.__array_finalize__(self, obj)

    # Surprisingly, MaskedArray is buggy and does not call __array_finalize__
    # after __getitem__, but instead does call _update_from, so we override
    # this instead to copy over the column metadata.
    def _update_from(self, obj):
        BaseColumn.__array_finalize__(self, obj)
        ma.MaskedArray._update_from(self, obj)

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

    def copy(self, order='C', data=None, copy_data=True):
        """
        Return a copy of the current MaskedColumn instance.  If ``data`` is supplied
        then a view (reference) of ``data`` is used, and ``copy_data`` is ignored.

        Parameters
        ----------
        data : array; optional
            Data to use when creating MaskedColumn copy.  If not supplied the
            column data array is used.
        copy_data : bool; optional
            Make a copy of input data instead of using a reference (default=True)

        Returns
        -------
        column : MaskedColumn
            A copy of ``self``
        """
        if data is None:
            data = self.view(ma.MaskedArray)
            if copy_data:
                data = data.copy(order)

        out = MaskedColumn(name=self.name, data=data, unit=self.unit, format=self.format,
                           # Do not include mask=self.mask since `data` has the mask
                           fill_value=self.fill_value,
                           description=self.description, meta=deepcopy(self.meta))
        self._copy_groups(out)
        return out

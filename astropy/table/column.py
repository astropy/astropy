# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six

import weakref

from copy import deepcopy

import numpy as np
from numpy import ma

from ..units import Unit, Quantity
from ..utils.compat import NUMPY_LT_1_8
from ..utils.console import color_print
from ..utils.metadata import MetaData
from . import groups
from . import pprint
from .np_utils import fix_column_name

from ..config import ConfigAlias


AUTO_COLNAME = ConfigAlias(
    '0.4', 'AUTO_COLNAME', 'auto_colname',
    'astropy.table.column', 'astropy.table')

# Create a generic TableFormatter object for use by bare columns with no
# parent table.
FORMATTER = pprint.TableFormatter()
INTEGER_TYPES = (int, long, np.integer) if six.PY2 else (int, np.integer)

def _auto_names(n_cols):
    from . import conf
    return [str(conf.auto_colname).format(i) for i in range(n_cols)]


# list of one and two-dimensional comparison functions, which sometimes return
# a Column class and sometimes a plain array. Used in __array_wrap__ to ensure
# they only return plain (masked) arrays (see #1446 and #1685)
_comparison_functions = set(
    [np.greater, np.greater_equal, np.less, np.less_equal,
     np.not_equal, np.equal,
     np.isfinite, np.isinf, np.isnan, np.sign, np.signbit])


COLUMN_ATTRS = set(['name', 'unit', 'dtype', 'format', 'description', 'meta', 'parent_table'])

def col_setattr(col, attr, value):
    """
    Set one of the column attributes.

    Warning: this function is subject to change or removal.
    """
    if attr not in COLUMN_ATTRS:
        raise AttributeError("attribute must be one of {0}".format(COLUMN_ATTRS))

    # The unit and dtype attributes are considered universal and do NOT get
    # stored in _astropy_column_attrs.  For BaseColumn instances use the usual setattr.
    if isinstance(col, BaseColumn):
        setattr(col, attr, value)
    else:
        # If no _astropy_column_attrs or it is None then convert to dict
        if getattr(col, '_astropy_column_attrs', None) is None:
            col._astropy_column_attrs = {}
        if attr == 'parent_table':
            value = None if value is None else weakref.ref(value)
        col._astropy_column_attrs[attr] = value

def col_getattr(col, attr, default=None):
    """
    Get one of the column attributes

    Warning: this function is subject to change or removal.
    """
    if attr not in COLUMN_ATTRS:
        raise AttributeError("attribute must be one of {0}".format(COLUMN_ATTRS))

    # The unit and dtype attributes are considered universal and do NOT get
    # stored in _astropy_column_attrs.  For BaseColumn instances use the usual setattr.
    if (isinstance(col, BaseColumn) or
            (isinstance(col, Quantity) and attr in ('dtype', 'unit'))):
        value = getattr(col, attr, default)
    else:
        # If col does not have _astropy_column_attrs or it is None (meaning
        # nothing has been set yet) then return default, otherwise look for
        # the attribute in the astropy_column_attrs dict.
        if getattr(col, '_astropy_column_attrs', None) is None:
            value = default
        else:
            value = col._astropy_column_attrs.get(attr, default)

        # Weak ref for parent table
        if attr == 'parent_table' and callable(value):
            value = value()

        # Mixins have a default dtype of Object if nothing else was set
        if attr == 'dtype' and value is None:
            value = np.dtype('O')

    return value

def _col_update_attrs_from(newcol, col, exclude_attrs=['name', 'parent_table']):
    """
    Update _astropy_column_attrs from mixin `col` to `newcol`.  Does nothing
    for BaseColumn cols

    Warning: this function is subject to change or removal.
    """
    if isinstance(newcol, BaseColumn):
        return

    attrs = COLUMN_ATTRS - set(exclude_attrs)
    for attr in attrs:
        val = col_getattr(col, attr)
        if val is not None:
            col_setattr(newcol, attr, deepcopy(val))

def col_iter_str_vals(col):
    """
    This is a mixin-safe version of Column.iter_str_vals.

    Warning: this function is subject to change or removal.
    """
    parent_table = col_getattr(col, 'parent_table')
    formatter = FORMATTER if parent_table is None else parent_table.formatter
    _pformat_col_iter = formatter._pformat_col_iter
    for str_val in _pformat_col_iter(col, -1, False, False, {}):
        yield str_val

def col_copy(col):
    """
    This is a mixin-safe version of Column.copy() (with copy_data=True).

    Warning: this function is subject to change or removal.
    """
    if isinstance(col, BaseColumn):
        return col.copy()

    if hasattr(col, '_astropy_column_attrs'):
        col_setattr(col, 'parent_table', None)  # Don't copy weakref to parent table
    newcol = col.copy() if hasattr(col, 'copy') else deepcopy(col)

    # Copy old attributes.  Even deepcopy above may not get this (e.g. pandas).
    if (not hasattr(newcol, '_astropy_column_attrs') or
            newcol._astropy_column_attrs is None):
        _column_attrs = deepcopy(getattr(col, '_astropy_column_attrs', {}))
        newcol._astropy_column_attrs = _column_attrs

    return newcol


class FalseArray(np.ndarray):
    def __new__(cls, shape):
        obj = np.zeros(shape, dtype=np.bool).view(cls)
        return obj

    def __setitem__(self, item, val):
        val = np.asarray(val)
        if np.any(val):
            raise ValueError('Cannot set any element of {0} class to True'
                             .format(self.__class__.__name__))

    def __setslice__(self, start, stop, val):
        val = np.asarray(val)
        if np.any(val):
            raise ValueError('Cannot set any element of {0} class to True'
                             .format(self.__class__.__name__))


class BaseColumn(np.ndarray):

    meta = MetaData()

    def __new__(cls, data=None, name=None,
                dtype=None, shape=(), length=0,
                description=None, unit=None, format=None, meta=None, copy=False):

        if data is None:
            dtype = (np.dtype(dtype).str, shape)
            self_data = np.zeros(length, dtype=dtype)
        elif isinstance(data, BaseColumn) and hasattr(data, '_name'):
            # When unpickling a MaskedColumn, ``data`` will be a bare
            # BaseColumn with none of the expected attributes.  In this case
            # do NOT execute this block which initializes from ``data``
            # attributes.
            self_data = np.array(data.data, dtype=dtype, copy=copy)
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
                self_data = np.array(data, dtype=dtype, copy=copy)
                unit = data.unit
            else:
                self_data = np.array(data.to(unit), dtype=dtype, copy=copy)
            if description is None:
                description = col_getattr(data, 'description')
            if format is None:
                format = col_getattr(data, 'format')
            if meta is None:
                meta = deepcopy(col_getattr(data, 'meta'))

        else:
            self_data = np.array(data, dtype=dtype, copy=copy)

        self = self_data.view(cls)
        self._name = fix_column_name(name)
        self.unit = unit
        self.format = format
        self.description = description
        self.meta = meta
        self._parent_table = None

        return self

    @property
    def data(self):
        return self.view(np.ndarray)

    @property
    def parent_table(self):
        if self._parent_table is None:
            return None
        else:
            return self._parent_table()

    @parent_table.setter
    def parent_table(self, table):
        if table is None:
            self._parent_table = None
        else:
            self._parent_table = weakref.ref(table)


    def copy(self, order='C', data=None, copy_data=True):
        """
        Return a copy of the current instance.

        If ``data`` is supplied then a view (reference) of ``data`` is used,
        and ``copy_data`` is ignored.

        Parameters
        ----------
        order : {'C', 'F', 'A', 'K'}, optional
            Controls the memory layout of the copy. 'C' means C-order,
            'F' means F-order, 'A' means 'F' if ``a`` is Fortran contiguous,
            'C' otherwise. 'K' means match the layout of ``a`` as closely
            as possible. (Note that this function and :func:numpy.copy are very
            similar, but have different default values for their order=
            arguments.)  Default is 'C'.
        data : array, optional
            If supplied then use a view of ``data`` instead of the instance
            data.  This allows copying the instance attributes and meta.
        copy_data : bool, optional
            Make a copy of the internal numpy array instead of using a
            reference.  Default is True.

        Returns
        -------
        col : Column or MaskedColumn
            Copy of the current column (same type as original)
        """
        if data is None:
            data = self.data
            if copy_data:
                data = data.copy(order)

        out = data.view(self.__class__)
        out.__array_finalize__(self)
        # for MaskedColumn, MaskedArray.__array_finalize__ also copies mask
        # from self, which is not the idea here, so undo
        if isinstance(self, MaskedColumn):
            out._mask = data._mask

        self._copy_groups(out)

        return out

    def __setstate__(self, state):
        """
        Restore the internal state of the Column/MaskedColumn for pickling
        purposes.  This requires that the last element of ``state`` is a
        5-tuple that has Column-specific state values.
        """
        # Get the Column attributes and meta
        name, unit, format, description, meta = state[-1]
        state = state[:-1]

        # Using super(type(self), self).__setstate__() gives an infinite
        # recursion.  Manually call the right super class to actually set up
        # the array object.
        super_class = ma.MaskedArray if isinstance(self, ma.MaskedArray) else np.ndarray
        super_class.__setstate__(self, state)

        # Set the Column attributes and meta
        self._name = name
        self.unit = unit
        self.format = format
        self.description = description
        self.meta = meta

    def __reduce__(self):
        """
        Return a 3-tuple for pickling a Column.  Use the super-class
        functionality but then add in a 5-tuple of Column-specific values
        that get used in __setstate__.
        """
        super_class = ma.MaskedArray if isinstance(self, ma.MaskedArray) else np.ndarray
        reconstruct_func, reconstruct_func_args, state = super_class.__reduce__(self)

        # Define Column-specific attrs and meta that gets added to state.
        column_state = (self.name, self.unit, self.format, self.description,
                        self.meta)
        state = state + (column_state,)

        return reconstruct_func, reconstruct_func_args, state

    def __getitem__(self, item):
        if isinstance(item, INTEGER_TYPES):
            return self.data[item]  # Return as plain ndarray or ma.MaskedArray
        else:
            return super(BaseColumn, self).__getitem__(item)

    # avoid == and != to be done based on type of subclass
    # (helped solve #1446; see also __array_wrap__)
    def __eq__(self, other):
        return self.data.__eq__(other)

    def __ne__(self, other):
        return self.data.__ne__(other)

    def __array_finalize__(self, obj):
        # Obj will be none for direct call to Column() creator
        if obj is None:
            return

        if six.callable(super(BaseColumn, self).__array_finalize__):
            super(BaseColumn, self).__array_finalize__(obj)

        # Self was created from template (e.g. obj[slice] or (obj * 2))
        # or viewcast e.g. obj.view(Column).  In either case we want to
        # init Column attributes for self from obj if possible.
        self.parent_table = None
        self._copy_attrs(obj)

    def __array_wrap__(self, out_arr, context=None):
        """
        __array_wrap__ is called at the end of every ufunc.

        Normally, we want a Column object back and do not have to do anything
        special. But there are two exceptions:

        1) If the output shape is different (e.g. for reduction ufuncs
           like sum() or mean()), a Column still linking to a parent_table
           makes little sense, so we return the output viewed as the
           column content (ndarray or MaskedArray).
           For this case, we use "[()]" to select everything, and to ensure we
           convert a zero rank array to a scalar. (For some reason np.sum()
           returns a zero rank scalar array while np.mean() returns a scalar;
           So the [()] is needed for this case.

        2) When the output is created by any function that returns a boolean
           we also want to consistently return an array rather than a column
           (see #1446 and #1685)
        """
        out_arr = super(BaseColumn, self).__array_wrap__(out_arr, context)
        if (self.shape != out_arr.shape or
            (isinstance(out_arr, BaseColumn) and
             (context is not None and context[0] in _comparison_functions))):
            return out_arr.data[()]
        else:
            return out_arr

    @property
    def name(self):
        """
        The name of this column.
        """
        return self._name

    @name.setter
    def name(self, val):
        val = fix_column_name(val)

        if self.parent_table is not None:
            table = self.parent_table
            table.columns._rename_column(self.name, val)

        self._name = val

    @property
    def descr(self):
        """Array-interface compliant full description of the column.

        This returns a 3-tuple (name, type, shape) that can always be
        used in a structured array dtype definition.
        """
        return (self.name, self.dtype.str, self.shape[1:])

    def iter_str_vals(self):
        """
        Return an iterator that yields the string-formatted values of this
        column.

        Returns
        -------
        str_vals : iterator
            Column values formatted as strings
        """
        # Iterate over formatted values with no max number of lines, no column
        # name, no unit, and ignoring the returned header info in outs.
        _pformat_col_iter = self._formatter._pformat_col_iter
        for str_val in _pformat_col_iter(self, -1, show_name=False, show_unit=False,
                                         show_dtype=False, outs={}):
            yield str_val

    def attrs_equal(self, col):
        """Compare the column attributes of ``col`` to this object.

        The comparison attributes are: ``name``, ``unit``, ``dtype``,
        ``format``, ``description``, and ``meta``.

        Parameters
        ----------
        col : Column
            Comparison column

        Returns
        -------
        equal : boolean
            True if all attributes are equal
        """
        if not isinstance(col, BaseColumn):
            raise ValueError('Comparison `col` must be a Column or '
                             'MaskedColumn object')

        attrs = ('name', 'unit', 'dtype', 'format', 'description', 'meta')
        equal = all(getattr(self, x) == getattr(col, x) for x in attrs)

        return equal

    @property
    def _formatter(self):
        return FORMATTER if (self.parent_table is None) else self.parent_table.formatter

    def pformat(self, max_lines=None, show_name=True, show_unit=False, show_dtype=False,
                html=False):
        """Return a list of formatted string representation of column values.

        If no value of ``max_lines`` is supplied then the height of the
        screen terminal is used to set ``max_lines``.  If the terminal
        height cannot be determined then the default will be
        determined using the ``astropy.conf.max_lines`` configuration
        item. If a negative value of ``max_lines`` is supplied then
        there is no line limit applied.

        Parameters
        ----------
        max_lines : int
            Maximum lines of output (header + data rows)

        show_name : bool
            Include column name (default=True)

        show_unit : bool
            Include a header row for unit (default=False)

        show_dtype : bool
            Include column dtype (default=False)

        html : bool
            Format the output as an HTML table (default=False)

        Returns
        -------
        lines : list
            List of lines with header and formatted column values

        """
        _pformat_col = self._formatter._pformat_col
        lines, outs = _pformat_col(self, max_lines, show_name=show_name,
                                   show_unit=show_unit, show_dtype=show_dtype,
                                   html=html)
        return lines

    def pprint(self, max_lines=None, show_name=True, show_unit=False, show_dtype=False):
        """Print a formatted string representation of column values.

        If no value of ``max_lines`` is supplied then the height of the
        screen terminal is used to set ``max_lines``.  If the terminal
        height cannot be determined then the default will be
        determined using the ``astropy.conf.max_lines`` configuration
        item. If a negative value of ``max_lines`` is supplied then
        there is no line limit applied.

        Parameters
        ----------
        max_lines : int
            Maximum number of values in output

        show_name : bool
            Include column name (default=True)

        show_unit : bool
            Include a header row for unit (default=False)

        show_dtype : bool
            Include column dtype (default=True)
        """
        _pformat_col = self._formatter._pformat_col
        lines, outs = _pformat_col(self, max_lines, show_name=show_name, show_unit=show_unit,
                                   show_dtype=show_dtype)

        n_header = outs['n_header']
        for i, line in enumerate(lines):
            if i < n_header:
                color_print(line, 'red')
            else:
                print(line)

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
        _more_tabcol = self._formatter._more_tabcol
        _more_tabcol(self, max_lines=max_lines, show_name=show_name,
                     show_unit=show_unit)

    @property
    def unit(self):
        """
        The unit associated with this column.  May be a string or a
        `astropy.units.UnitBase` instance.

        Setting the ``unit`` property does not change the values of the
        data.  To perform a unit conversion, use ``convert_unit_to``.
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

    def convert_unit_to(self, new_unit, equivalencies=[]):
        """
        Converts the values of the column in-place from the current
        unit to the given unit.

        To change the unit associated with this column without
        actually changing the data values, simply set the ``unit``
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

    @property
    def groups(self):
        if not hasattr(self, '_groups'):
            self._groups = groups.ColumnGroups(self)
        return self._groups

    def group_by(self, keys):
        """
        Group this column by the specified ``keys``

        This effectively splits the column into groups which correspond to
        unique values of the ``keys`` grouping object.  The output is a new
        `Column` or `MaskedColumn` which contains a copy of this column but
        sorted by row according to ``keys``.

        The ``keys`` input to ``group_by`` must be a numpy array with the
        same length as this column.

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

    # Strip off the BaseColumn-ness for repr and str so that
    # MaskedColumn.data __repr__ does not include masked_BaseColumn(data =
    # [1 2], ...).
    def __repr__(self):
        return np.asarray(self).__repr__()

    @property
    def quantity(self):
        """
        A view of this table column as a `~astropy.units.Quantity` object with
        units given by the Column's `unit` parameter.
        """
        # the Quantity initializer is used here because it correctly fails
        # if the column's values are non-numeric (like strings), while .view
        # will happily return a quantity with gibberish for numerical values
        return Quantity(self, copy=False, dtype=self.dtype, order='A')

    def to(self, unit, equivalencies=[], **kwargs):
        """
        Converts this table column to a `~astropy.units.Quantity` object with
        the requested units.

        Parameters
        ----------
        unit : `~astropy.units.Unit` or str
            The unit to convert to (i.e., a valid argument to the
            :meth:`astropy.units.Quantity.to` method).
        equivalencies : list of equivalence pairs, optional
            Equivalencies to use for this conversion.  See
            :meth:`astropy.units.Quantity.to` for more details.

        Returns
        -------
        quantity : `~astropy.units.Quantity`
            A quantity object with the contents of this column in the units
            ``unit``.
        """
        return self.quantity.to(unit, equivalencies)

    def _copy_attrs(self, obj):
        """
        Copy key column attributes from ``obj`` to self
        """
        for attr in ('name', 'unit', 'format', 'description'):
            val = getattr(obj, attr, None)
            setattr(self, attr, val)
        self.meta = deepcopy(getattr(obj, 'meta', {}))


class Column(BaseColumn):
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
    format : str or None or function or callable
        Format string for outputting column values.  This can be an
        "old-style" (``format % value``) or "new-style" (`str.format`)
        format specification string or a function or any callable object that
        accepts a single value and returns a string.
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
    """

    def __new__(cls, data=None, name=None,
                dtype=None, shape=(), length=0,
                description=None, unit=None, format=None, meta=None, copy=False):

        if isinstance(data, MaskedColumn) and np.any(data.mask):
            raise TypeError("Cannot convert a MaskedColumn with masked value to a Column")

        self = super(Column, cls).__new__(cls, data=data, name=name, dtype=dtype,
                                          shape=shape, length=length, description=description,
                                          unit=unit, format=format, meta=meta, copy=copy)
        return self

    def _base_repr_(self, html=False):
        # If scalar then just convert to correct numpy type and use numpy repr
        if self.ndim == 0:
            return repr(self.item())

        descr_vals = [self.__class__.__name__]
        unit = None if self.unit is None else str(self.unit)
        shape = None if self.ndim <= 1 else self.shape[1:]
        for attr, val in (('name', self.name),
                          ('dtype', self.dtype.name),
                          ('shape', shape),
                          ('unit', unit),
                          ('format', self.format),
                          ('description', self.description),
                          ('length', len(self))):

            if val is not None:
                descr_vals.append('{0}={1}'.format(attr, repr(val)))

        descr = '<' + ' '.join(descr_vals) + '>\n'

        if html:
            from ..utils.xml.writer import xml_escape
            descr = xml_escape(descr)

        data_lines, outs = self._formatter._pformat_col(
            self, show_name=False, show_unit=False, show_length=False, html=html)

        out = descr + '\n'.join(data_lines)
        if six.PY2 and isinstance(out, six.text_type):
            out = out.encode('utf-8')

        return out

    def _repr_html_(self):
        return self._base_repr_(html=True)

    def __repr__(self):
        return self._base_repr_(html=False)

    def __unicode__(self):
        # If scalar then just convert to correct numpy type and use numpy repr
        if self.ndim == 0:
            return str(self.item())

        lines, outs = self._formatter._pformat_col(self)
        return '\n'.join(lines)
    if six.PY3:
        __str__ = __unicode__

    def __bytes__(self):
        return six.text_type(self).encode('utf-8')
    if six.PY2:
        __str__ = __bytes__

    # Set items using a view of the underlying data, as it gives an
    # order-of-magnitude speed-up. [#2994]
    def __setitem__(self, index, value):
        self.data[index] = value

    # # Set slices using a view of the underlying data, as it gives an
    # # order-of-magnitude speed-up.  Only gets called in Python 2.  [#3020]
    def __setslice__(self, start, stop, value):
        self.data.__setslice__(start, stop, value)

    def insert(self, obj, values):
        """
        Insert values before the given indices in the column and return
        a new `~astropy.table.Column` object.

        Parameters
        ----------
        obj : int, slice or sequence of ints
            Object that defines the index or indices before which ``values`` is
            inserted.
        values : array_like
            Value(s) to insert.  If the type of ``values`` is different
            from that of quantity, ``values`` is converted to the matching type.
            ``values`` should be shaped so that it can be broadcast appropriately

        Returns
        -------
        out : `~astropy.table.Column`
            A copy of column with ``values`` and ``mask`` inserted.  Note that the
            insertion does not occur in-place: a new column is returned.
        """
        if self.dtype.kind == 'O':
            # Even if values is array-like (e.g. [1,2,3]), insert as a single
            # object.  Numpy.insert instead inserts each element in an array-like
            # input individually.
            data = np.insert(self, obj, None, axis=0)
            data[obj] = values
        else:
            # Explicitly convert to dtype of this column.  Needed because numpy 1.7
            # enforces safe casting by default, so .  This isn't the case for 1.6 or 1.8+.
            values = np.asarray(values, dtype=self.dtype)
            data = np.insert(self, obj, values, axis=0)
        out = data.view(self.__class__)
        out.__array_finalize__(self)
        return out

    # We do this to make the methods show up in the API docs
    name = BaseColumn.name
    unit = BaseColumn.unit
    copy = BaseColumn.copy
    more = BaseColumn.more
    pprint = BaseColumn.pprint
    pformat = BaseColumn.pformat
    convert_unit_to = BaseColumn.convert_unit_to
    quantity = BaseColumn.quantity
    to = BaseColumn.to


class MaskedColumn(Column, ma.MaskedArray):
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
    format : str or None or function or callable
        Format string for outputting column values.  This can be an
        "old-style" (``format % value``) or "new-style" (`str.format`)
        format specification string or a function or any callable object that
        accepts a single value and returns a string.
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
    """

    def __new__(cls, data=None, name=None, mask=None, fill_value=None,
                dtype=None, shape=(), length=0,
                description=None, unit=None, format=None, meta=None, copy=False):

        if mask is None and hasattr(data, 'mask'):
            mask = data.mask
        else:
            mask = deepcopy(mask)

        # Create self using MaskedArray as a wrapper class, following the example of
        # class MSubArray in
        # https://github.com/numpy/numpy/blob/maintenance/1.8.x/numpy/ma/tests/test_subclassing.py
        # This pattern makes it so that __array_finalize__ is called as expected (e.g. #1471 and
        # https://github.com/astropy/astropy/commit/ff6039e8)

        # First just pass through all args and kwargs to BaseColumn, then wrap that object
        # with MaskedArray.
        self_data = BaseColumn(data, dtype=dtype, shape=shape, length=length, name=name,
                               unit=unit, format=format, description=description, meta=meta, copy=copy)
        self = ma.MaskedArray.__new__(cls, data=self_data, mask=mask)

        # Note: do not set fill_value in the MaskedArray constructor because this does not
        # go through the fill_value workarounds (see _fix_fill_value below).
        if fill_value is None and hasattr(data, 'fill_value'):
            fill_value = data.fill_value
        self.fill_value = fill_value

        self.parent_table = None

        return self

    def _fix_fill_value(self, val):
        """Fix a fill value (if needed) to work around a bug with setting the fill
        value of a string array in MaskedArray with Python 3.x.  See
        https://github.com/numpy/numpy/pull/2733.  This mimics the check in
        numpy.ma.core._check_fill_value() (version < 1.8) which incorrectly sets
        fill_value to a default if self.dtype.char is 'U' (which is the case for Python
        3).  Here we change the string to a byte string so that in Python 3 the
        isinstance(val, basestring) part fails.
        """

        if (NUMPY_LT_1_8 and isinstance(val, six.string_types) and
                (self.dtype.char not in 'SV')):
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
            The value to use for invalid entries (`None` by default).  If
            `None`, the ``fill_value`` attribute of the array is used
            instead.

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
        # Use parent table definition of Column if available
        column_cls = self.parent_table.Column if (self.parent_table is not None) else Column
        out = column_cls(name=self.name, data=data, unit=self.unit,
                         format=self.format, description=self.description,
                         meta=deepcopy(self.meta))
        return out

    def insert(self, obj, values, mask=None):
        """
        Insert values along the given axis before the given indices and return
        a new `~astropy.table.MaskedColumn` object.

        Parameters
        ----------
        obj : int, slice or sequence of ints
            Object that defines the index or indices before which ``values`` is
            inserted.
        values : array_like
            Value(s) to insert.  If the type of ``values`` is different
            from that of quantity, ``values`` is converted to the matching type.
            ``values`` should be shaped so that it can be broadcast appropriately
        mask : boolean array_like
            Mask value(s) to insert.  If not supplied then False is used.

        Returns
        -------
        out : `~astropy.table.MaskedColumn`
            A copy of column with ``values`` and ``mask`` inserted.  Note that the
            insertion does not occur in-place: a new masked column is returned.
        """
        self_ma = self.data  # self viewed as MaskedArray

        if self.dtype.kind == 'O':
            # Even if values is array-like (e.g. [1,2,3]), insert as a single
            # object.  Numpy.insert instead inserts each element in an array-like
            # input individually.
            new_data = np.insert(self_ma.data, obj, None, axis=0)
            new_data[obj] = values
        else:
            # Explicitly convert to dtype of this column.  Needed because numpy 1.7
            # enforces safe casting by default, so .  This isn't the case for 1.6 or 1.8+.
            values = np.asarray(values, dtype=self.dtype)
            new_data = np.insert(self_ma.data, obj, values, axis=0)

        if mask is None:
            if self.dtype.kind == 'O':
                mask = False
            else:
                mask = np.zeros(values.shape, dtype=np.bool)
        new_mask = np.insert(self_ma.mask, obj, mask, axis=0)
        new_ma = np.ma.array(new_data, mask=new_mask, copy=False)

        out = new_ma.view(self.__class__)
        out.parent_table = None
        out._copy_attrs(self)

        return out

    def __getitem__(self, item):
        out = super(MaskedColumn, self).__getitem__(item)

        # Fixes issue #3023: when calling getitem with a MaskedArray subclass
        # the original object attributes are not copied.
        if out.__class__ is self.__class__:
            out.parent_table = None
            out._copy_attrs(self)

        return out

    # Set items and slices using MaskedArray method, instead of falling through
    # to the (faster) Column version which uses an ndarray view.  This doesn't
    # copy the mask properly. See test_setting_from_masked_column test.
    def __setitem__(self, index, value):
        ma.MaskedArray.__setitem__(self, index, value)

    def __setslice__(self, start, stop, value):
        ma.MaskedArray.__setslice__(self, start, stop, value)

    # We do this to make the methods show up in the API docs
    name = BaseColumn.name
    copy = BaseColumn.copy
    more = BaseColumn.more
    pprint = BaseColumn.pprint
    pformat = BaseColumn.pformat
    convert_unit_to = BaseColumn.convert_unit_to

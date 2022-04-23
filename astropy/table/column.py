# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools
import warnings
import weakref

from copy import deepcopy

import numpy as np
from numpy import ma

from astropy.units import Unit, Quantity, StructuredUnit
from astropy.utils.console import color_print
from astropy.utils.metadata import MetaData
from astropy.utils.data_info import BaseColumnInfo, dtype_info_name
from astropy.utils.misc import dtype_bytes_or_chars
from . import groups
from . import pprint

# These "shims" provide __getitem__ implementations for Column and MaskedColumn
from ._column_mixins import _ColumnGetitemShim, _MaskedColumnGetitemShim

# Create a generic TableFormatter object for use by bare columns with no
# parent table.
FORMATTER = pprint.TableFormatter()


class StringTruncateWarning(UserWarning):
    """
    Warning class for when a string column is assigned a value
    that gets truncated because the base (numpy) string length
    is too short.

    This does not inherit from AstropyWarning because we want to use
    stacklevel=2 to show the user where the issue occurred in their code.
    """
    pass


# Always emit this warning, not just the first instance
warnings.simplefilter('always', StringTruncateWarning)


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


def col_copy(col, copy_indices=True):
    """
    Mixin-safe version of Column.copy() (with copy_data=True).

    Parameters
    ----------
    col : Column or mixin column
        Input column
    copy_indices : bool
        Copy the column ``indices`` attribute

    Returns
    -------
    col : Copy of input column
    """
    if isinstance(col, BaseColumn):
        return col.copy()

    newcol = col.copy() if hasattr(col, 'copy') else deepcopy(col)
    # If the column has info defined, we copy it and adjust any indices
    # to point to the copied column.  By guarding with the if statement,
    # we avoid side effects (of creating the default info instance).
    if 'info' in col.__dict__:
        newcol.info = col.info
        if copy_indices and col.info.indices:
            newcol.info.indices = deepcopy(col.info.indices)
            for index in newcol.info.indices:
                index.replace_col(col, newcol)

    return newcol


class FalseArray(np.ndarray):
    """
    Boolean mask array that is always False.

    This is used to create a stub ``mask`` property which is a boolean array of
    ``False`` used by default for mixin columns and corresponding to the mixin
    column data shape.  The ``mask`` looks like a normal numpy array but an
    exception will be raised if ``True`` is assigned to any element.  The
    consequences of the limitation are most obvious in the high-level table
    operations.

    Parameters
    ----------
    shape : tuple
        Data shape
    """
    def __new__(cls, shape):
        obj = np.zeros(shape, dtype=bool).view(cls)
        return obj

    def __setitem__(self, item, val):
        val = np.asarray(val)
        if np.any(val):
            raise ValueError('Cannot set any element of {} class to True'
                             .format(self.__class__.__name__))


def _expand_string_array_for_values(arr, values):
    """
    For string-dtype return a version of ``arr`` that is wide enough for ``values``.
    If ``arr`` is not string-dtype or does not need expansion then return ``arr``.

    Parameters
    ----------
    arr : np.ndarray
        Input array
    values : scalar or array-like
        Values for width comparison for string arrays

    Returns
    -------
    arr_expanded : np.ndarray

    """
    if arr.dtype.kind in ('U', 'S') and values is not np.ma.masked:
        # Find the length of the longest string in the new values.
        values_str_len = np.char.str_len(values).max()

        # Determine character repeat count of arr.dtype.  Returns a positive
        # int or None (something like 'U0' is not possible in numpy).  If new values
        # are longer than current then make a new (wider) version of arr.
        arr_str_len = dtype_bytes_or_chars(arr.dtype)
        if arr_str_len and values_str_len > arr_str_len:
            arr_dtype = arr.dtype.byteorder + arr.dtype.kind + str(values_str_len)
            arr = arr.astype(arr_dtype)

    return arr


def _convert_sequence_data_to_array(data, dtype=None):
    """Convert N-d sequence-like data to ndarray or MaskedArray.

    This is the core function for converting Python lists or list of lists to a
    numpy array. This handles embedded np.ma.masked constants in ``data`` along
    with the special case of an homogeneous list of MaskedArray elements.

    Considerations:

    - np.ma.array is about 50 times slower than np.array for list input. This
      function avoids using np.ma.array on list input.
    - np.array emits a UserWarning for embedded np.ma.masked, but only for int
      or float inputs. For those it converts to np.nan and forces float dtype.
      For other types np.array is inconsistent, for instance converting
      np.ma.masked to "0.0" for str types.
    - Searching in pure Python for np.ma.masked in ``data`` is comparable in
      speed to calling ``np.array(data)``.
    - This function may end up making two additional copies of input ``data``.

    Parameters
    ----------
    data : N-d sequence
        Input data, typically list or list of lists
    dtype : None or dtype-like
        Output datatype (None lets np.array choose)

    Returns
    -------
    np_data : np.ndarray or np.ma.MaskedArray

    """
    np_ma_masked = np.ma.masked  # Avoid repeated lookups of this object

    # Special case of an homogeneous list of MaskedArray elements (see #8977).
    # np.ma.masked is an instance of MaskedArray, so exclude those values.
    if (hasattr(data, '__len__')
        and len(data) > 0
        and all(isinstance(val, np.ma.MaskedArray)
                and val is not np_ma_masked for val in data)):
        np_data = np.ma.array(data, dtype=dtype)
        return np_data

    # First convert data to a plain ndarray. If there are instances of np.ma.masked
    # in the data this will issue a warning for int and float.
    with warnings.catch_warnings(record=True) as warns:
        # Ensure this warning from numpy is always enabled and that it is not
        # converted to an error (which can happen during pytest).
        warnings.filterwarnings('always', category=UserWarning,
                                message='.*converting a masked element.*')
        # FutureWarning in numpy 1.21. See https://github.com/astropy/astropy/issues/11291
        # and https://github.com/numpy/numpy/issues/18425.
        warnings.filterwarnings('always', category=FutureWarning,
                                message='.*Promotion of numbers and bools to strings.*')
        try:
            np_data = np.array(data, dtype=dtype)
        except np.ma.MaskError:
            # Catches case of dtype=int with masked values, instead let it
            # convert to float
            np_data = np.array(data)
        except Exception:
            # Conversion failed for some reason, e.g. [2, 1*u.m] gives TypeError in Quantity.
            # First try to interpret the data as Quantity. If that still fails then fall
            # through to object
            try:
                np_data = Quantity(data, dtype)
            except Exception:
                dtype = object
                np_data = np.array(data, dtype=dtype)

    if np_data.ndim == 0 or (np_data.ndim > 0 and len(np_data) == 0):
        # Implies input was a scalar or an empty list (e.g. initializing an
        # empty table with pre-declared names and dtypes but no data).  Here we
        # need to fall through to initializing with the original data=[].
        return data

    # If there were no warnings and the data are int or float, then we are done.
    # Other dtypes like string or complex can have masked values and the
    # np.array() conversion gives the wrong answer (e.g. converting np.ma.masked
    # to the string "0.0").
    if len(warns) == 0 and np_data.dtype.kind in ('i', 'f'):
        return np_data

    # Now we need to determine if there is an np.ma.masked anywhere in input data.

    # Make a statement like below to look for np.ma.masked in a nested sequence.
    # Because np.array(data) succeeded we know that `data` has a regular N-d
    # structure. Find ma_masked:
    #   any(any(any(d2 is ma_masked for d2 in d1) for d1 in d0) for d0 in data)
    # Using this eval avoids creating a copy of `data` in the more-usual case of
    # no masked elements.
    any_statement = 'd0 is ma_masked'
    for ii in reversed(range(np_data.ndim)):
        if ii == 0:
            any_statement = f'any({any_statement} for d0 in data)'
        elif ii == np_data.ndim - 1:
            any_statement = f'any(d{ii} is ma_masked for d{ii} in d{ii-1})'
        else:
            any_statement = f'any({any_statement} for d{ii} in d{ii-1})'
    context = {'ma_masked': np.ma.masked, 'data': data}
    has_masked = eval(any_statement, context)

    # If there are any masks then explicitly change each one to a fill value and
    # set a mask boolean array. If not has_masked then we're done.
    if has_masked:
        mask = np.zeros(np_data.shape, dtype=bool)
        data_filled = np.array(data, dtype=object)

        # Make type-appropriate fill value based on initial conversion.
        if np_data.dtype.kind == 'U':
            fill = ''
        elif np_data.dtype.kind == 'S':
            fill = b''
        else:
            # Zero works for every numeric type.
            fill = 0

        ranges = [range(dim) for dim in np_data.shape]
        for idxs in itertools.product(*ranges):
            val = data_filled[idxs]
            if val is np_ma_masked:
                data_filled[idxs] = fill
                mask[idxs] = True
            elif isinstance(val, bool) and dtype is None:
                # If we see a bool and dtype not specified then assume bool for
                # the entire array. Not perfect but in most practical cases OK.
                # Unfortunately numpy types [False, 0] as int, not bool (and
                # [False, np.ma.masked] => array([0.0, np.nan])).
                dtype = bool

        # If no dtype is provided then need to convert back to list so np.array
        # does type autodetection.
        if dtype is None:
            data_filled = data_filled.tolist()

        # Use np.array first to convert `data` to ndarray (fast) and then make
        # masked array from an ndarray with mask (fast) instead of from `data`.
        np_data = np.ma.array(np.array(data_filled, dtype=dtype), mask=mask)

    return np_data


def _make_compare(oper):
    """
    Make Column comparison methods which encode the ``other`` object to utf-8
    in the case of a bytestring dtype for Py3+.

    Parameters
    ----------
    oper : str
        Operator name
    """
    swapped_oper = {'__eq__': '__eq__',
                    '__ne__': '__ne__',
                    '__gt__': '__lt__',
                    '__lt__': '__gt__',
                    '__ge__': '__le__',
                    '__le__': '__ge__'}[oper]

    def _compare(self, other):
        op = oper  # copy enclosed ref to allow swap below

        # Special case to work around #6838.  Other combinations work OK,
        # see tests.test_column.test_unicode_sandwich_compare().  In this
        # case just swap self and other.
        #
        # This is related to an issue in numpy that was addressed in np 1.13.
        # However that fix does not make this problem go away, but maybe
        # future numpy versions will do so.  NUMPY_LT_1_13 to get the
        # attention of future maintainers to check (by deleting or versioning
        # the if block below).  See #6899 discussion.
        # 2019-06-21: still needed with numpy 1.16.
        if (isinstance(self, MaskedColumn) and self.dtype.kind == 'U'
                and isinstance(other, MaskedColumn) and other.dtype.kind == 'S'):
            self, other = other, self
            op = swapped_oper

        if self.dtype.char == 'S':
            other = self._encode_str(other)

        # Now just let the regular ndarray.__eq__, etc., take over.
        result = getattr(super(Column, self), op)(other)
        # But we should not return Column instances for this case.
        return result.data if isinstance(result, Column) else result

    return _compare


class ColumnInfo(BaseColumnInfo):
    """
    Container for meta information like name, description, format.

    This is required when the object is used as a mixin column within a table,
    but can be used as a general way to store meta information.
    """
    attr_names = BaseColumnInfo.attr_names | {'groups'}
    _attrs_no_copy = BaseColumnInfo._attrs_no_copy | {'groups'}
    attrs_from_parent = attr_names
    _supports_indexing = True
    # For structured columns, data is used to store a dict of columns.
    # Store entries in that dict as name.key instead of name.data.key.
    _represent_as_dict_primary_data = 'data'

    def _represent_as_dict(self):
        result = super()._represent_as_dict()
        names = self._parent.dtype.names
        # For a regular column, we are done, but for a structured
        # column, we use a SerializedColumns to store the pieces.
        if names is None:
            return result

        from .serialize import SerializedColumn

        data = SerializedColumn()
        # If this column has a StructuredUnit, we split it and store
        # it on the corresponding part. Otherwise, we just store it
        # as an attribute below.  All other attributes we remove from
        # the parts, so that we do not store them multiple times.
        # (Note that attributes are not linked to the parent, so it
        # is safe to reset them.)
        # TODO: deal with (some of) this in Column.__getitem__?
        # Alternatively: should we store info on the first part?
        # TODO: special-case format somehow? Can we have good formats
        # for structured columns?
        unit = self.unit
        if isinstance(unit, StructuredUnit) and len(unit) == len(names):
            units = unit.values()
            unit = None  # No need to store as an attribute as well.
        else:
            units = [None] * len(names)
        for name, part_unit in zip(names, units):
            part = self._parent[name]
            part.unit = part_unit
            part.description = None
            part.meta = {}
            part.format = None
            data[name] = part

        # Create the attributes required to reconstruct the column.
        result['data'] = data
        # Store the shape if needed. Just like scalar data, a structured data
        # column (e.g. with dtype `f8,i8`) can be multidimensional within each
        # row and have a shape, and that needs to be distinguished from the
        # case that each entry in the structure has the same shape (e.g.,
        # distinguist a column with dtype='f8,i8' and 2 elements per row from
        # one with dtype '2f8,2i8' and just one element per row).
        if shape := self._parent.shape[1:]:
            result['shape'] = list(shape)
        # Also store the standard info attributes since these are
        # stored on the parent and can thus just be passed on as
        # arguments.  TODO: factor out with essentially the same
        # code in serialize._represent_mixin_as_column.
        if unit is not None and unit != '':
            result['unit'] = unit
        if self.format is not None:
            result['format'] = self.format
        if self.description is not None:
            result['description'] = self.description
        if self.meta:
            result['meta'] = self.meta

        return result

    def _construct_from_dict(self, map):
        if not isinstance(map.get('data'), dict):
            return super()._construct_from_dict(map)

        # Reconstruct a structured Column, by first making an empty column
        # and then filling it with the structured data.
        data = map.pop('data')
        shape = tuple(map.pop('shape', ()))
        # There are three elements in the shape of `part`:
        # (table length, shape of structured column, shape of part like '3f8')
        # The column `shape` only includes the second, so by adding one to its
        # length to include the table length, we pick off a possible last bit.
        dtype = np.dtype([(name, part.dtype, part.shape[len(shape)+1:])
                          for name, part in data.items()])
        units = tuple(col.info.unit for col in data.values())
        if all(unit is not None for unit in units):
            map['unit'] = StructuredUnit(units, dtype)
        map.update(dtype=dtype, shape=shape, length=len(data[dtype.names[0]]))
        # Construct the empty column from `map` (note: 'data' removed above).
        result = super()._construct_from_dict(map)
        # Fill it with the structured data.
        for name in dtype.names:
            result[name] = data[name]
        return result

    def new_like(self, cols, length, metadata_conflicts='warn', name=None):
        """
        Return a new Column instance which is consistent with the
        input ``cols`` and has ``length`` rows.

        This is intended for creating an empty column object whose elements can
        be set in-place for table operations like join or vstack.

        Parameters
        ----------
        cols : list
            List of input columns
        length : int
            Length of the output column object
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output column name

        Returns
        -------
        col : Column (or subclass)
            New instance of this class consistent with ``cols``

        """
        attrs = self.merge_cols_attributes(cols, metadata_conflicts, name,
                                           ('meta', 'unit', 'format', 'description'))

        return self._parent_cls(length=length, **attrs)

    def get_sortable_arrays(self):
        """
        Return a list of arrays which can be lexically sorted to represent
        the order of the parent column.

        For Column this is just the column itself.

        Returns
        -------
        arrays : list of ndarray
        """
        return [self._parent]


class BaseColumn(_ColumnGetitemShim, np.ndarray):

    meta = MetaData()

    def __new__(cls, data=None, name=None,
                dtype=None, shape=(), length=0,
                description=None, unit=None, format=None, meta=None,
                copy=False, copy_indices=True):
        if data is None:
            self_data = np.zeros((length,)+shape, dtype=dtype)
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
                meta = data.meta
            if name is None:
                name = data.name
        elif isinstance(data, Quantity):
            if unit is None:
                self_data = np.array(data, dtype=dtype, copy=copy)
                unit = data.unit
            else:
                self_data = Quantity(data, unit, dtype=dtype, copy=copy).value
            # If 'info' has been defined, copy basic properties (if needed).
            if 'info' in data.__dict__:
                if description is None:
                    description = data.info.description
                if format is None:
                    format = data.info.format
                if meta is None:
                    meta = data.info.meta

        else:
            if np.dtype(dtype).char == 'S':
                data = cls._encode_str(data)
            self_data = np.array(data, dtype=dtype, copy=copy)

        self = self_data.view(cls)
        self._name = None if name is None else str(name)
        self._parent_table = None
        self.unit = unit
        self._format = format
        self.description = description
        self.meta = meta
        self.indices = deepcopy(getattr(data, 'indices', [])) if copy_indices else []
        for index in self.indices:
            index.replace_col(data, self)

        return self

    @property
    def data(self):
        return self.view(np.ndarray)

    @property
    def value(self):
        """
        An alias for the existing ``data`` attribute.
        """
        return self.data

    @property
    def parent_table(self):
        # Note: It seems there are some cases where _parent_table is not set,
        # such after restoring from a pickled Column.  Perhaps that should be
        # fixed, but this is also okay for now.
        if getattr(self, '_parent_table', None) is None:
            return None
        else:
            return self._parent_table()

    @parent_table.setter
    def parent_table(self, table):
        if table is None:
            self._parent_table = None
        else:
            self._parent_table = weakref.ref(table)

    info = ColumnInfo()

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

        # If there is meta on the original column then deepcopy (since "copy" of column
        # implies complete independence from original).  __array_finalize__ will have already
        # made a light copy.  I'm not sure how to avoid that initial light copy.
        if self.meta is not None:
            out.meta = self.meta  # MetaData descriptor does a deepcopy here

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
        # Get the Column attributes
        names = ('_name', '_unit', '_format', 'description', 'meta', 'indices')
        attrs = {name: val for name, val in zip(names, state[-1])}

        state = state[:-1]

        # Using super().__setstate__(state) gives
        # "TypeError 'int' object is not iterable", raised in
        # astropy.table._column_mixins._ColumnGetitemShim.__setstate_cython__()
        # Previously, it seems to have given an infinite recursion.
        # Hence, manually call the right super class to actually set up
        # the array object.
        super_class = ma.MaskedArray if isinstance(self, ma.MaskedArray) else np.ndarray
        super_class.__setstate__(self, state)

        # Set the Column attributes
        for name, val in attrs.items():
            setattr(self, name, val)
        self._parent_table = None

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
                        self.meta, self.indices)
        state = state + (column_state,)

        return reconstruct_func, reconstruct_func_args, state

    def __array_finalize__(self, obj):
        # Obj will be none for direct call to Column() creator
        if obj is None:
            return

        if callable(super().__array_finalize__):
            super().__array_finalize__(obj)

        # Self was created from template (e.g. obj[slice] or (obj * 2))
        # or viewcast e.g. obj.view(Column).  In either case we want to
        # init Column attributes for self from obj if possible.
        self.parent_table = None
        if not hasattr(self, 'indices'):  # may have been copied in __new__
            self.indices = []
        self._copy_attrs(obj)
        if 'info' in getattr(obj, '__dict__', {}):
            self.info = obj.info

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
        out_arr = super().__array_wrap__(out_arr, context)
        if (self.shape != out_arr.shape
            or (isinstance(out_arr, BaseColumn)
                and (context is not None
                     and context[0] in _comparison_functions))):
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
        if val is not None:
            val = str(val)

        if self.parent_table is not None:
            table = self.parent_table
            table.columns._rename_column(self.name, val)

        self._name = val

    @property
    def format(self):
        """
        Format string for displaying values in this column.
        """

        return self._format

    @format.setter
    def format(self, format_string):

        prev_format = getattr(self, '_format', None)

        self._format = format_string  # set new format string

        try:
            # test whether it formats without error exemplarily
            self.pformat(max_lines=1)
        except Exception as err:
            # revert to restore previous format if there was one
            self._format = prev_format
            raise ValueError(
                "Invalid format for column '{}': could not display "
                "values in this column using this format".format(
                    self.name)) from err

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
        equal : bool
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
            Include column name. Default is True.

        show_unit : bool
            Include a header row for unit. Default is False.

        show_dtype : bool
            Include column dtype. Default is False.

        html : bool
            Format the output as an HTML table. Default is False.

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
            Include column name. Default is True.

        show_unit : bool
            Include a header row for unit. Default is False.

        show_dtype : bool
            Include column dtype. Default is True.
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
            Maximum number of lines in table output.

        show_name : bool
            Include a header row for column names. Default is True.

        show_unit : bool
            Include a header row for unit. Default is False.

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

    def searchsorted(self, v, side='left', sorter=None):
        # For bytes type data, encode the `v` value as UTF-8 (if necessary) before
        # calling searchsorted. This prevents a factor of 1000 slowdown in
        # searchsorted in this case.
        a = self.data
        if a.dtype.kind == 'S' and not isinstance(v, bytes):
            v = np.asarray(v)
            if v.dtype.kind == 'U':
                v = np.char.encode(v, 'utf-8')
        return np.searchsorted(a, v, side=side, sorter=sorter)
    searchsorted.__doc__ = np.ndarray.searchsorted.__doc__

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

        equivalencies : list of tuple
           A list of equivalence pairs to try if the unit are not
           directly convertible.  See :ref:`astropy:unit_equivalencies`.

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
        return Quantity(self, self.unit, copy=False, dtype=self.dtype, order='A', subok=True)

    def to(self, unit, equivalencies=[], **kwargs):
        """
        Converts this table column to a `~astropy.units.Quantity` object with
        the requested units.

        Parameters
        ----------
        unit : unit-like
            The unit to convert to (i.e., a valid argument to the
            :meth:`astropy.units.Quantity.to` method).
        equivalencies : list of tuple
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
        for attr in ('name', 'unit', '_format', 'description'):
            val = getattr(obj, attr, None)
            setattr(self, attr, val)

        # Light copy of meta if it is not empty
        obj_meta = getattr(obj, 'meta', None)
        if obj_meta:
            self.meta = obj_meta.copy()

    @staticmethod
    def _encode_str(value):
        """
        Encode anything that is unicode-ish as utf-8.  This method is only
        called for Py3+.
        """
        if isinstance(value, str):
            value = value.encode('utf-8')
        elif isinstance(value, bytes) or value is np.ma.masked:
            pass
        else:
            arr = np.asarray(value)
            if arr.dtype.char == 'U':
                arr = np.char.encode(arr, encoding='utf-8')
                if isinstance(value, np.ma.MaskedArray):
                    arr = np.ma.array(arr, mask=value.mask, copy=False)
            value = arr

        return value

    def tolist(self):
        if self.dtype.kind == 'S':
            return np.chararray.decode(self, encoding='utf-8').tolist()
        else:
            return super().tolist()


class Column(BaseColumn):
    """Define a data column for use in a Table object.

    Parameters
    ----------
    data : list, ndarray, or None
        Column data values
    name : str
        Column name and key for reference within Table
    dtype : `~numpy.dtype`-like
        Data type for column
    shape : tuple or ()
        Dimensions of a single row element in the column data
    length : int or 0
        Number of row elements in column data
    description : str or None
        Full description of column
    unit : str or None
        Physical unit
    format : str, None, or callable
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
      `<https://numpy.org/doc/stable/reference/arrays.dtypes.html>`_.
      Examples include:

      - Python non-string type (float, int, bool)
      - Numpy non-string type (e.g. np.float32, np.int64, np.bool\\_)
      - Numpy.dtype array-protocol type strings (e.g. 'i4', 'f8', 'S15')

      If no ``dtype`` value is provide then the type is inferred using
      ``np.array(data)``.

    - Provide ``length`` and optionally ``shape``, but not ``data``

      Examples::

        col = Column(name='name', length=5)
        col = Column(name='name', dtype=int, length=10, shape=(3,4))

      The default ``dtype`` is ``np.float64``.  The ``shape`` argument is the
      array shape of a single cell in the column.

    To access the ``Column`` data as a raw `numpy.ndarray` object, you can use
    one of the ``data`` or ``value`` attributes (which are equivalent)::

        col.data
        col.value
    """

    def __new__(cls, data=None, name=None,
                dtype=None, shape=(), length=0,
                description=None, unit=None, format=None, meta=None,
                copy=False, copy_indices=True):

        if isinstance(data, MaskedColumn) and np.any(data.mask):
            raise TypeError("Cannot convert a MaskedColumn with masked value to a Column")

        self = super().__new__(
            cls, data=data, name=name, dtype=dtype, shape=shape, length=length,
            description=description, unit=unit, format=format, meta=meta,
            copy=copy, copy_indices=copy_indices)
        return self

    def __setattr__(self, item, value):
        if not isinstance(self, MaskedColumn) and item == "mask":
            raise AttributeError("cannot set mask value to a column in non-masked Table")
        super().__setattr__(item, value)

        if item == 'unit' and issubclass(self.dtype.type, np.number):
            try:
                converted = self.parent_table._convert_col_for_table(self)
            except AttributeError:  # Either no parent table or parent table is None
                pass
            else:
                if converted is not self:
                    self.parent_table.replace_column(self.name, converted)

    def _base_repr_(self, html=False):
        # If scalar then just convert to correct numpy type and use numpy repr
        if self.ndim == 0:
            return repr(self.item())

        descr_vals = [self.__class__.__name__]
        unit = None if self.unit is None else str(self.unit)
        shape = None if self.ndim <= 1 else self.shape[1:]
        for attr, val in (('name', self.name),
                          ('dtype', dtype_info_name(self.dtype)),
                          ('shape', shape),
                          ('unit', unit),
                          ('format', self.format),
                          ('description', self.description),
                          ('length', len(self))):

            if val is not None:
                descr_vals.append(f'{attr}={val!r}')

        descr = '<' + ' '.join(descr_vals) + '>\n'

        if html:
            from astropy.utils.xml.writer import xml_escape
            descr = xml_escape(descr)

        data_lines, outs = self._formatter._pformat_col(
            self, show_name=False, show_unit=False, show_length=False, html=html)

        out = descr + '\n'.join(data_lines)

        return out

    def _repr_html_(self):
        return self._base_repr_(html=True)

    def __repr__(self):
        return self._base_repr_(html=False)

    def __str__(self):
        # If scalar then just convert to correct numpy type and use numpy repr
        if self.ndim == 0:
            return str(self.item())

        lines, outs = self._formatter._pformat_col(self)
        return '\n'.join(lines)

    def __bytes__(self):
        return str(self).encode('utf-8')

    def _check_string_truncate(self, value):
        """
        Emit a warning if any elements of ``value`` will be truncated when
        ``value`` is assigned to self.
        """
        # Convert input ``value`` to the string dtype of this column and
        # find the length of the longest string in the array.
        value = np.asanyarray(value, dtype=self.dtype.type)
        if value.size == 0:
            return
        value_str_len = np.char.str_len(value).max()

        # Parse the array-protocol typestring (e.g. '|U15') of self.dtype which
        # has the character repeat count on the right side.
        self_str_len = dtype_bytes_or_chars(self.dtype)

        if value_str_len > self_str_len:
            warnings.warn('truncated right side string(s) longer than {} '
                          'character(s) during assignment'
                          .format(self_str_len),
                          StringTruncateWarning,
                          stacklevel=3)

    def __setitem__(self, index, value):
        if self.dtype.char == 'S':
            value = self._encode_str(value)

        # Issue warning for string assignment that truncates ``value``
        if issubclass(self.dtype.type, np.character):
            self._check_string_truncate(value)

        # update indices
        self.info.adjust_indices(index, value, len(self))

        # Set items using a view of the underlying data, as it gives an
        # order-of-magnitude speed-up. [#2994]
        self.data[index] = value

    __eq__ = _make_compare('__eq__')
    __ne__ = _make_compare('__ne__')
    __gt__ = _make_compare('__gt__')
    __lt__ = _make_compare('__lt__')
    __ge__ = _make_compare('__ge__')
    __le__ = _make_compare('__le__')

    def insert(self, obj, values, axis=0):
        """
        Insert values before the given indices in the column and return
        a new `~astropy.table.Column` object.

        Parameters
        ----------
        obj : int, slice or sequence of int
            Object that defines the index or indices before which ``values`` is
            inserted.
        values : array-like
            Value(s) to insert.  If the type of ``values`` is different from
            that of the column, ``values`` is converted to the matching type.
            ``values`` should be shaped so that it can be broadcast appropriately.
        axis : int, optional
            Axis along which to insert ``values``.  If ``axis`` is None then
            the column array is flattened before insertion.  Default is 0,
            which will insert a row.

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
            data = np.insert(self, obj, None, axis=axis)
            data[obj] = values
        else:
            self_for_insert = _expand_string_array_for_values(self, values)
            data = np.insert(self_for_insert, obj, values, axis=axis)

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


class MaskedColumnInfo(ColumnInfo):
    """
    Container for meta information like name, description, format.

    This is required when the object is used as a mixin column within a table,
    but can be used as a general way to store meta information.  In this case
    it just adds the ``mask_val`` attribute.
    """
    # Add `serialize_method` attribute to the attrs that MaskedColumnInfo knows
    # about.  This allows customization of the way that MaskedColumn objects
    # get written to file depending on format.  The default is to use whatever
    # the writer would normally do, which in the case of FITS or ECSV is to use
    # a NULL value within the data itself.  If serialize_method is 'data_mask'
    # then the mask is explicitly written out as a separate column if there
    # are any masked values.  See also code below.
    attr_names = ColumnInfo.attr_names | {'serialize_method'}

    # When `serialize_method` is 'data_mask', and data and mask are being written
    # as separate columns, use column names <name> and <name>.mask (instead
    # of default encoding as <name>.data and <name>.mask).
    _represent_as_dict_primary_data = 'data'

    mask_val = np.ma.masked

    def __init__(self, bound=False):
        super().__init__(bound)

        # If bound to a data object instance then create the dict of attributes
        # which stores the info attribute values.
        if bound:
            # Specify how to serialize this object depending on context.
            self.serialize_method = {'fits': 'null_value',
                                     'ecsv': 'null_value',
                                     'hdf5': 'data_mask',
                                     'parquet': 'data_mask',
                                     None: 'null_value'}

    def _represent_as_dict(self):
        out = super()._represent_as_dict()
        # If we are a structured masked column, then our parent class,
        # ColumnInfo, will already have set up a dict with masked parts,
        # which will be serialized later, so no further work needed here.
        if self._parent.dtype.names is not None:
            return out

        col = self._parent

        # If the serialize method for this context (e.g. 'fits' or 'ecsv') is
        # 'data_mask', that means to serialize using an explicit mask column.
        method = self.serialize_method[self._serialize_context]

        if method == 'data_mask':
            # Note: a driver here is a performance issue in #8443 where repr() of a
            # np.ma.MaskedArray value is up to 10 times slower than repr of a normal array
            # value.  So regardless of whether there are masked elements it is useful to
            # explicitly define this as a serialized column and use col.data.data (ndarray)
            # instead of letting it fall through to the "standard" serialization machinery.
            out['data'] = col.data.data

            if np.any(col.mask):
                # Only if there are actually masked elements do we add the ``mask`` column
                out['mask'] = col.mask

        elif method == 'null_value':
            pass

        else:
            raise ValueError('serialize method must be either "data_mask" or "null_value"')

        return out


class MaskedColumn(Column, _MaskedColumnGetitemShim, ma.MaskedArray):
    """Define a masked data column for use in a Table object.

    Parameters
    ----------
    data : list, ndarray, or None
        Column data values
    name : str
        Column name and key for reference within Table
    mask : list, ndarray or None
        Boolean mask for which True indicates missing or invalid data
    fill_value : float, int, str, or None
        Value used when filling masked column elements
    dtype : `~numpy.dtype`-like
        Data type for column
    shape : tuple or ()
        Dimensions of a single row element in the column data
    length : int or 0
        Number of row elements in column data
    description : str or None
        Full description of column
    unit : str or None
        Physical unit
    format : str, None, or callable
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
      `<https://numpy.org/doc/stable/reference/arrays.dtypes.html>`_.
      Examples include:

      - Python non-string type (float, int, bool)
      - Numpy non-string type (e.g. np.float32, np.int64, np.bool\\_)
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

    To access the ``Column`` data as a raw `numpy.ma.MaskedArray` object, you can
    use one of the ``data`` or ``value`` attributes (which are equivalent)::

        col.data
        col.value
    """
    info = MaskedColumnInfo()

    def __new__(cls, data=None, name=None, mask=None, fill_value=None,
                dtype=None, shape=(), length=0,
                description=None, unit=None, format=None, meta=None,
                copy=False, copy_indices=True):

        if mask is None:
            # If mask is None then we need to determine the mask (if any) from the data.
            # The naive method is looking for a mask attribute on data, but this can fail,
            # see #8816.  Instead use ``MaskedArray`` to do the work.
            mask = ma.MaskedArray(data).mask
            if mask is np.ma.nomask:
                # Handle odd-ball issue with np.ma.nomask (numpy #13758), and see below.
                mask = False
            elif copy:
                mask = mask.copy()

        elif mask is np.ma.nomask:
            # Force the creation of a full mask array as nomask is tricky to
            # use and will fail in an unexpected manner when setting a value
            # to the mask.
            mask = False
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
                               unit=unit, format=format, description=description,
                               meta=meta, copy=copy, copy_indices=copy_indices)
        self = ma.MaskedArray.__new__(cls, data=self_data, mask=mask)
        # The above process preserves info relevant for Column, but this does
        # not include serialize_method (and possibly other future attributes)
        # relevant for MaskedColumn, so we set info explicitly.
        if 'info' in getattr(data, '__dict__', {}):
            self.info = data.info

        # Note: do not set fill_value in the MaskedArray constructor because this does not
        # go through the fill_value workarounds.
        if fill_value is None and getattr(data, 'fill_value', None) is not None:
            # Coerce the fill_value to the correct type since `data` may be a
            # different dtype than self.
            fill_value = np.array(data.fill_value, self.dtype)[()]
        self.fill_value = fill_value

        self.parent_table = None

        # needs to be done here since self doesn't come from BaseColumn.__new__
        for index in self.indices:
            index.replace_col(self_data, self)

        return self

    @property
    def fill_value(self):
        return self.get_fill_value()  # defer to native ma.MaskedArray method

    @fill_value.setter
    def fill_value(self, val):
        """Set fill value both in the masked column view and in the parent table
        if it exists.  Setting one or the other alone doesn't work."""

        # another ma bug workaround: If the value of fill_value for a string array is
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
        """The plain MaskedArray data held by this column."""
        out = self.view(np.ma.MaskedArray)
        # By default, a MaskedArray view will set the _baseclass to be the
        # same as that of our own class, i.e., BaseColumn.  Since we want
        # to return a plain MaskedArray, we reset the baseclass accordingly.
        out._baseclass = np.ndarray
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

        data = super().filled(fill_value)
        # Use parent table definition of Column if available
        column_cls = self.parent_table.Column if (self.parent_table is not None) else Column

        out = column_cls(name=self.name, data=data, unit=self.unit,
                         format=self.format, description=self.description,
                         meta=deepcopy(self.meta))
        return out

    def insert(self, obj, values, mask=None, axis=0):
        """
        Insert values along the given axis before the given indices and return
        a new `~astropy.table.MaskedColumn` object.

        Parameters
        ----------
        obj : int, slice or sequence of int
            Object that defines the index or indices before which ``values`` is
            inserted.
        values : array-like
            Value(s) to insert.  If the type of ``values`` is different from
            that of the column, ``values`` is converted to the matching type.
            ``values`` should be shaped so that it can be broadcast appropriately.
        mask : bool or array-like
            Mask value(s) to insert.  If not supplied, and values does not have
            a mask either, then False is used.
        axis : int, optional
            Axis along which to insert ``values``.  If ``axis`` is None then
            the column array is flattened before insertion.  Default is 0,
            which will insert a row.

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
            new_data = np.insert(self_ma.data, obj, None, axis=axis)
            new_data[obj] = values
        else:
            self_ma = _expand_string_array_for_values(self_ma, values)
            new_data = np.insert(self_ma.data, obj, values, axis=axis)

        if mask is None:
            mask = getattr(values, 'mask', np.ma.nomask)
            if mask is np.ma.nomask:
                if self.dtype.kind == 'O':
                    mask = False
                else:
                    mask = np.zeros(np.shape(values), dtype=bool)

        new_mask = np.insert(self_ma.mask, obj, mask, axis=axis)
        new_ma = np.ma.array(new_data, mask=new_mask, copy=False)

        out = new_ma.view(self.__class__)
        out.parent_table = None
        out.indices = []
        out._copy_attrs(self)
        out.fill_value = self.fill_value

        return out

    def _copy_attrs_slice(self, out):
        # Fixes issue #3023: when calling getitem with a MaskedArray subclass
        # the original object attributes are not copied.
        if out.__class__ is self.__class__:
            # TODO: this part is essentially the same as what is done in
            # __array_finalize__ and could probably be called directly in our
            # override of __getitem__ in _columns_mixins.pyx). Refactor?
            if 'info' in self.__dict__:
                out.info = self.info
            out.parent_table = None
            # we need this because __getitem__ does a shallow copy of indices
            if out.indices is self.indices:
                out.indices = []
            out._copy_attrs(self)
        return out

    def __setitem__(self, index, value):
        # Issue warning for string assignment that truncates ``value``
        if self.dtype.char == 'S':
            value = self._encode_str(value)

        if issubclass(self.dtype.type, np.character):
            # Account for a bug in np.ma.MaskedArray setitem.
            # https://github.com/numpy/numpy/issues/8624
            value = np.ma.asanyarray(value, dtype=self.dtype.type)

            # Check for string truncation after filling masked items with
            # empty (zero-length) string.  Note that filled() does not make
            # a copy if there are no masked items.
            self._check_string_truncate(value.filled(''))

        # update indices
        self.info.adjust_indices(index, value, len(self))

        ma.MaskedArray.__setitem__(self, index, value)

    # We do this to make the methods show up in the API docs
    name = BaseColumn.name
    copy = BaseColumn.copy
    more = BaseColumn.more
    pprint = BaseColumn.pprint
    pformat = BaseColumn.pformat
    convert_unit_to = BaseColumn.convert_unit_to

# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""This module contains functions and methods that relate to the DataInfo class
which provides a container for informational attributes as well as summary info
methods.

A DataInfo object is attached to the Quantity, SkyCoord, and Time classes in
astropy.  Here it allows those classes to be used in Tables and uniformly carry
table column attributes such as name, format, dtype, meta, and description.
"""

# Note: these functions and classes are tested extensively in astropy table
# tests via their use in providing mixin column info, and in
# astropy/tests/test_info for providing table and column info summary data.


import os
import re
import sys
import weakref
import warnings
from io import StringIO
from copy import deepcopy
from functools import partial
from collections import OrderedDict
from contextlib import contextmanager

import numpy as np

from . import metadata


__all__ = ['data_info_factory', 'dtype_info_name', 'BaseColumnInfo',
           'DataInfo', 'MixinInfo', 'ParentDtypeInfo']

# Tuple of filterwarnings kwargs to ignore when calling info
IGNORE_WARNINGS = (dict(category=RuntimeWarning, message='All-NaN|'
                        'Mean of empty slice|Degrees of freedom <= 0|'
                        'invalid value encountered in sqrt'),)


@contextmanager
def serialize_context_as(context):
    """Set context for serialization.

    This will allow downstream code to understand the context in which a column
    is being serialized.  Objects like Time or SkyCoord will have different
    default serialization representations depending on context.

    Parameters
    ----------
    context : str
        Context name, e.g. 'fits', 'hdf5', 'parquet', 'ecsv', 'yaml'
    """
    old_context = BaseColumnInfo._serialize_context
    BaseColumnInfo._serialize_context = context
    try:
        yield
    finally:
        BaseColumnInfo._serialize_context = old_context


def dtype_info_name(dtype):
    """Return a human-oriented string name of the ``dtype`` arg.
    This can be use by astropy methods that present type information about
    a data object.

    The output is mostly equivalent to ``dtype.name`` which takes the form
    <type_name>[B] where <type_name> is like ``int`` or ``bool`` and [B] is an
    optional number of bits which gets included only for numeric types.

    The output is shown below for ``bytes`` and ``str`` types, with <N> being
    the number of characters. This representation corresponds to the Python
    type that matches the dtype::

      Numpy          S<N>      U<N>
      Python      bytes<N>   str<N>

    Parameters
    ----------
    dtype : str, `~numpy.dtype`, type
        Input as an object that can be converted via :class:`numpy.dtype`.

    Returns
    -------
    dtype_info_name : str
        String name of ``dtype``
    """
    dtype = np.dtype(dtype)
    if dtype.names is not None:
        return '({})'.format(', '.join(dtype_info_name(dt[0])
                                       for dt in dtype.fields.values()))
    if dtype.subdtype is not None:
        dtype, shape = dtype.subdtype
    else:
        shape = ()

    if dtype.kind in ('S', 'U'):
        type_name = 'bytes' if dtype.kind == 'S' else 'str'
        length = re.search(r'(\d+)', dtype.str).group(1)
        out = type_name + length
    else:
        out = dtype.name

    if shape:
        out += f"[{','.join(str(n) for n in shape)}]"

    return out


def data_info_factory(names, funcs):
    """
    Factory to create a function that can be used as an ``option``
    for outputting data object summary information.

    Examples
    --------
    >>> from astropy.utils.data_info import data_info_factory
    >>> from astropy.table import Column
    >>> c = Column([4., 3., 2., 1.])
    >>> mystats = data_info_factory(names=['min', 'median', 'max'],
    ...                             funcs=[np.min, np.median, np.max])
    >>> c.info(option=mystats)
    min = 1
    median = 2.5
    max = 4
    n_bad = 0
    length = 4

    Parameters
    ----------
    names : list
        List of information attribute names
    funcs : list
        List of functions that compute the corresponding information attribute

    Returns
    -------
    func : function
        Function that can be used as a data info option
    """
    def func(dat):
        outs = []
        for name, func in zip(names, funcs):
            try:
                if isinstance(func, str):
                    out = getattr(dat, func)()
                else:
                    out = func(dat)
            except Exception:
                outs.append('--')
            else:
                try:
                    outs.append(f'{out:g}')
                except (TypeError, ValueError):
                    outs.append(str(out))

        return OrderedDict(zip(names, outs))
    return func


def _get_obj_attrs_map(obj, attrs):
    """
    Get the values for object ``attrs`` and return as a dict.  This
    ignores any attributes that are None.  In the context of serializing
    the supported core astropy classes this conversion will succeed and
    results in more succinct and less python-specific YAML.
    """
    out = {}
    for attr in attrs:
        val = getattr(obj, attr, None)

        if val is not None:
            out[attr] = val
    return out


def _get_data_attribute(dat, attr=None):
    """
    Get a data object attribute for the ``attributes`` info summary method
    """
    if attr == 'class':
        val = type(dat).__name__
    elif attr == 'dtype':
        val = dtype_info_name(dat.info.dtype)
    elif attr == 'shape':
        datshape = dat.shape[1:]
        val = datshape if datshape else ''
    else:
        val = getattr(dat.info, attr)
    if val is None:
        val = ''
    return str(val)


class InfoAttribute:
    def __init__(self, attr, default=None):
        self.attr = attr
        self.default = default

    def __get__(self, instance, owner_cls):
        if instance is None:
            return self

        return instance._attrs.get(self.attr, self.default)

    def __set__(self, instance, value):
        if instance is None:
            # This is an unbound descriptor on the class
            raise ValueError('cannot set unbound descriptor')

        instance._attrs[self.attr] = value


class ParentAttribute:
    def __init__(self, attr):
        self.attr = attr

    def __get__(self, instance, owner_cls):
        if instance is None:
            return self

        return getattr(instance._parent, self.attr)

    def __set__(self, instance, value):
        if instance is None:
            # This is an unbound descriptor on the class
            raise ValueError('cannot set unbound descriptor')

        setattr(instance._parent, self.attr, value)


class DataInfoMeta(type):
    def __new__(mcls, name, bases, dct):
        # Ensure that we do not gain a __dict__, which would mean
        # arbitrary attributes could be set.
        dct.setdefault('__slots__', [])
        return super().__new__(mcls, name, bases, dct)

    def __init__(cls, name, bases, dct):
        super().__init__(name, bases, dct)

        # Define default getters/setters for attributes, if needed.
        for attr in cls.attr_names:
            if attr not in dct:
                # If not defined explicitly for this class, did any of
                # its superclasses define it, and, if so, was this an
                # automatically defined look-up-on-parent attribute?
                cls_attr = getattr(cls, attr, None)
                if attr in cls.attrs_from_parent:
                    # If the attribute is supposed to be stored on the parent,
                    # and that is stated by this class yet it was not the case
                    # on the superclass, override it.
                    if 'attrs_from_parent' in dct and not isinstance(cls_attr, ParentAttribute):
                        setattr(cls, attr, ParentAttribute(attr))
                elif not cls_attr or isinstance(cls_attr, ParentAttribute):
                    # If the attribute is not meant to be stored on the parent,
                    # and if it was not defined already or was previously defined
                    # as an attribute on the parent, define a regular
                    # look-up-on-info attribute
                    setattr(cls, attr,
                            InfoAttribute(attr, cls._attr_defaults.get(attr)))


class DataInfo(metaclass=DataInfoMeta):
    """
    Descriptor that data classes use to add an ``info`` attribute for storing
    data attributes in a uniform and portable way.  Note that it *must* be
    called ``info`` so that the DataInfo() object can be stored in the
    ``instance`` using the ``info`` key.  Because owner_cls.x is a descriptor,
    Python doesn't use __dict__['x'] normally, and the descriptor can safely
    store stuff there.  Thanks to
    https://nbviewer.jupyter.org/urls/gist.github.com/ChrisBeaumont/5758381/raw/descriptor_writeup.ipynb
    for this trick that works for non-hashable classes.

    Parameters
    ----------
    bound : bool
        If True this is a descriptor attribute in a class definition, else it
        is a DataInfo() object that is bound to a data object instance. Default is False.
    """
    _stats = ['mean', 'std', 'min', 'max']
    attrs_from_parent = set()
    attr_names = set(['name', 'unit', 'dtype', 'format', 'description', 'meta'])
    _attr_defaults = {'dtype': np.dtype('O')}
    _attrs_no_copy = set()
    _info_summary_attrs = ('dtype', 'shape', 'unit', 'format', 'description', 'class')
    __slots__ = ['_parent_cls', '_parent_ref', '_attrs']
    # This specifies the list of object attributes which must be stored in
    # order to re-create the object after serialization.  This is independent
    # of normal `info` attributes like name or description.  Subclasses will
    # generally either define this statically (QuantityInfo) or dynamically
    # (SkyCoordInfo).  These attributes may be scalars or arrays.  If arrays
    # that match the object length they will be serialized as an independent
    # column.
    _represent_as_dict_attrs = ()

    # This specifies attributes which are to be provided to the class
    # initializer as ordered args instead of keyword args.  This is needed
    # for Quantity subclasses where the keyword for data varies (e.g.
    # between Quantity and Angle).
    _construct_from_dict_args = ()

    # This specifies the name of an attribute which is the "primary" data.
    # Then when representing as columns
    # (table.serialize._represent_mixin_as_column) the output for this
    # attribute will be written with the just name of the mixin instead of the
    # usual "<name>.<attr>".
    _represent_as_dict_primary_data = None

    def __init__(self, bound=False):
        # If bound to a data object instance then create the dict of attributes
        # which stores the info attribute values. Default of None for "unset"
        # except for dtype where the default is object.
        if bound:
            self._attrs = {}

    @property
    def _parent(self):
        try:
            parent = self._parent_ref()
        except AttributeError:
            return None

        if parent is None:
            raise AttributeError("""\
failed to access "info" attribute on a temporary object.

It looks like you have done something like ``col[3:5].info`` or
``col.quantity.info``, i.e.  you accessed ``info`` from a temporary slice
object that only exists momentarily.  This has failed because the reference to
that temporary object is now lost.  Instead force a permanent reference (e.g.
``c = col[3:5]`` followed by ``c.info``).""")

        return parent

    def __get__(self, instance, owner_cls):
        if instance is None:
            # This is an unbound descriptor on the class
            self._parent_cls = owner_cls
            return self

        info = instance.__dict__.get('info')
        if info is None:
            info = instance.__dict__['info'] = self.__class__(bound=True)
        # We set _parent_ref on every call, since if one makes copies of
        # instances, 'info' will be copied as well, which will lose the
        # reference.
        info._parent_ref = weakref.ref(instance)
        return info

    def __set__(self, instance, value):
        if instance is None:
            # This is an unbound descriptor on the class
            raise ValueError('cannot set unbound descriptor')

        if isinstance(value, DataInfo):
            info = instance.__dict__['info'] = self.__class__(bound=True)
            attr_names = info.attr_names
            if value.__class__ is self.__class__:
                # For same class, attributes are guaranteed to be stored in
                # _attrs, so speed matters up by not accessing defaults.
                # Doing this before difference in for loop helps speed.
                attr_names = attr_names & set(value._attrs)  # NOT in-place!
            else:
                # For different classes, copy over the attributes in common.
                attr_names = attr_names & (value.attr_names - value._attrs_no_copy)

            for attr in attr_names - info.attrs_from_parent - info._attrs_no_copy:
                info._attrs[attr] = deepcopy(getattr(value, attr))

        else:
            raise TypeError('info must be set with a DataInfo instance')

    def __getstate__(self):
        return self._attrs

    def __setstate__(self, state):
        self._attrs = state

    def _represent_as_dict(self, attrs=None):
        """Get the values for the parent ``attrs`` and return as a dict.

        By default, uses '_represent_as_dict_attrs'.
        """
        if attrs is None:
            attrs = self._represent_as_dict_attrs
        return _get_obj_attrs_map(self._parent, attrs)

    def _construct_from_dict(self, map):
        args = [map.pop(attr) for attr in self._construct_from_dict_args]
        return self._parent_cls(*args, **map)

    info_summary_attributes = staticmethod(
        data_info_factory(names=_info_summary_attrs,
                          funcs=[partial(_get_data_attribute, attr=attr)
                                 for attr in _info_summary_attrs]))

    # No nan* methods in numpy < 1.8
    info_summary_stats = staticmethod(
        data_info_factory(names=_stats,
                          funcs=[getattr(np, 'nan' + stat)
                                 for stat in _stats]))

    def __call__(self, option='attributes', out=''):
        """
        Write summary information about data object to the ``out`` filehandle.
        By default this prints to standard output via sys.stdout.

        The ``option`` argument specifies what type of information
        to include.  This can be a string, a function, or a list of
        strings or functions.  Built-in options are:

        - ``attributes``: data object attributes like ``dtype`` and ``format``
        - ``stats``: basic statistics: min, mean, and max

        If a function is specified then that function will be called with the
        data object as its single argument.  The function must return an
        OrderedDict containing the information attributes.

        If a list is provided then the information attributes will be
        appended for each of the options, in order.

        Examples
        --------

        >>> from astropy.table import Column
        >>> c = Column([1, 2], unit='m', dtype='int32')
        >>> c.info()
        dtype = int32
        unit = m
        class = Column
        n_bad = 0
        length = 2

        >>> c.info(['attributes', 'stats'])
        dtype = int32
        unit = m
        class = Column
        mean = 1.5
        std = 0.5
        min = 1
        max = 2
        n_bad = 0
        length = 2

        Parameters
        ----------
        option : str, callable, list of (str or callable)
            Info option, defaults to 'attributes'.
        out : file-like, None
            Output destination, defaults to sys.stdout.  If None then the
            OrderedDict with information attributes is returned

        Returns
        -------
        info : `~collections.OrderedDict` or None
            `~collections.OrderedDict` if out==None else None
        """
        if out == '':
            out = sys.stdout

        dat = self._parent
        info = OrderedDict()
        name = dat.info.name
        if name is not None:
            info['name'] = name

        options = option if isinstance(option, (list, tuple)) else [option]
        for option in options:
            if isinstance(option, str):
                if hasattr(self, 'info_summary_' + option):
                    option = getattr(self, 'info_summary_' + option)
                else:
                    raise ValueError('option={} is not an allowed information type'
                                     .format(option))

            with warnings.catch_warnings():
                for ignore_kwargs in IGNORE_WARNINGS:
                    warnings.filterwarnings('ignore', **ignore_kwargs)
                info.update(option(dat))

        if hasattr(dat, 'mask'):
            n_bad = np.count_nonzero(dat.mask)
        else:
            try:
                n_bad = np.count_nonzero(np.isinf(dat) | np.isnan(dat))
            except Exception:
                n_bad = 0
        info['n_bad'] = n_bad

        try:
            info['length'] = len(dat)
        except (TypeError, IndexError):
            pass

        if out is None:
            return info

        for key, val in info.items():
            if val != '':
                out.write(f'{key} = {val}' + os.linesep)

    def __repr__(self):
        if self._parent is None:
            return super().__repr__()

        out = StringIO()
        self.__call__(out=out)
        return out.getvalue()


class BaseColumnInfo(DataInfo):
    """
    Base info class for anything that can be a column in an astropy
    Table.  There are at least two classes that inherit from this:

      ColumnInfo: for native astropy Column / MaskedColumn objects
      MixinInfo: for mixin column objects

    Note that this class is defined here so that mixins can use it
    without importing the table package.
    """
    attr_names = DataInfo.attr_names.union(['parent_table', 'indices'])
    _attrs_no_copy = set(['parent_table', 'indices'])

    # Context for serialization.  This can be set temporarily via
    # ``serialize_context_as(context)`` context manager to allow downstream
    # code to understand the context in which a column is being serialized.
    # Typical values are 'fits', 'hdf5', 'parquet', 'ecsv', 'yaml'.  Objects
    # like Time or SkyCoord will have different default serialization
    # representations depending on context.
    _serialize_context = None
    __slots__ = ['_format_funcs', '_copy_indices']

    @property
    def parent_table(self):
        value = self._attrs.get('parent_table')
        if callable(value):
            value = value()
        return value

    @parent_table.setter
    def parent_table(self, parent_table):
        if parent_table is None:
            self._attrs.pop('parent_table', None)
        else:
            parent_table = weakref.ref(parent_table)
            self._attrs['parent_table'] = parent_table

    def __init__(self, bound=False):
        super().__init__(bound=bound)

        # If bound to a data object instance then add a _format_funcs dict
        # for caching functions for print formatting.
        if bound:
            self._format_funcs = {}

    def __set__(self, instance, value):
        # For Table columns do not set `info` when the instance is a scalar.
        try:
            if not instance.shape:
                return
        except AttributeError:
            pass

        super().__set__(instance, value)

    def iter_str_vals(self):
        """
        This is a mixin-safe version of Column.iter_str_vals.
        """
        col = self._parent
        if self.parent_table is None:
            from astropy.table.column import FORMATTER as formatter
        else:
            formatter = self.parent_table.formatter

        _pformat_col_iter = formatter._pformat_col_iter
        for str_val in _pformat_col_iter(col, -1, False, False, {}):
            yield str_val

    @property
    def indices(self):
        # Implementation note: the auto-generation as an InfoAttribute cannot
        # be used here, since on access, one should not just return the
        # default (empty list is this case), but set _attrs['indices'] so that
        # if the list is appended to, it is registered here.
        return self._attrs.setdefault('indices', [])

    @indices.setter
    def indices(self, indices):
        self._attrs['indices'] = indices

    def adjust_indices(self, index, value, col_len):
        '''
        Adjust info indices after column modification.

        Parameters
        ----------
        index : slice, int, list, or ndarray
            Element(s) of column to modify. This parameter can
            be a single row number, a list of row numbers, an
            ndarray of row numbers, a boolean ndarray (a mask),
            or a column slice.
        value : int, list, or ndarray
            New value(s) to insert
        col_len : int
            Length of the column
        '''
        if not self.indices:
            return

        if isinstance(index, slice):
            # run through each key in slice
            t = index.indices(col_len)
            keys = list(range(*t))
        elif isinstance(index, np.ndarray) and index.dtype.kind == 'b':
            # boolean mask
            keys = np.where(index)[0]
        else:  # single int
            keys = [index]

        value = np.atleast_1d(value)  # turn array(x) into array([x])
        if value.size == 1:
            # repeat single value
            value = list(value) * len(keys)

        for key, val in zip(keys, value):
            for col_index in self.indices:
                col_index.replace(key, self.name, val)

    def slice_indices(self, col_slice, item, col_len):
        '''
        Given a sliced object, modify its indices
        to correctly represent the slice.

        Parameters
        ----------
        col_slice : `~astropy.table.Column` or mixin
            Sliced object. If not a column, it must be a valid mixin, see
            https://docs.astropy.org/en/stable/table/mixin_columns.html
        item : slice, list, or ndarray
            Slice used to create col_slice
        col_len : int
            Length of original object
        '''
        from astropy.table.sorted_array import SortedArray
        if not getattr(self, '_copy_indices', True):
            # Necessary because MaskedArray will perform a shallow copy
            col_slice.info.indices = []
            return col_slice
        elif isinstance(item, slice):
            col_slice.info.indices = [x[item] for x in self.indices]
        elif self.indices:
            if isinstance(item, np.ndarray) and item.dtype.kind == 'b':
                # boolean mask
                item = np.where(item)[0]
            # Empirical testing suggests that recreating a BST/RBT index is
            # more effective than relabelling when less than ~60% of
            # the total number of rows are involved, and is in general
            # more effective for SortedArray.
            small = len(item) <= 0.6 * col_len
            col_slice.info.indices = []
            for index in self.indices:
                if small or isinstance(index, SortedArray):
                    new_index = index.get_slice(col_slice, item)
                else:
                    new_index = deepcopy(index)
                    new_index.replace_rows(item)
                col_slice.info.indices.append(new_index)

        return col_slice

    @staticmethod
    def merge_cols_attributes(cols, metadata_conflicts, name, attrs):
        """
        Utility method to merge and validate the attributes ``attrs`` for the
        input table columns ``cols``.

        Note that ``dtype`` and ``shape`` attributes are handled specially.
        These should not be passed in ``attrs`` but will always be in the
        returned dict of merged attributes.

        Parameters
        ----------
        cols : list
            List of input Table column objects
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output column name
        attrs : list
            List of attribute names to be merged

        Returns
        -------
        attrs : dict
            Of merged attributes.

        """
        from astropy.table.np_utils import TableMergeError

        def warn_str_func(key, left, right):
            out = ("In merged column '{}' the '{}' attribute does not match "
                   "({} != {}).  Using {} for merged output"
                   .format(name, key, left, right, right))
            return out

        def getattrs(col):
            return {attr: getattr(col.info, attr) for attr in attrs
                    if getattr(col.info, attr, None) is not None}

        out = getattrs(cols[0])
        for col in cols[1:]:
            out = metadata.merge(out, getattrs(col), metadata_conflicts=metadata_conflicts,
                                 warn_str_func=warn_str_func)

        # Output dtype is the superset of all dtypes in in_cols
        out['dtype'] = metadata.common_dtype(cols)

        # Make sure all input shapes are the same
        uniq_shapes = set(col.shape[1:] for col in cols)
        if len(uniq_shapes) != 1:
            raise TableMergeError('columns have different shapes')
        out['shape'] = uniq_shapes.pop()

        # "Merged" output name is the supplied name
        if name is not None:
            out['name'] = name

        return out

    def get_sortable_arrays(self):
        """
        Return a list of arrays which can be lexically sorted to represent
        the order of the parent column.

        The base method raises NotImplementedError and must be overridden.

        Returns
        -------
        arrays : list of ndarray
        """
        raise NotImplementedError(f'column {self.name} is not sortable')


class MixinInfo(BaseColumnInfo):

    @property
    def name(self):
        return self._attrs.get('name')

    @name.setter
    def name(self, name):
        # For mixin columns that live within a table, rename the column in the
        # table when setting the name attribute.  This mirrors the same
        # functionality in the BaseColumn class.
        if self.parent_table is not None:
            new_name = None if name is None else str(name)
            self.parent_table.columns._rename_column(self.name, new_name)

        self._attrs['name'] = name


class ParentDtypeInfo(MixinInfo):
    """Mixin that gets info.dtype from parent"""

    attrs_from_parent = set(['dtype'])  # dtype and unit taken from parent

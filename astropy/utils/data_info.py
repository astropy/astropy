# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Note: these functions and classes are tested extensively in astropy table
# tests via their use in providing mixin column info, and in
# astropy/tests/test_info for providing table and column info summary data.

from __future__ import absolute_import, division, print_function

import os
import sys
import weakref
from copy import deepcopy
import numpy as np
from functools import partial
import warnings
import re

from ..extern import six
from ..utils import OrderedDict
from ..utils.compat import NUMPY_LT_1_8


# Tuple of filterwarnings kwargs to ignore when calling info
IGNORE_WARNINGS = (dict(category=RuntimeWarning,
                        module=r'numpy\.lib\.nanfunctions'),)

STRING_TYPE_NAMES = {(False, 'S'): 'str',  # not PY3
                     (False, 'U'): 'unicode',
                     (True, 'S'): 'bytes', # PY3
                     (True, 'U'): 'str'}

def dtype_info_name(dtype):
    """Return a human-oriented string name of the ``dtype`` arg.
    This can be use by astropy methods that present type information about
    a data object.

    The output is mostly equivalent to ``dtype.name`` which takes the form
    <type_name>[B] where <type_name> is like ``int`` or ``bool`` and [B] is an
    optional number of bits which gets included only for numeric types.

    For bytes, string and unicode types, the output is shown below, where <N>
    is the number of characters.  This representation corresponds to the Python
    type that matches the dtype::

      Numpy          S<N>      U<N>
      Python 2      str<N>  unicode<N>
      Python 3    bytes<N>   str<N>

    Parameters
    ----------
    dtype: str, np.dtype, type
        Input dtype as an object that can be converted via np.dtype()

    Returns
    -------
    dtype_info_name: str
        String name of ``dtype``
    """
    dtype = np.dtype(dtype)
    if dtype.kind in ('S', 'U'):
        length = re.search(r'(\d+)', dtype.str).group(1)
        type_name = STRING_TYPE_NAMES[(six.PY3, dtype.kind)]
        out = type_name + length
    else:
        out = dtype.name

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
    min = 1.0
    median = 2.5
    max = 4.0
    n_bad = 0
    length = 4

    Parameters
    ----------
    names: list
        List of information attribute names
    funcs: list
        List of functions that compute the corresponding information attribute

    Returns
    -------
    func: function
        Function that can be used as a data info option
    """
    def func(dat):
        outs = []
        for name, func in zip(names, funcs):
            try:
                if isinstance(func, six.string_types):
                    out = getattr(dat, func)()
                else:
                    out = func(dat)
            except:
                outs.append('--')
            else:
                outs.append(str(out))

        return OrderedDict(zip(names, outs))
    return func


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


class DataInfo(object):
    """
    Descriptor that data classes use to add an ``info`` attribute for storing
    data attributes in a uniform and portable way.  Note that it *must* be
    called ``info`` so that the DataInfo() object can be stored in the
    ``instance`` using the ``info`` key.  Because owner_cls.x is a descriptor,
    Python doesn't use __dict__['x'] normally, and the descriptor can safely
    store stuff there.  Thanks to http://nbviewer.ipython.org/urls/
    gist.github.com/ChrisBeaumont/5758381/raw/descriptor_writeup.ipynb for
    this trick that works for non-hashable classes.

    Parameters
    ----------
    bound : bool, default=False
        If True this is a descriptor attribute in a class definition, else it
        is a DataInfo() object that is bound to a data object instance.
    """
    _stats = ['mean', 'std', 'min', 'max']
    attrs_from_parent = set()
    attr_names = set(['name', 'unit', 'dtype', 'format', 'description', 'meta'])
    _attrs_no_copy = set()
    _info_summary_attrs = ('dtype', 'shape', 'unit', 'format', 'description', 'class')
    _parent = None

    def __init__(self, bound=False):
        # If bound to a data object instance then create the dict of attributes
        # which stores the info attribute values.
        if bound:
            self._attrs = dict((attr, None) for attr in self.attr_names)

    def __get__(self, instance, owner_cls):
        if instance is None:
            # This is an unbound descriptor on the class
            info = self
        else:
            info = instance.__dict__.get('info')
            if info is None:
                info = instance.__dict__['info'] = self.__class__(bound=True)
            info._parent = instance
        return info

    def __set__(self, instance, value):
        if instance is None:
            # This is an unbound descriptor on the class
            raise ValueError('cannot set unbound descriptor')

        if isinstance(value, DataInfo):
            info = instance.__dict__['info'] = self.__class__(bound=True)
            for attr in info.attr_names - info.attrs_from_parent - info._attrs_no_copy:
                info._attrs[attr] = deepcopy(getattr(value, attr))

        else:
            raise TypeError('info must be set with a DataInfo instance')

    def __getstate__(self):
        return self._attrs

    def __setstate__(self, state):
        self._attrs = state

    def __getattr__(self, attr):
        if attr.startswith('_'):
            return super(DataInfo, self).__getattribute__(attr)

        if attr in self.attrs_from_parent:
            return getattr(self._parent, attr)

        try:
            value = self._attrs[attr]
        except KeyError:
            super(DataInfo, self).__getattribute__(attr)  # Generate AttributeError

        # Weak ref for parent table
        if attr == 'parent_table' and callable(value):
            value = value()

        # Mixins have a default dtype of Object if nothing else was set
        if attr == 'dtype' and value is None:
            value = np.dtype('O')

        return value

    def __setattr__(self, attr, value):
        propobj = getattr(self.__class__, attr, None)

        # If attribute is taken from parent properties and there is not a
        # class property (getter/setter) for this attribute then set
        # attribute directly in parent.
        if attr in self.attrs_from_parent and not isinstance(propobj, property):
            setattr(self._parent, attr, value)
            return

        # Check if there is a property setter and use it if possible.
        if isinstance(propobj, property):
            if propobj.fset is None:
                raise AttributeError("can't set attribute")
            propobj.fset(self, value)
            return

        # Private attr names get directly set
        if attr.startswith('_'):
            super(DataInfo, self).__setattr__(attr, value)
            return

        # Finally this must be an actual data attribute that this class is handling.
        if attr not in self.attr_names:
            raise AttributeError("attribute must be one of {0}".format(self.attr_names))

        if attr == 'parent_table':
            value = None if value is None else weakref.ref(value)

        self._attrs[attr] = value

    info_summary_attributes = staticmethod(
        data_info_factory(names=_info_summary_attrs,
                          funcs=[partial(_get_data_attribute, attr=attr)
                                 for attr in _info_summary_attrs]))

    # No nan* methods in numpy < 1.8
    info_summary_stats = staticmethod(
        data_info_factory(names=_stats,
                          funcs=[getattr(np, ('' if NUMPY_LT_1_8 else 'nan') + stat)
                                 for stat in _stats]))

    def __call__(self, option='attributes', out=''):
        """
        Write summary information about data object to the ``out`` filehandle.
        By default this prints to standard output via sys.stdout.

        The ``option` argument specifies what type of information
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
        option: str, function, list of (str or function)
            Info option (default='attributes')
        out: file-like object, None
            Output destination (default=sys.stdout).  If None then the
            OrderedDict with information attributes is returned

        Returns
        -------
        info: OrderedDict if out==None else None
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
            if isinstance(option, six.string_types):
                if hasattr(self, 'info_summary_' + option):
                    option = getattr(self, 'info_summary_' + option)
                else:
                    raise ValueError('option={0} is not an allowed information type'
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
            except:
                n_bad = 0
        info['n_bad'] = n_bad

        try:
            info['length'] = len(dat)
        except TypeError:
            pass

        if out is None:
            return info

        for key, val in info.items():
            if val != '':
                out.write('{0} = {1}'.format(key, val) + os.linesep)

    def __repr__(self):
        if self._parent is None:
            return super(DataInfo, self).__repr__()

        out = six.moves.cStringIO()
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
    _attrs_no_copy = set(['parent_table'])

    def iter_str_vals(self):
        """
        This is a mixin-safe version of Column.iter_str_vals.
        """
        col = self._parent
        if self.parent_table is None:
            from ..table.column import FORMATTER as formatter
        else:
            formatter = self.parent_table.formatter

        _pformat_col_iter = formatter._pformat_col_iter
        for str_val in _pformat_col_iter(col, -1, False, False, {}):
            yield str_val

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
        else: # single int
            keys = [index]

        value = np.atleast_1d(value) # turn array(x) into array([x])
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
        col_slice : Column or mixin
            Sliced object
        item : slice, list, or ndarray
            Slice used to create col_slice
        col_len : int
            Length of original object
        '''
        from ..table.index import Index
        from ..table.sorted_array import SortedArray
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
            threshold = 0.6
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


class MixinInfo(BaseColumnInfo):
    pass

class ParentDtypeInfo(MixinInfo):
    """Mixin that gets info.dtype from parent"""

    attrs_from_parent = set(['dtype'])  # dtype and unit taken from parent

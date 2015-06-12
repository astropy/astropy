# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Note: these functions and classes are tested extensively in astropy table
# tests via their use in providing mixin column info, and in
# astropy/tests/test_info for providing table and column info summary data.

import os
import sys
from copy import deepcopy
import weakref
import numpy as np
from functools import partial

from ..extern import six
from ..utils import OrderedDict
from ..utils.compat import NUMPY_LT_1_8

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
        val = dat.info.dtype.name
    elif attr == 'shape':
        datshape = dat.shape[1:]
        val = datshape if datshape else ''
    else:
        val = getattr(dat.info, attr)
    if val is None:
        val = ''
    return str(val)


class InfoDescriptor(object):
    """
    Descriptor that data classes use to add an ``info`` attribute for storing
    data attributes in a uniform and portable way.  Note that it *must* be
    called ``info`` so that the ``info_cls()`` object can be stored in the
    ``instance`` using the ``info`` key.  Because owner_cls.x is a descriptor,
    Python doesn't use __dict__['x'] normally, and the descriptor can safely
    store stuff there.  Thanks to http://nbviewer.ipython.org/urls/
    gist.github.com/ChrisBeaumont/5758381/raw/descriptor_writeup.ipynb for
    this trick that works for non-hashable classes.

    Parameters
    ----------
    info_cls : class
        Class reference for the DataInfo subclass that actually stores info
    """
    def __init__(self, info_cls):
        self.info_cls = info_cls

    def __get__(self, instance, owner_cls):
        if 'info' not in instance.__dict__:
            instance.__dict__['info'] = self.info_cls()
        instance.__dict__['info']._parent_ref = weakref.ref(instance)
        return instance.__dict__['info']

    def __set__(self, instance, value):
        if isinstance(value, DataInfo):
            instance.__dict__['info'] = value
        else:
            raise TypeError('info must be set with a DataInfo instance')


class DataInfo(object):

    _stats = ['mean', 'std', 'min', 'max']
    attrs_from_parent = set()
    attr_names = set(['name', 'unit', 'dtype', 'format', 'description',
                      'meta', 'parent_table'])
    _info_summary_attrs = ('dtype', 'shape', 'unit', 'format', 'description', 'class')
    _parent_ref = None

    def __init__(self):
        self._attrs = dict((attr, None) for attr in self.attr_names)

    def __getstate__(self):
        return self._attrs

    def __setstate__(self, state):
        self._attrs = state

    def __getattr__(self, attr):
        if attr.startswith('_'):
            return super(DataInfo, self).__getattribute__(attr)

        if attr in self.attrs_from_parent:
            return getattr(self._parent_ref(), attr)

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
            setattr(self._parent_ref(), attr, value)
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

    def copy(self):
        out = self.__class__()
        for attr in self.attr_names - self.attrs_from_parent:
            setattr(out, attr, deepcopy(getattr(self, attr)))

        return out

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

        dat = self._parent_ref()
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
        out = six.moves.cStringIO()
        self.__call__(out=out)
        return out.getvalue()

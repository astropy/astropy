import os
import sys
from copy import deepcopy
import weakref
import numpy as np
from functools import partial

from ..extern import six
from ..utils import OrderedDict

COLUMN_ATTRS = set(['name', 'unit', 'dtype', 'format', 'description', 'meta', 'parent_table'])
INFO_SUMMARY_ATTRS = ('dtype', 'shape', 'unit', 'format', 'description', 'class')


def column_info_factory(names, funcs):
    """
    Factory to create a function that can be used as an ``option``
    for outputting table or column summary information.

    Examples
    --------
    >>> from astropy.table.info import column_info_factory
    >>> from astropy.table import Column
    >>> c = Column([4., 3., 2., 1.])
    >>> mystats = column_info_factory(names=['min', 'median', 'max'],
    ...                               funcs=[np.min, np.median, np.max])
    >>> c.info(option=mystats)
    min = 1.0
    median = 2.5
    max = 4.0

    Parameters
    ----------
    names: list
        List of information attribute names
    funcs: list
        List of functions that compute the corresponding information attribute

    Returns
    -------
    func: function
        Function that can be used as a table or column info option
    """
    def func(col):
        outs = []
        for name, func in zip(names, funcs):
            try:
                if isinstance(func, six.string_types):
                    out = getattr(col, func)()
                else:
                    out = func(col)
            except:
                outs.append('--')
            else:
                outs.append(str(out))

        return OrderedDict(zip(names, outs))
    return func


def _get_column_attribute(col, attr=None):
    """
    Get a column attribute for the ``attributes`` info summary method
    """
    from ..table.column import Column, MaskedColumn

    if attr == 'class':
        val = '' if type(col) in (Column, MaskedColumn) else type(col).__name__
    elif attr == 'dtype':
        val = col.info.dtype.name
    elif attr == 'shape':
        colshape = col.shape[1:]
        val = colshape if colshape else ''
    else:
        val = getattr(col.info, attr)
    if val is None:
        val = ''
    return str(val)


class DataInfo(object):

    def __init__(self, info_cls):
        self.info_cls = info_cls

    def __get__(self, instance, owner_cls):
        if 'info' not in instance.__dict__:
            instance.__dict__['info'] = self.info_cls()
        instance.__dict__['info']._parent_col = weakref.ref(instance)
        return instance.__dict__['info']

    def __set__(self, instance, value):
        if isinstance(value, BaseDataInfo):
            instance.__dict__['info'] = value
        else:
            raise TypeError('info must be set with a BaseDataInfo instance')


class BaseDataInfo(object):
    _parent_col = None
    # By default the unit and dtype attributes should refer to the parent properties
    _attrs_from_parent = set()
    _stats = ['mean', 'std', 'min', 'max']

    def __init__(self):
        self._attrs = dict((attr, None) for attr in COLUMN_ATTRS)

    def __getstate__(self):
        return (self._attrs, self._attrs_from_parent)

    def __setstate__(self, state):
        self._attrs, self._attrs_from_parent = state

    def __getattr__(self, attr):
        if attr.startswith('_'):
            return super(BaseDataInfo, self).__getattribute__(attr)

        if attr in self._attrs_from_parent:
            return getattr(self._parent_col(), attr)

        try:
            value = self._attrs[attr]
        except KeyError:
            super(BaseDataInfo, self).__getattribute__(attr)  # Generate AttributeError

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
        # ColumnInfo property (getter/setter) for this attribute then set
        # attribute directly in parent.
        if attr in self._attrs_from_parent and not isinstance(propobj, property):
            setattr(self._parent_col(), attr, value)
            return

        # Check if there is a property setter
        if isinstance(propobj, property):
            if propobj.fset is None:
                raise AttributeError("can't set attribute")
            propobj.fset(self, value)
            return

        if attr.startswith('_'):
            super(BaseDataInfo, self).__setattr__(attr, value)
            return

        if attr not in COLUMN_ATTRS:
            raise AttributeError("attribute must be one of {0}".format(COLUMN_ATTRS))

        if attr == 'parent_table':
            value = None if value is None else weakref.ref(value)

        self._attrs[attr] = value

    def copy(self):
        out = self.__class__()
        for attr in COLUMN_ATTRS - self._attrs_from_parent:
            setattr(out, attr, deepcopy(getattr(self, attr)))

        return out

    info_summary_attributes = staticmethod(
        column_info_factory(names=INFO_SUMMARY_ATTRS,
                            funcs=[partial(_get_column_attribute, attr=attr)
                                   for attr in INFO_SUMMARY_ATTRS]))

    info_summary_stats = staticmethod(
        column_info_factory(names=_stats,
                            funcs=[getattr(np, 'nan' + stat) for stat in _stats]))

    def __call__(self, option='attributes', out=''):
        """
        Write summary information about column to the ``out`` filehandle.
        By default this prints to standard output via sys.stdout.

        The ``option` argument specifies what type of information
        to include.  This can be a string, a function, or a list of
        strings or functions.  Built-in options are:

        - ``attributes``: column attributes like ``dtype`` and ``format``
        - ``stats``: basic statistics: min, mean, and max

        If a function is specified then that function will be called with the
        column as its single argument.  The function must return an OrderedDict
        containing the information attributes.

        If a list is provided then the information attributes will be
        appended for each of the options, in order.

        Examples
        --------

        >>> from astropy.table import Column
        >>> c = Column([1, 2, 3], unit='m', dtype='int32')
        >>> c.info()
        dtype = int64
        unit = m
        >>> c.info(['attributes', 'stats'])
        dtype = int32
        unit = m
        min = 1
        mean = 2.0
        max = 3

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

        col = self._parent_col()
        info = OrderedDict()
        name = col.info.name
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
            info.update(option(col))

        if hasattr(col, 'mask'):
            n_bad = np.count_nonzero(col.mask)
        else:
            try:
                n_bad = np.count_nonzero(np.isinf(col) | np.isnan(col))
            except:
                n_bad = 0
        info['n_bad'] = n_bad

        info['length'] = len(col)

        if out is None:
            return info

        for key, val in info.items():
            if val != '':
                out.write('{0} = {1}'.format(key, val) + os.linesep)

    def __repr__(self):
        with six.moves.cStringIO() as out:
            self.__call__(out=out)
            return out.getvalue()

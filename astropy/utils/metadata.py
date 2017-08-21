# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains helper functions and classes for handling metadata.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..extern import six
from ..utils import wraps

import warnings

import collections
from collections import OrderedDict
from copy import deepcopy

import numpy as np
from ..utils.exceptions import AstropyWarning
from ..utils.misc import dtype_bytes_or_chars


__all__ = ['MergeConflictError', 'MergeConflictWarning', 'MERGE_STRATEGIES',
           'common_dtype', 'MergePlus', 'MergeNpConcatenate', 'MergeStrategy',
           'MergeStrategyMeta', 'enable_merge_strategies', 'merge', 'MetaData']


class MergeConflictError(TypeError):
    pass


class MergeConflictWarning(AstropyWarning):
    pass


MERGE_STRATEGIES = []


def common_dtype(arrs):
    """
    Use numpy to find the common dtype for a list of ndarrays.

    Only allow arrays within the following fundamental numpy data types:
    ``np.bool``, ``np.object``, ``np.number``, ``np.character``, ``np.void``

    Parameters
    ----------
    arrs : list of ndarray objects
        Arrays for which to find the common dtype

    Returns
    -------
    dtype_str : str
        String representation of dytpe (dtype ``str`` attribute)
    """
    def dtype(arr):
        return getattr(arr, 'dtype', np.dtype('O'))

    np_types = (np.bool_, np.object_, np.number, np.character, np.void)
    uniq_types = set(tuple(issubclass(dtype(arr).type, np_type) for np_type in np_types)
                     for arr in arrs)
    if len(uniq_types) > 1:
        # Embed into the exception the actual list of incompatible types.
        incompat_types = [dtype(arr).name for arr in arrs]
        tme = MergeConflictError('Arrays have incompatible types {0}'
                                 .format(incompat_types))
        tme._incompat_types = incompat_types
        raise tme

    arrs = [np.empty(1, dtype=dtype(arr)) for arr in arrs]

    # For string-type arrays need to explicitly fill in non-zero
    # values or the final arr_common = .. step is unpredictable.
    for i, arr in enumerate(arrs):
        if arr.dtype.kind in ('S', 'U'):
            arrs[i] = [(u'0' if arr.dtype.kind == 'U' else b'0') *
                       dtype_bytes_or_chars(arr.dtype)]

    arr_common = np.array([arr[0] for arr in arrs])
    return arr_common.dtype.str


class MergeStrategyMeta(type):
    """
    Metaclass that registers MergeStrategy subclasses into the
    MERGE_STRATEGIES registry.
    """

    def __new__(mcls, name, bases, members):
        cls = super(MergeStrategyMeta, mcls).__new__(mcls, name, bases, members)

        # Wrap ``merge`` classmethod to catch any exception and re-raise as
        # MergeConflictError.
        if 'merge' in members and isinstance(members['merge'], classmethod):
            orig_merge = members['merge'].__func__

            @wraps(orig_merge)
            def merge(cls, left, right):
                try:
                    return orig_merge(cls, left, right)
                except Exception as err:
                    raise MergeConflictError(err)

            cls.merge = classmethod(merge)

        # Register merging class (except for base MergeStrategy class)
        if 'types' in members:
            types = members['types']
            if isinstance(types, tuple):
                types = [types]
            for left, right in reversed(types):
                MERGE_STRATEGIES.insert(0, (left, right, cls))

        return cls


@six.add_metaclass(MergeStrategyMeta)
class MergeStrategy(object):
    """
    Base class for defining a strategy for merging metadata from two
    sources, left and right, into a single output.

    The primary functionality for the class is the ``merge(cls, left, right)``
    class method.  This takes ``left`` and ``right`` side arguments and
    returns a single merged output.

    The first class attribute is ``types``.  This is defined as a list of
    (left_types, right_types) tuples that indicate for which input types the
    merge strategy applies.  In determining whether to apply this merge
    strategy to a pair of (left, right) objects, a test is done:
    ``isinstance(left, left_types) and isinstance(right, right_types)``.  For
    example::

      types = [(np.ndarray, np.ndarray),  # Two ndarrays
               (np.ndarray, (list, tuple)),  # ndarray and (list or tuple)
               ((list, tuple), np.ndarray)]  # (list or tuple) and ndarray

    As a convenience, ``types`` can be defined as a single two-tuple instead of
    a list of two-tuples, e.g. ``types = (np.ndarray, np.ndarray)``.

    The other class attribute is ``enabled``, which defaults to ``False`` in
    the base class.  By defining a subclass of ``MergeStrategy`` the new merge
    strategy is automatically registered to be available for use in
    merging. However, by default the new merge strategy is *not enabled*.  This
    prevents inadvertently changing the behavior of unrelated code that is
    performing metadata merge operations.

    In most cases (particularly in library code that others might use) it is
    recommended to leave custom strategies disabled and use the
    `~astropy.utils.metadata.enable_merge_strategies` context manager to locally
    enable the desired strategies.  However, if one is confident that the
    new strategy will not produce unexpected behavior, then one can globally
    enable it by setting the ``enabled`` class attribute to ``True``.

    Examples
    --------
    Here we define a custom merge strategy that takes an int or float on
    the left and right sides and returns a list with the two values.

      >>> from astropy.utils.metadata import MergeStrategy
      >>> class MergeNumbersAsList(MergeStrategy):
      ...     types = ((int, float), (int, float))  # (left_types, right_types)
      ...
      ...     @classmethod
      ...     def merge(cls, left, right):
      ...         return [left, right]

    """
    # Set ``enabled = True`` to globally enable applying this merge strategy.
    # This is not generally recommended.
    enabled = False

    # types = [(left_types, right_types), ...]


class MergePlus(MergeStrategy):
    """
    Merge ``left`` and ``right`` objects using the plus operator.  This
    merge strategy is globally enabled by default.
    """
    types = [(list, list), (tuple, tuple)]
    enabled = True

    @classmethod
    def merge(cls, left, right):
        return left + right


class MergeNpConcatenate(MergeStrategy):
    """
    Merge ``left`` and ``right`` objects using np.concatenate.  This
    merge strategy is globally enabled by default.

    This will upcast a list or tuple to np.ndarray and the output is
    always ndarray.
    """
    types = [(np.ndarray, np.ndarray),
             (np.ndarray, (list, tuple)),
             ((list, tuple), np.ndarray)]
    enabled = True

    @classmethod
    def merge(cls, left, right):
        left, right = np.asanyarray(left), np.asanyarray(right)
        common_dtype([left, right])  # Ensure left and right have compatible dtype
        return np.concatenate([left, right])


def _both_isinstance(left, right, cls):
    return isinstance(left, cls) and isinstance(right, cls)


def _not_equal(left, right):
    try:
        return bool(left != right)
    except Exception:
        return True


class _EnableMergeStrategies(object):
    def __init__(self, *merge_strategies):
        self.merge_strategies = merge_strategies
        self.orig_enabled = {}
        for left_type, right_type, merge_strategy in MERGE_STRATEGIES:
            if issubclass(merge_strategy, merge_strategies):
                self.orig_enabled[merge_strategy] = merge_strategy.enabled
                merge_strategy.enabled = True

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        for merge_strategy, enabled in self.orig_enabled.items():
            merge_strategy.enabled = enabled


def enable_merge_strategies(*merge_strategies):
    """
    Context manager to temporarily enable one or more custom metadata merge
    strategies.

    Examples
    --------
    Here we define a custom merge strategy that takes an int or float on
    the left and right sides and returns a list with the two values.

      >>> from astropy.utils.metadata import MergeStrategy
      >>> class MergeNumbersAsList(MergeStrategy):
      ...     types = ((int, float),  # left side types
      ...              (int, float))  # right side types
      ...     @classmethod
      ...     def merge(cls, left, right):
      ...         return [left, right]

    By defining this class the merge strategy is automatically registered to be
    available for use in merging. However, by default new merge strategies are
    *not enabled*.  This prevents inadvertently changing the behavior of
    unrelated code that is performing metadata merge operations.

    In order to use the new merge strategy, use this context manager as in the
    following example::

      >>> from astropy.table import Table, vstack
      >>> from astropy.utils.metadata import enable_merge_strategies
      >>> t1 = Table([[1]], names=['a'])
      >>> t2 = Table([[2]], names=['a'])
      >>> t1.meta = {'m': 1}
      >>> t2.meta = {'m': 2}
      >>> with enable_merge_strategies(MergeNumbersAsList):
      ...    t12 = vstack([t1, t2])
      >>> t12.meta['m']
      [1, 2]

    One can supply further merge strategies as additional arguments to the
    context manager.

    As a convenience, the enabling operation is actually done by checking
    whether the registered strategies are subclasses of the context manager
    arguments.  This means one can define a related set of merge strategies and
    then enable them all at once by enabling the base class.  As a trivial
    example, *all* registered merge strategies can be enabled with::

      >>> with enable_merge_strategies(MergeStrategy):
      ...    t12 = vstack([t1, t2])

    Parameters
    ----------
    merge_strategies : one or more `~astropy.utils.metadata.MergeStrategy` args
        Merge strategies that will be enabled.

    """

    return _EnableMergeStrategies(*merge_strategies)


def _warn_str_func(key, left, right):
    out = ('Cannot merge meta key {0!r} types {1!r}'
           ' and {2!r}, choosing {0}={3!r}'
           .format(key, type(left), type(right), right))
    return out


def _error_str_func(key, left, right):
    out = ('Cannot merge meta key {0!r} '
           'types {1!r} and {2!r}'
           .format(key, type(left), type(right)))
    return out


def merge(left, right, merge_func=None, metadata_conflicts='warn',
          warn_str_func=_warn_str_func,
          error_str_func=_error_str_func):
    """
    Merge the ``left`` and ``right`` metadata objects.

    This is a simplistic and limited implementation at this point.
    """
    if not _both_isinstance(left, right, dict):
        raise MergeConflictError('Can only merge two dict-based objects')

    out = deepcopy(left)

    for key, val in six.iteritems(right):
        # If no conflict then insert val into out dict and continue
        if key not in out:
            out[key] = deepcopy(val)
            continue

        # There is a conflict that must be resolved
        if _both_isinstance(left[key], right[key], dict):
            out[key] = merge(left[key], right[key], merge_func,
                             metadata_conflicts=metadata_conflicts)

        else:
            try:
                if merge_func is None:
                    for left_type, right_type, merge_cls in MERGE_STRATEGIES:
                        if not merge_cls.enabled:
                            continue
                        if (isinstance(left[key], left_type) and
                                isinstance(right[key], right_type)):
                            out[key] = merge_cls.merge(left[key], right[key])
                            break
                    else:
                        raise MergeConflictError
                else:
                    out[key] = merge_func(left[key], right[key])
            except MergeConflictError:

                # Pick the metadata item that is not None, or they are both not
                # None, then if they are equal, there is no conflict, and if
                # they are different, there is a conflict and we pick the one
                # on the right (or raise an error).

                if left[key] is None:
                    # This may not seem necessary since out[key] gets set to
                    # right[key], but not all objects support != which is
                    # needed for one of the if clauses.
                    out[key] = right[key]
                elif right[key] is None:
                    out[key] = left[key]
                elif _not_equal(left[key], right[key]):
                    if metadata_conflicts == 'warn':
                        warnings.warn(warn_str_func(key, left[key], right[key]),
                                      MergeConflictWarning)
                    elif metadata_conflicts == 'error':
                        raise MergeConflictError(error_str_func(key, left[key], right[key]))
                    elif metadata_conflicts != 'silent':
                        raise ValueError('metadata_conflicts argument must be one '
                                         'of "silent", "warn", or "error"')
                    out[key] = right[key]
                else:
                    out[key] = right[key]

    return out


class MetaData(object):
    """
    A descriptor for classes that have a ``meta`` property.

    This can be set to any valid `~collections.Mapping`.

    Parameters
    ----------
    doc : `str`, optional
        Documentation for the attribute of the class.
        Default is ``""``.

        .. versionadded:: 1.2

    copy : `bool`, optional
        If ``True`` the the value is deepcopied before setting, otherwise it
        is saved as reference.
        Default is ``True``.

        .. versionadded:: 1.2
    """

    def __init__(self, doc="", copy=True):
        self.__doc__ = doc
        self.copy = copy

    def __get__(self, instance, owner):
        if instance is None:
            return self
        if not hasattr(instance, '_meta'):
            instance._meta = OrderedDict()
        return instance._meta

    def __set__(self, instance, value):
        if value is None:
            instance._meta = OrderedDict()
        else:
            if isinstance(value, collections.Mapping):
                if self.copy:
                    instance._meta = deepcopy(value)
                else:
                    instance._meta = value
            else:
                raise TypeError("meta attribute must be dict-like")

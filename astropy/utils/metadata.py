# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains helper functions for handling metadata.
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


class MergeConflictError(TypeError):
    pass


class MergeConflictWarning(AstropyWarning):
    pass


MERGE_STRATEGY_CLASSES = []

def common_dtype(arrs):
    """
    Use numpy to find the common dtype for a list of ndarrays.

    Only allow arrays within the following fundamental numpy data types:
    np.bool_, np.object_, np.number, np.character, np.void

    Parameters
    ----------
    arrs: list of ndarray objects

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
    for arr in arrs:
        if arr.dtype.kind in ('S', 'U'):
            arr[0] = '0' * arr.itemsize

    arr_common = np.array([arr[0] for arr in arrs])
    return arr_common.dtype.str

class MergeStrategyMeta(type):
    """
    Metaclass that registers MergeStrategy subclasses into the
    MERGE_STRATEGY_CLASSES registry.
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
                MERGE_STRATEGY_CLASSES.insert(0, (left, right, cls))

        return cls


@six.add_metaclass(MergeStrategyMeta)
class MergeStrategy(object):
    # Set ``enabled = True`` to globally enable applying this merge strategy
    enabled = False


class MergePlus(MergeStrategy):
    """
    Merge ``left`` and ``right`` objects using the plus operator.  This
    merge strategy is globally enabled by default.
    """
    types = [(list, list),
             (tuple, tuple)]
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
    except:
        return True


def merge(left, right, merge_func=None, metadata_conflicts='warn'):
    """
    Merge the ``left`` and ``right`` metadata objects.

    This is a simplistic and limited implementation at this point.
    """
    if not _both_isinstance(left, right, dict):
        raise MergeConflictError('Can only merge two dict-based objects')

    out = deepcopy(left)

    for key, val in list(six.iteritems(right)):
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
                    for left_type, right_type, merge_cls in MERGE_STRATEGY_CLASSES:
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
                        warnings.warn('Cannot merge meta key {0!r} types {1!r}'
                                      ' and {2!r}, choosing {0}={3!r}'
                                      .format(key, type(left[key]), type(right[key]), right[key]),
                                      MergeConflictWarning)
                    elif metadata_conflicts == 'error':
                        raise MergeConflictError('Cannot merge meta key {0!r} '
                                                 'types {1!r} and {2!r}'
                                                 .format(key, type(left[key]), type(right[key])))
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

    This can be set to any valid mapping.
    """

    def __get__(self, instance, owner):
        if not hasattr(instance, '_meta'):
            instance._meta = OrderedDict()
        return instance._meta

    def __set__(self, instance, value):
        if value is None:
            instance._meta = OrderedDict()
        else:
            if isinstance(value, collections.Mapping):
                instance._meta = deepcopy(value)
            else:
                raise TypeError("meta attribute must be dict-like")

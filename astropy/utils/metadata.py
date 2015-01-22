# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains helper functions for handling metadata.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six

import warnings

import collections
from copy import deepcopy

from .compat.odict import OrderedDict
from ..utils.exceptions import AstropyWarning


class MergeConflictError(TypeError):
    pass


class MergeConflictWarning(AstropyWarning):
    pass


def take_left(left, right):
    return left


def take_right(left, right):
    return left


def concat(left, right):
    """
    Concatenate the ``left`` and ``right`` objects

    Currently only lists and tuples are allowed, but this could change.
    """
    if type(left) != type(right):
        raise MergeConflictError()
    if any(isinstance(left, cls) for cls in (tuple, list)):
        return left + right
    else:
        raise MergeConflictError()


def raise_error(left, right):
    """
    Raise a MergeConflictError if left and right sequences
    are being merged.
    """
    raise MergeConflictError()


def _both_isinstance(left, right, cls):
    return isinstance(left, cls) and isinstance(right, cls)


def merge(left, right, merge_func=concat, metadata_conflicts='warn'):
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
                elif left[key] != right[key]:
                    if metadata_conflicts == 'warn':
                        warnings.warn('Cannot merge meta key {0!r} types {1!r} and {2!r}, choosing {0}={3!r}'
                                                 .format(key, type(left[key]), type(right[key]), right[key]), MergeConflictWarning)
                    elif metadata_conflicts == 'error':
                        raise MergeConflictError('Cannot merge meta key {0!r} types {1!r} and {2!r}'
                                                 .format(key, type(left[key]), type(right[key])))
                    elif metadata_conflicts != 'silent':
                        raise ValueError('metadata_conflicts argument must be one of "silent", "warn", or "error"')
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

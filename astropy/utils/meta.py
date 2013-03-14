# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains helper functions for handling metadata.
"""

from copy import deepcopy
from . import OrderedDict

class MergeConflictError(TypeError):
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


def merge(left, right, merge_func=concat):
    """
    Merge the ``left`` and ``right`` metadata objects.

    This is a simplistic and limited implemenation at this point.
    """
    if not _both_isinstance(left, right, dict): 
        raise MergeConflictError('Can only merge two dict-based objects')

    out = left.__class__()
    
    for key, val in left.items():
        out[key] = deepcopy(val)

    for key, val in right.items():
        # If no conflict then insert val into out dict and continue
        if key not in out:
            out[key] = deepcopy(val)
            continue

        # There is a conflict that must be resolved
        if _both_isinstance(left[key], right[key], dict):
            out[key] = merge(left[key], right[key], merge_func)

        else:
            try:
                out[key] = merge_func(left[key], right[key])
            except MergeConflictError:
                raise MergeConflictError('Cannot merge meta key {0!r} types {1!r} and {2!r}'
                                         .format(key, type(left[key]), type(right[key])))

    return out

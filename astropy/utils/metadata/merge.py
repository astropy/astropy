# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Metadata merging."""

import warnings
from contextlib import contextmanager
from copy import deepcopy
from functools import wraps

import numpy as np

from .exceptions import MergeConflictError, MergeConflictWarning
from .utils import common_dtype

__all__ = [
    "MERGE_STRATEGIES",
    "MergeNpConcatenate",
    "MergePlus",
    "MergeStrategy",
    "enable_merge_strategies",
    "merge",
]

MERGE_STRATEGIES = []


class MergeStrategy:
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

    def __init_subclass__(cls):
        members = vars(cls)
        # Wrap ``merge`` classmethod to catch any exception and re-raise as
        # MergeConflictError.
        if isinstance((merge_ := members.get("merge")), classmethod):
            orig_merge = merge_.__func__

            @wraps(orig_merge)
            def merge(cls, left, right):
                try:
                    return orig_merge(cls, left, right)
                except Exception as err:
                    raise MergeConflictError(err)

            cls.merge = classmethod(merge)

        # Register merging class (except for base MergeStrategy class)
        if (types := members.get("types")) is not None:
            if isinstance(types, tuple):
                types = [types]
            for left, right in reversed(types):
                MERGE_STRATEGIES.insert(0, (left, right, cls))


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

    types = [
        (np.ndarray, np.ndarray),
        (np.ndarray, (list, tuple)),
        ((list, tuple), np.ndarray),
    ]
    enabled = True

    @classmethod
    def merge(cls, left, right):
        left, right = np.asanyarray(left), np.asanyarray(right)
        common_dtype([left, right])  # Ensure left and right have compatible dtype
        return np.concatenate([left, right])


# ============================================================================


@contextmanager
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
    *merge_strategies : :class:`~astropy.utils.metadata.MergeStrategy` class
        Merge strategies that will be enabled.
    """
    orig_enabled = {}
    for _, _, merge_strategy in MERGE_STRATEGIES:
        if issubclass(merge_strategy, merge_strategies):
            orig_enabled[merge_strategy] = merge_strategy.enabled
            merge_strategy.enabled = True

    yield

    for merge_strategy, enabled in orig_enabled.items():
        merge_strategy.enabled = enabled


# =============================================================================


def _both_isinstance(left, right, cls):
    return isinstance(left, cls) and isinstance(right, cls)


def _not_equal(left, right):
    try:
        return bool(left != right)
    except Exception:
        return True


def _warn_str_func(key, left, right):
    out = (
        f"Cannot merge meta key {key!r} types {type(left)!r}"
        f" and {type(right)!r}, choosing {key}={right!r}"
    )
    return out


def _error_str_func(key, left, right):
    out = f"Cannot merge meta key {key!r} types {type(left)!r} and {type(right)!r}"
    return out


def _raise_merge_conflict_error(left, right):
    raise MergeConflictError


def merge(
    left,
    right,
    merge_func=None,
    metadata_conflicts="warn",
    warn_str_func=_warn_str_func,
    error_str_func=_error_str_func,
):
    """
    Merge the ``left`` and ``right`` metadata objects.

    This is a simplistic and limited implementation at this point.
    """
    if not _both_isinstance(left, right, dict):
        raise MergeConflictError("Can only merge two dict-based objects")

    out = deepcopy(left)

    for key, right_val in right.items():
        # If no conflict then insert val into out dict and continue
        if key not in out:
            out[key] = deepcopy(right_val)
            continue

        # There is a conflict that must be resolved
        left_val = out[key]
        if _both_isinstance(left_val, right_val, dict):
            out[key] = merge(
                left_val, right_val, merge_func, metadata_conflicts=metadata_conflicts
            )
            continue

        local_merge_func = merge_func
        if local_merge_func is None:
            for left_type, right_type, merge_cls in MERGE_STRATEGIES:
                if (
                    merge_cls.enabled
                    and isinstance(left_val, left_type)
                    and isinstance(right_val, right_type)
                ):
                    local_merge_func = merge_cls.merge
                    break
            else:
                local_merge_func = _raise_merge_conflict_error
        try:
            out[key] = local_merge_func(left_val, right_val)
        except MergeConflictError:
            # If right_val is None or equal to left_val then we are happy with what's
            # already in out. Otherwise right_val takes priority (or there's an error).
            if right_val is not None and _not_equal(left_val, right_val):
                if metadata_conflicts == "warn":
                    warnings.warn(
                        warn_str_func(key, left_val, right_val), MergeConflictWarning
                    )
                elif metadata_conflicts == "error":
                    raise MergeConflictError(error_str_func(key, left_val, right_val))
                elif metadata_conflicts != "silent":
                    raise ValueError(
                        "metadata_conflicts argument must be one "
                        'of "silent", "warn", or "error"'
                    )
                out[key] = right_val
    return out

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module containing specialized collection classes.
"""

from __future__ import annotations

import warnings
from functools import wraps
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Any, Self, SupportsIndex


class HomogeneousList(list):
    """
    A subclass of list that contains only elements of a given type or
    types.  If an item that is not of the specified type is added to
    the list, a `TypeError` is raised.
    """

    def __init__(self, types, values=[]):
        """
        Parameters
        ----------
        types : sequence of types
            The types to accept.

        values : sequence, optional
            An initial set of values.
        """
        self._types = types
        super().__init__()
        for v in values:
            self.append(v)

    def _assert(self, x):
        if not isinstance(x, self._types):
            raise TypeError(
                f"homogeneous list must contain only objects of type '{self._types}'"
            )

    def __iadd__(self, other):
        self.extend(other)
        return self

    def __setitem__(self, idx, value):
        if isinstance(idx, slice):
            value = list(value)
            for item in value:
                self._assert(item)
        else:
            self._assert(value)
        return super().__setitem__(idx, value)

    def append(self, x):
        self._assert(x)
        return super().append(x)

    def insert(self, i, x):
        self._assert(x)
        return super().insert(i, x)

    def extend(self, x):
        for item in x:
            self._assert(item)
            super().append(item)


def _locked_mutator(instance_method):
    @wraps(instance_method)
    def owned_method(self, *args: Any, _owned: bool = False, **kwargs):
        if _owned:
            return instance_method(self, *args, **kwargs)

        cls_name = self._owner_class.__name__
        mtd_name = instance_method.__name__
        msg = (
            f"Direct mutations of {cls_name}.{self._attr_name} "
            f"via {mtd_name} are deprecated since astropy 7.1 and "
            "will raise an exception in the future. "
        )

        if (repl := self._replacements.get(mtd_name)) is not None:
            msg += f"Please use {cls_name}.{repl} instead"
        else:
            msg += (
                "There is currently no planned replacement for this method. "
                "If you need one, please open a feature request at "
                "https://github.com/astropy/astropy/issues/new/choose"
            )

        warnings.warn(
            msg,
            category=DeprecationWarning,
            stacklevel=2,
        )
        return instance_method(self, *args, **kwargs)

    return owned_method


class _OwnedMixin:
    """
    A collection that can only be mutated by some owner object.
    This is meant as a helper class to prevent users from bypassing dedicated
    APIs by overriding mutating methods from the list API.

    As of astropy 7.1, direct mutations are considered deprecated so a warning
    is emitted. In a future major release of astropy, this class is meant to be
    removed, and its instances should be made completely private, possibly
    re-exposed as immutable properties (tuple).

    Within multiple inheritance, this class is assumed to appear first on the
    inheritance chain, and be combined with (a subclass of) list e.g.

    >>> class _OwnedHomogeneousList(_OwnedMixin, HomogeneousList):
    ...     pass

    otherwise, deprecation messages' might not point to the correct line.

    """

    def __init__(
        self,
        *args: Any,
        owner_class: type,
        public_attr_name: str,
        replacements: dict[str, str],
        **kwargs,
    ):
        # storing metadata to make warnings and errors more useful
        self._owner_class = owner_class
        self._attr_name = public_attr_name
        self._replacements = replacements
        super().__init__(*args, **kwargs)

    @_locked_mutator
    def append(self, object, /):
        return super().append(object)

    @_locked_mutator
    def clear(self, /):
        return super().clear()

    @_locked_mutator
    def extend(self, iterable, /):
        return super().extend(iterable)

    @_locked_mutator
    def insert(self, index, object, /):
        return super().insert(index, object)

    @_locked_mutator
    def pop(self, index=-1, /):
        return super().pop(index)

    @_locked_mutator
    def remove(self, value, /):
        return super().remove(value)

    @_locked_mutator
    def reverse(self, /):
        return super().reverse()

    @_locked_mutator
    def sort(self, /, *, key=None, reverse=False):
        return super().sort(key=key, reverse=reverse)

    @_locked_mutator
    def __iadd__(self, value: Iterable) -> Self:
        return super().__iadd__(value)

    @_locked_mutator
    def __imul__(self, value: SupportsIndex) -> Self:
        return super().__imul__(value)

    @_locked_mutator
    def __setitem__(self, key, value, /) -> Self:
        return super().__setitem__(key, value)

    @_locked_mutator
    def __delitem__(self, key: SupportsIndex | slice) -> None:
        return super().__delitem__(key)


class _OwnedHomogeneousList(_OwnedMixin, HomogeneousList):
    pass

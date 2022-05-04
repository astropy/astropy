# -*- coding: utf-8 -*-

"""Descriptors."""

# STDLIB
from copy import deepcopy
import weakref
from typing import Any, Optional, Type, TypeVar

__all__ = ["InstanceDescriptor", "SlotsInstanceDescriptor"]


_EnclType = TypeVar("_EnclType")
_SIDT = TypeVar("_SIDT", bound="SlotsInstanceDescriptor")


class SlotsInstanceDescriptor:
    """Descriptor that provides access to its enclosing instance."""

    __slots__ = ("__weakref__", "_parent_attr", "_parent_ref")
    # __weakref__ so a weakref can be made to a SlotsInstanceDescriptor

    _parent_attr: str
    """Return the name of the attribute on the enclosing object"""

    _parent_ref: Optional[weakref.ReferenceType]
    """Weak reference to the parent object."""

    _parent_ref_missing_msg = "no reference exists to the original enclosing object"

    _get_deepcopy: bool = False
    """
    Whether ``__get__`` uses `copy.deepcopy` or the normal class constructor
    to make descriptor instances.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        # Depending on the placement in an inheritance chain __init__
        # may or may not take arguments. Both options are supported.
        try:
            super().__init__(*args, **kwargs)
        except TypeError:
            super().__init__()

        self._parent_ref = None

    @property
    def _parent(self) -> _EnclType:
        """Enclosing instance."""
        if isinstance(self._parent_ref, weakref.ReferenceType):
            enclosing: Optional[_EnclType] = self._parent_ref()
        else:
            enclosing = None

        if enclosing is None:
            raise AttributeError(self._parent_ref_missing_msg)

        return enclosing

    def __set_name__(self, objcls: Any, name: str) -> None:
        # Depending on the placement in an inheritance chain __set_name__
        # may or may not exist. Both options are supported.
        try:
            super().__set_name__(objcls, name)
        except AttributeError:
            pass

        # The name of the attribute on the enclosing object
        self._parent_attr = name

    def __get__(self: _SIDT, obj: Optional[_EnclType], objcls: Optional[Type[_EnclType]], **kwargs: Any) -> _SIDT:
        # When called without an instance, return self to allow access
        # to descriptor attributes.
        if obj is None:
            return self

        # accessed from an obj
        descriptor: Optional[_SIDT] = obj.__dict__.get(self._parent_attr)  # get from obj
        if descriptor is None:  # hasn't been created on the obj
            # Make a new ``self``-like descriptor
            descriptor = type(self)(**kwargs) if not self._get_deepcopy else deepcopy(self)
            # Copy over attributes from ``__set_name__``.
            # If a subclass has others, either set them in that ``__get__``
            # or use the ``deepcopy`` option.
            descriptor._parent_attr = self._parent_attr
            # put descriptor in enclosing instance's dictionary
            obj.__dict__[self._parent_attr] = descriptor

        # We set `_parent_ref` on every call, since if one makes copies of objs,
        # 'descriptor' will be copied as well, which will lose the reference.
        descriptor._parent_ref = weakref.ref(obj)

        return descriptor


class InstanceDescriptor(SlotsInstanceDescriptor):
    __slots__ = ("__dict__", "_parent_attr", "_parent_ref")

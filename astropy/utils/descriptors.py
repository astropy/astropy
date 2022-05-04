# -*- coding: utf-8 -*-

"""Descriptors."""

# STDLIB
import weakref
from typing import Any, Optional, Type, TypeVar

__all__ = ["InstanceDescriptor", "SlotsInstanceDescriptor"]


_EnclType = TypeVar("_EnclType")
_SIDT = TypeVar("_SIDT", bound="SlotsInstanceDescriptor")


class SlotsInstanceDescriptor:
    """Descriptor that provides access to its parent instance."""

    __slots__ = ("_parent_attr", "_parent_ref")

    _parent_attr: str
    _parent_ref: weakref.ReferenceType

    def __init__(self, **kwargs) -> None:
        super().__init__()

    @property
    def _parent(self) -> _EnclType:
        """Parent instance."""
        if isinstance(getattr(self, "_parent_ref", None), weakref.ReferenceType):
            parent = self._parent_ref()
        else:
            parent = None

        if parent is None:
            raise ValueError("no reference exists to the original parent object")

        return parent

    def __set_name__(self, _: Any, name: str) -> None:
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
            descriptor = type(self)(**kwargs)
            # Copy over attributes from ``__set_name__``.
            # If a subclass has others, either set them in that ``__get__``
            # or use the ``deepcopy`` option.
            descriptor._parent_attr = self._parent_attr
            # put descriptor in parent instance's dictionary
            obj.__dict__[self._parent_attr] = descriptor

        # We set `_parent_ref` on every call, since if one makes copies of objs,
        # 'descriptor' will be copied as well, which will lose the reference.
        descriptor._parent_ref = weakref.ref(obj)  # type: ignore

        return descriptor


class InstanceDescriptor(SlotsInstanceDescriptor):
    __slots__ = ("__dict__", "_parent_attr", "_parent_ref")

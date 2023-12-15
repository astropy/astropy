"""
A simple class to manage a piece of global science state.  See
:ref:`astropy:config-developer` for more details.
"""


__all__ = ["replace"]


import dataclasses
from functools import singledispatch
from typing import Any, ClassVar, Protocol, TypeVar, runtime_checkable

T = TypeVar("T")


@singledispatch
def replace(obj: T, /, **kwargs: Any) -> T:
    """Replace attributes on an object with new values.

    Parameters
    ----------
    obj : object, positional-only
        The object to replace attributes on.
    **kwargs : dict
        The attributes to replace and their new values.

    Returns
    -------
    obj : object
        The object with the replaced attributes.
    """
    msg = f"cannot replace attributes on {obj!r}"
    raise TypeError(msg)


@runtime_checkable
class _DataclassInstance(Protocol):
    """Runtime protocol for dataclass instances."""

    __dataclass_fields__: ClassVar[dict[str, dataclasses.Field[Any]]]


@replace.register
def _replace_dataclass(obj: _DataclassInstance, /, **kwargs: Any) -> _DataclassInstance:  # type: ignore[misc]
    return dataclasses.replace(obj, **kwargs)

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Static typing for :mod:`astropy.cosmology`."""

from __future__ import annotations

__all__ = ["PathLike", "ReadableFileLike", "WriteableFileLike"]

import os
from typing import TYPE_CHECKING, Protocol, TypeVar, runtime_checkable

if TYPE_CHECKING:
    from typing import TypeAlias

T_co = TypeVar("T_co", covariant=True)
T_contra = TypeVar("T_contra", contravariant=True)

PathLike: TypeAlias = str | bytes | os.PathLike
"""Type alias for a path-like object.

This is a union of :class:`str`, :class:`bytes`, and :class:`~os.PathLike`.
"""


@runtime_checkable
class ReadableFileLike(Protocol[T_co]):
    """A file-like object that supports reading with a method ``read``.

    This is a :class:`~typing.Protocol` that can be used to annotate file-like
    objects. It is also runtime-checkable and can be used with :func:`isinstance`.
    See :func:`~typing.runtime_checkable` for more information about how runtime
    checking with Protocols works.
    """

    def read(self) -> T_co:
        ...


@runtime_checkable
class WriteableFileLike(Protocol[T_contra]):
    """Protocol for file-like objects that are writeable.

    This is a :class:`~typing.Protocol` that can be used to annotate file-like
    objects. It is also runtime-checkable and can be used with :func:`isinstance`.
    See :func:`~typing.runtime_checkable` for more information about how runtime
    checking with Protocols works.
    """

    def write(self, data: T_contra) -> None:
        ...

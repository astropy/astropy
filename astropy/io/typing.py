# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Type annotations for ``astropy.io``.

These are type annotations for I/O-related functions and classes. Some of the type
objects can also be used as runtime-checkable :class:`~typing.Protocol` objects.
"""

from __future__ import annotations

__all__ = ["PathLike", "ReadableFileLike", "WriteableFileLike"]

import os
from typing import TYPE_CHECKING, Protocol, TypeVar, runtime_checkable

if TYPE_CHECKING:
    from typing import TypeAlias

_T_co = TypeVar("_T_co", covariant=True)
_T_contra = TypeVar("_T_contra", contravariant=True)

PathLike: TypeAlias = str | bytes | os.PathLike
"""Type alias for a path-like object.

This is a union of :class:`str`, :class:`bytes`, and :class:`~os.PathLike`.
"""


@runtime_checkable
class ReadableFileLike(Protocol[_T_co]):
    """A file-like object that supports reading with a method ``read``.

    This is a :class:`~typing.Protocol` that can be used to annotate file-like
    objects. It is also runtime-checkable and can be used with :func:`isinstance`.
    See :func:`~typing.runtime_checkable` for more information about how runtime
    checking with Protocols works.
    """

    def read(self) -> _T_co: ...


@runtime_checkable
class WriteableFileLike(Protocol[_T_contra]):
    """A file-like object that supports writing with a method ``write``.

    This is a :class:`~typing.Protocol` that can be used to annotate file-like
    objects. It is also runtime-checkable and can be used with :func:`isinstance`.
    See :func:`~typing.runtime_checkable` for more information about how runtime
    checking with Protocols works.
    """

    def write(self, data: _T_contra) -> None: ...

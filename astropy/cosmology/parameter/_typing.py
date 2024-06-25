"""Type hints for the parameter module."""

__all__: list[str] = []

from collections.abc import Callable
from typing import Any, TypeAlias, TypeVar

from astropy.cosmology.core import Cosmology

from ._core import Parameter

VT = TypeVar("VT")
_FValidateCallable: TypeAlias = Callable[[Cosmology, Parameter, Any], VT]  # noqa: PYI047

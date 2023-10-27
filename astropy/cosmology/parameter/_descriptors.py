# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__ = []

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, NoReturn

from astropy.utils.compat.misc import PYTHON_LT_3_10

if TYPE_CHECKING:
    from astropy.cosmology.core import Cosmology


_dataclass_kwargs = {} if PYTHON_LT_3_10 else {"slots": True}


@dataclass(frozen=True, **_dataclass_kwargs)
class ParametersAttribute:
    """Immutable mapping of the Parameters.

    If accessed from the class, this returns a mapping of the Parameter
    objects themselves.  If accessed from an instance, this returns a
    mapping of the values of the Parameters.
    """

    attr_name: str = field(init=False)

    def __set_name__(self, owner: type[Cosmology], name: str) -> None:
        object.__setattr__(self, "attr_name", f"_{name}")

    def __get__(
        self, instance: Cosmology | None, owner: type[Cosmology] | None
    ) -> MappingProxyType[str, Any]:
        # called from the class
        if instance is None:
            return getattr(owner, self.attr_name)
        # called from the instance
        return MappingProxyType(
            {n: getattr(instance, n) for n in getattr(instance, self.attr_name)}
        )

    def __set__(self, instance: Any, value: Any) -> NoReturn:
        msg = f"cannot set {self.attr_name.lstrip('_').rstrip('_')!r} of {instance!r}."
        raise AttributeError(msg)

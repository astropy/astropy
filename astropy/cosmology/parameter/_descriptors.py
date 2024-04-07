# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__: list[str] = []

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, NoReturn

if TYPE_CHECKING:
    from astropy.cosmology.core import Cosmology


@dataclass(frozen=True, slots=True)
class ParametersAttribute:
    """Immutable mapping of the :class:`~astropy.cosmology.Parameter` objects or values.

    If accessed from the :class:`~astropy.cosmology.Cosmology` class, this returns a
    mapping of the :class:`~astropy.cosmology.Parameter` objects themselves.  If
    accessed from an instance, this returns a mapping of the values of the Parameters.

    This class is used to implement :obj:`astropy.cosmology.Cosmology.parameters`.

    Parameters
    ----------
    attr_name : str
        The name of the class attribute that is a `~types.MappingProxyType[str,
        astropy.cosmology.Parameter]` of all the cosmology's parameters. When accessed
        from the class, this attribute is returned. When accessed from an instance, a
        mapping of the cosmology instance's values for each key is returned.

    Examples
    --------
    The normal usage of this class is the ``parameters`` attribute of
    :class:`~astropy.cosmology.Cosmology`.

        >>> from astropy.cosmology import FlatLambdaCDM, Planck18

        >>> FlatLambdaCDM.parameters
        mappingproxy({'H0': Parameter(...), ...})

        >>> Planck18.parameters
        mappingproxy({'H0': <Quantity 67.66 km / (Mpc s)>, ...})
    """

    attr_name: str
    """Class attribute name on Cosmology for the mapping of Parameter objects."""

    _name: str = field(init=False)
    """The name of the descriptor on the containing class."""

    def __set_name__(self, owner: Any, name: str) -> None:
        object.__setattr__(self, "_name", name)

    def __get__(
        self, instance: Cosmology | None, owner: type[Cosmology] | None
    ) -> MappingProxyType[str, Any]:
        # Called from the class
        if instance is None:
            return getattr(owner, self.attr_name)
        # Called from the instance
        return MappingProxyType(
            {n: getattr(instance, n) for n in getattr(instance, self.attr_name)}
        )

    def __set__(self, instance: Any, value: Any) -> NoReturn:
        msg = f"cannot set {self._name!r} of {instance!r}."
        raise AttributeError(msg)

# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__: list[str] = []

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, NoReturn

from astropy.utils.compat.misc import PYTHON_LT_3_10

if TYPE_CHECKING:
    from astropy.cosmology.core import Cosmology


_dataclass_kwargs = {} if PYTHON_LT_3_10 else {"slots": True}


@dataclass(frozen=True, **_dataclass_kwargs)
class ParametersAttribute:
    """Immutable mapping of the :class:`~astropy.cosmology.Parameter` objects or values.

    If accessed from the :class:`~astropy.cosmology.Cosmology` class, this returns a
    mapping of the :class:`~astropy.cosmology.Parameter` objects themselves.  If
    accessed from an instance, this returns a mapping of the values of the Parameters.

    This class is used to implement the ``parameters`` attribute of
    :class:`~astropy.cosmology.Cosmology`.

    Parameters
    ----------
    attr_name : str
        The name of the class attribute containing the Parameter objects.

    Examples
    --------
    This class works as a descriptor on the containing class. When accessed from the
    class, it returns the attribute specified by ``attr_name``. When accessed from an
    instance, it returns a proxy mapping of the value for each key in ``attr_name``.

    In the following example ``attr_name`` is set to ``_attr_map``, which is not a
    dictionary, but a tuple of strings. When accessed from the class, the tuple is
    returned. When accessed from an instance, a mapping of the values of the attributes
    named in the tuple is returned.

    >>> class Obj:
    ...     a = 1
    ...     b = 2
    ...     c = 3
    ...     _attr_map = ("a", "b", "c")
    ...     attr = ParametersAttribute(attr_name="_attr_map")

    >>> Obj.attr
    ("a", "b", "c")

    >>> obj = Obj()
    >>> obj.attr
    mappingproxy({'a': 1, 'b': 2, 'c': 3})
    """

    attr_name: str
    """The name of the class attribute containing the Parameter objects.

    ``attr_name`` is different from the name of the descriptor on the containing class,
    which is set by :meth:`__set_name__`.
    """

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

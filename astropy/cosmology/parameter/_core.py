# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__ = ["Parameter"]

import copy
from collections.abc import Callable
from dataclasses import KW_ONLY, dataclass, field, fields, is_dataclass, replace
from enum import Enum, auto
from typing import TYPE_CHECKING, Any, Generic, TypeAlias, TypeVar, overload

import astropy.units as u

from ._converter import _REGISTRY_FVALIDATORS, _register_validator

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self

    from astropy.cosmology import Cosmology


_VT = TypeVar("_VT")  # the type of the VParameter value
_FValidateCallable: TypeAlias = Callable[["Cosmology", "Parameter", Any], _VT]


class Sentinel(Enum):
    """Sentinel values for Parameter fields."""

    MISSING = auto()
    """A sentinel value signifying a missing default."""

    def __repr__(self) -> str:
        return f"<{self.name}>"


MISSING = Sentinel.MISSING


@dataclass(frozen=True)
class _UnitField:
    # TODO: rm this class when py3.13+ allows for `field(converter=...)`

    def __get__(
        self, obj: Parameter | None, objcls: type[Parameter] | None
    ) -> u.Unit | None:
        if obj is None:  # calling `Parameter.unit` from the class
            return None
        return getattr(obj, "_unit", None)

    def __set__(self, obj: Parameter, value: Any) -> None:
        object.__setattr__(obj, "_unit", u.Unit(value) if value is not None else None)


@dataclass(frozen=True)
class _FValidateField(Generic[_VT]):
    default: _FValidateCallable[_VT] | str = "default"

    def __get__(
        self, obj: Parameter | None, objcls: type[Parameter] | None
    ) -> _FValidateCallable[_VT] | str:
        if obj is None:  # calling `Parameter.fvalidate` from the class
            return self.default
        return obj._fvalidate  # calling `Parameter.fvalidate` from an instance

    def __set__(self, obj: Parameter, value: Any) -> None:
        # Always store input fvalidate.
        object.__setattr__(obj, "_fvalidate_in", value)

        # Process to the callable.
        if value in _REGISTRY_FVALIDATORS:
            value = _REGISTRY_FVALIDATORS[value]
        elif isinstance(value, str):
            msg = f"`fvalidate`, if str, must be in {_REGISTRY_FVALIDATORS.keys()}"
            raise ValueError(msg)
        elif not callable(value):
            msg = f"`fvalidate` must be a function or {_REGISTRY_FVALIDATORS.keys()}"
            raise TypeError(msg)
        object.__setattr__(obj, "_fvalidate", value)


@dataclass(frozen=True)
class Parameter(Generic[_VT]):
    r"""Cosmological parameter (descriptor).

    Should only be used with a :class:`~astropy.cosmology.Cosmology` subclass.

    Parameters
    ----------
    default : Any (optional, keyword-only)
        Default value of the Parameter. If not given the
        Parameter must be set when initializing the cosmology.
    derived : bool (optional, keyword-only)
        Whether the Parameter is 'derived', default `False`.
        Derived parameters behave similarly to normal parameters, but are not
        sorted by the |Cosmology| signature (probably not there) and are not
        included in all methods. For reference, see ``Ode0`` in
        ``FlatFLRWMixin``, which removes :math:`\Omega_{de,0}`` as an
        independent parameter (:math:`\Omega_{de,0} \equiv 1 - \Omega_{tot}`).
    unit : unit-like or None (optional, keyword-only)
        The `~astropy.units.Unit` for the Parameter. If None (default) no
        unit as assumed.
    equivalencies : `~astropy.units.Equivalency` or sequence thereof
        Unit equivalencies for this Parameter.
    fvalidate : callable[[object, object, Any], Any] or str (optional, keyword-only)
        Function to validate the Parameter value from instances of the
        cosmology class. If "default", uses default validator to assign units
        (with equivalencies), if Parameter has units.
        For other valid string options, see ``Parameter._registry_validators``.
        'fvalidate' can also be set through a decorator with
        :meth:`~astropy.cosmology.Parameter.validator`.
    doc : str or None (optional, keyword-only)
        Parameter description.

    Examples
    --------
    For worked examples see :class:`~astropy.cosmology.FLRW`.
    """

    _: KW_ONLY

    default: Any = MISSING
    """Default value of the Parameter.

    By default set to ``MISSING``, which indicates the parameter must be set
    when initializing the cosmology.
    """

    derived: bool = False
    """Whether the Parameter can be set, or is derived, on the cosmology."""

    # Units
    unit: _UnitField = _UnitField()
    """The unit of the Parameter (can be `None` for unitless)."""

    equivalencies: u.Equivalency | Sequence[u.Equivalency] = field(default_factory=list)
    """Unit equivalencies available when setting the parameter."""

    # Setting
    fvalidate: _FValidateField[_VT] = _FValidateField(default="default")
    """Function to validate/convert values when setting the Parameter."""

    # Info
    doc: str | None = None
    """Parameter description."""

    name: str = field(init=False, compare=True, default=None, repr=False)
    """The name of the Parameter on the Cosmology.

    Cannot be set directly.
    """

    def __post_init__(self) -> None:
        self._fvalidate_in: _FValidateCallable[_VT] | str
        self._fvalidate: _FValidateCallable[_VT]
        object.__setattr__(self, "__doc__", self.doc)
        # Now setting a dummy attribute name. The cosmology class will call
        # `__set_name__`, passing the real attribute name. However, if Parameter is not
        # init'ed as a descriptor then this ensures that all declared fields exist.
        self.__set_name__(None, "name not initialized")

    def __set_name__(self, cosmo_cls: type, name: str) -> None:
        # attribute name on container cosmology class
        object.__setattr__(self, "name", name)

    # -------------------------------------------
    # descriptor and property-like methods

    @overload
    def __get__(self, cosmology: None, cosmo_cls: Any) -> Parameter: ...

    @overload
    def __get__(self, cosmology: Cosmology, cosmo_cls: Any) -> _VT: ...

    def __get__(
        self, cosmology: Cosmology | None, cosmo_cls: Any = None
    ) -> Parameter | _VT:
        # Get from class
        if cosmology is None:
            # If the Parameter is being set as part of a dataclass constructor, then we
            # raise an AttributeError if the default is MISSING. This is to prevent the
            # Parameter from being set as the default value of the dataclass field and
            # erroneously included in the class' __init__ signature.
            if self.default is MISSING and (
                not is_dataclass(cosmo_cls)
                or self.name not in cosmo_cls.__dataclass_fields__
            ):
                raise AttributeError
            return self

        # Get from instance
        return cosmology.__dict__[self.name]

    def __set__(self, cosmology: Cosmology, value: Any) -> None:
        """Allows attribute setting once.

        Raises AttributeError subsequently.
        """
        # Raise error if setting 2nd time. The built-in Cosmology objects are frozen
        # dataclasses and this is redundant, however user defined cosmology classes do
        # not have to be frozen.
        if self.name in cosmology.__dict__:
            raise AttributeError(f"cannot assign to field {self.name!r}")

        # Change `self` to the default value if default is MISSING.
        # This is done for backwards compatibility only - so that Parameter can be used
        # in a dataclass and still return `self` when accessed from a class.
        # Accessing the Parameter object via `cosmo_cls.param_name` will be removed
        # in favor of `cosmo_cls.parameters["param_name"]`.
        if value is self:
            value = self.default

        # Validate value, generally setting units if present
        value = self.validate(cosmology, copy.deepcopy(value))

        # Make the value read-only, if ndarray-like
        if hasattr(value, "setflags"):
            value.setflags(write=False)

        # Set the value on the cosmology
        cosmology.__dict__[self.name] = value

    # -------------------------------------------
    # validate value

    def validator(self, fvalidate: _FValidateCallable[_VT]) -> Self:
        """Make new Parameter with custom ``fvalidate``.

        Note: ``Parameter.fvalidator`` must be the top-most descriptor decorator.

        Parameters
        ----------
        fvalidate : callable[[type, type, Any], Any]

        Returns
        -------
        `~astropy.cosmology.Parameter`
            Copy of this Parameter but with custom ``fvalidate``.
        """
        return self.clone(fvalidate=fvalidate)

    def validate(self, cosmology: Cosmology, value: Any) -> _VT:
        """Run the validator on this Parameter.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.Cosmology` instance
        value : Any
            The object to validate.

        Returns
        -------
        Any
            The output of calling ``fvalidate(cosmology, self, value)``
            (yes, that parameter order).
        """
        return self._fvalidate(cosmology, self, value)

    @overload
    @staticmethod
    def register_validator(
        key: str, fvalidate: _FValidateCallable[_VT]
    ) -> _FValidateCallable[_VT]: ...

    @overload
    @staticmethod
    def register_validator(
        key: str, fvalidate: None = None
    ) -> Callable[[_FValidateCallable[_VT]], _FValidateCallable[_VT]]: ...

    @staticmethod
    def register_validator(
        key: str, fvalidate: _FValidateCallable[_VT] | None = None
    ) -> (
        _FValidateCallable[_VT]
        | Callable[[_FValidateCallable[_VT]], _FValidateCallable[_VT]]
    ):
        """Decorator to register a new kind of validator function.

        Parameters
        ----------
        key : str
        fvalidate : callable[[object, object, Any], Any] or None, optional
            Value validation function.

        Returns
        -------
        ``validator`` or callable[``validator``]
            if validator is None returns a function that takes and registers a
            validator. This allows ``register_validator`` to be used as a
            decorator.
        """
        return _register_validator(key, fvalidate=fvalidate)

    # -------------------------------------------

    def clone(self, **kw: Any) -> Self:
        """Clone this `Parameter`, changing any constructor argument.

        Parameters
        ----------
        **kw : Any
            Passed to constructor. The current values, eg. ``fvalidate`` are
            used as the default values, so an empty ``**kw`` is an exact copy.

        Examples
        --------
        >>> p = Parameter()
        >>> p
        Parameter(derived=False, unit=None, equivalencies=[],
                  fvalidate='default', doc=None)

        >>> p.clone(unit="km")
        Parameter(derived=False, unit=Unit("km"), equivalencies=[],
                  fvalidate='default', doc=None)
        """
        kw.setdefault("fvalidate", self._fvalidate_in)  # prefer the input fvalidate
        cloned = replace(self, **kw)
        # Transfer over the __set_name__ stuff. If `clone` is used to make a
        # new descriptor, __set_name__ will be called again, overwriting this.
        cloned.__set_name__(None, self.name)

        return cloned

    def __repr__(self) -> str:
        """Return repr(self)."""
        fields_repr = (
            # Get the repr, using the input fvalidate over the processed value
            f"{f.name}={(getattr(self, f.name if f.name != 'fvalidate' else '_fvalidate_in'))!r}"
            for f in fields(self)
            # Only show fields that should be displayed and are not sentinel values
            if f.repr and (f.name != "default" or self.default is not MISSING)
        )
        return f"{self.__class__.__name__}({', '.join(fields_repr)})"

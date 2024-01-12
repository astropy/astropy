# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any, TypeAlias, TypeVar, overload

import astropy.units as u

if TYPE_CHECKING:
    from astropy.cosmology import Cosmology, Parameter

__all__: list[str] = []

# Callable[[Cosmology, Parameter, Any], Any]
_FValidateCallable: TypeAlias = Callable[["Cosmology", "Parameter", Any], Any]
T = TypeVar("T")

_REGISTRY_FVALIDATORS: dict[str, _FValidateCallable] = {}


@overload
def _register_validator(
    key: str, fvalidate: _FValidateCallable
) -> _FValidateCallable: ...


@overload
def _register_validator(
    key: str, fvalidate: None = None
) -> Callable[[_FValidateCallable], _FValidateCallable]: ...


def _register_validator(
    key: str, fvalidate: _FValidateCallable | None = None
) -> _FValidateCallable | Callable[[_FValidateCallable], _FValidateCallable]:
    """Decorator to register a new kind of validator function.

    Parameters
    ----------
    key : str
    fvalidate : callable[[Cosmology, Parameter, Any], Any] or None, optional
        Value validation function.

    Returns
    -------
    ``validator`` or callable[``validator``]
        if validator is None returns a function that takes and registers a
        validator. This allows ``register_validator`` to be used as a
        decorator.
    """
    if key in _REGISTRY_FVALIDATORS:
        raise KeyError(f"validator {key!r} already registered with Parameter.")

    # fvalidate directly passed
    if fvalidate is not None:
        _REGISTRY_FVALIDATORS[key] = fvalidate
        return fvalidate

    # for use as a decorator
    def register(fvalidate: _FValidateCallable) -> _FValidateCallable:
        """Register validator function.

        Parameters
        ----------
        fvalidate : callable[[Cosmology, Parameter, Any], Any]
            Validation function.

        Returns
        -------
        ``validator``
        """
        _REGISTRY_FVALIDATORS[key] = fvalidate
        return fvalidate

    return register


# ======================================================================


@_register_validator("default")
def _validate_with_unit(cosmology: Cosmology, param: Parameter, value: Any) -> Any:
    """Default Parameter value validator.

    Adds/converts units if Parameter has a unit.
    """
    if param.unit is not None:
        with u.add_enabled_equivalencies(param.equivalencies):
            value = u.Quantity(value, param.unit)
    return value


@_register_validator("float")
def _validate_to_float(cosmology: Cosmology, param: Parameter, value: Any) -> Any:
    """Parameter value validator with units, and converted to float."""
    value = _validate_with_unit(cosmology, param, value)
    return float(value)


@_register_validator("scalar")
def _validate_to_scalar(cosmology: Cosmology, param: Parameter, value: Any) -> Any:
    """"""
    value = _validate_with_unit(cosmology, param, value)
    if not value.isscalar:
        raise ValueError(f"{param.name} is a non-scalar quantity")
    return value


@_register_validator("non-negative")
def _validate_non_negative(cosmology: Cosmology, param: Parameter, value: Any) -> Any:
    """Parameter value validator where value is a positive float."""
    value = _validate_to_float(cosmology, param, value)
    if value < 0.0:
        raise ValueError(f"{param.name} cannot be negative.")
    return value

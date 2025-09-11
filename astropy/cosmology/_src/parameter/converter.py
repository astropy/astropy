# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = [
    "validate_non_negative",
    "validate_to_float",
    "validate_to_scalar",
    "validate_with_unit",
]

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

from numpy.typing import NDArray

import astropy.units as u

if TYPE_CHECKING:
    import astropy.cosmology

FValidateCallable = Callable[[object, object, Any], Any]
_REGISTRY_FVALIDATORS: dict[str, FValidateCallable] = {}


def _register_validator(
    key: str, fvalidate: FValidateCallable | None = None
) -> FValidateCallable | Callable[[FValidateCallable], FValidateCallable]:
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
    if key in _REGISTRY_FVALIDATORS:
        raise KeyError(f"validator {key!r} already registered with Parameter.")

    # fvalidate directly passed
    if fvalidate is not None:
        _REGISTRY_FVALIDATORS[key] = fvalidate
        return fvalidate

    # for use as a decorator
    def register(fvalidate):
        """Register validator function.

        Parameters
        ----------
        fvalidate : callable[[object, object, Any], Any]
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
def validate_with_unit(
    cosmology: "astropy.cosmology.Cosmology",
    param: "astropy.cosmology.Parameter",
    value: Any,
) -> Any:
    """Default Parameter value validator.

    Adds/converts units if Parameter has a unit.
    """
    if param.unit is not None:
        with u.add_enabled_equivalencies(param.equivalencies):
            value = u.Quantity(value, param.unit)
    return value


@_register_validator("float")
def validate_to_float(
    cosmology: "astropy.cosmology.Cosmology",
    param: "astropy.cosmology.Parameter",
    value: Any,
) -> float:
    """Parameter value validator with units, and converted to float."""
    value = validate_with_unit(cosmology, param, value)
    return float(value)


@_register_validator("scalar")
def validate_to_scalar(
    cosmology: "astropy.cosmology.Cosmology",
    param: "astropy.cosmology.Parameter",
    value: Any,
) -> NDArray:
    """Parameter value validator where value is a scalar."""
    value = validate_with_unit(cosmology, param, value)
    if not value.isscalar:
        raise ValueError(f"{param.name} is a non-scalar quantity")
    return value


@_register_validator("non-negative")
def validate_non_negative(
    cosmology: "astropy.cosmology.Cosmology",
    param: "astropy.cosmology.Parameter",
    value: Any,
) -> float:
    """Parameter value validator where value is a positive float."""
    value = validate_to_float(cosmology, param, value)
    if value < 0.0:
        raise ValueError(f"{param.name} cannot be negative.")
    return value

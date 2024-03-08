# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import astropy.units as u

__all__ = []

FValidateCallable = Callable[[object, object, Any], Any]
_REGISTRY_FVALIDATORS: dict[str, FValidateCallable] = {}


def _register_validator(key, fvalidate=None):
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
def _validate_with_unit(cosmology, param, value):
    """Default Parameter value validator.

    Adds/converts units if Parameter has a unit.
    """
    if param.unit is not None:
        with u.add_enabled_equivalencies(param.equivalencies):
            value = u.Quantity(value, param.unit)
    return value


@_register_validator("float")
def _validate_to_float(cosmology, param, value):
    """Parameter value validator with units, and converted to float."""
    value = _validate_with_unit(cosmology, param, value)
    return float(value)


@_register_validator("scalar")
def _validate_to_scalar(cosmology, param, value):
    """"""
    value = _validate_with_unit(cosmology, param, value)
    if not value.isscalar:
        raise ValueError(f"{param.name} is a non-scalar quantity")
    return value


@_register_validator("non-negative")
def _validate_non_negative(cosmology, param, value):
    """Parameter value validator where value is a positive float."""
    value = _validate_to_float(cosmology, param, value)
    if value < 0.0:
        raise ValueError(f"{param.name} cannot be negative.")
    return value

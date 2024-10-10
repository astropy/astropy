# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Custom errors and exceptions for astropy.units."""

__all__ = [
    "UnitsError",
    "UnitConversionError",
    "UnitScaleError",
    "UnitTypeError",
    "UnitsWarning",
    "UnitParserWarning",
    "OGIPParserWarning",
    "OGIPInvalidMultiplicationWarning",
]

from astropy.utils.exceptions import AstropyWarning


class UnitsError(Exception):
    """
    The base class for unit-specific exceptions.
    """


class UnitConversionError(UnitsError, ValueError):
    """
    Used specifically for errors related to converting between units or
    interpreting units in terms of other units.
    """


class UnitScaleError(UnitsError, ValueError):
    """
    Used to catch the errors involving scaled units,
    which are not recognized by FITS format.
    """


class UnitTypeError(UnitsError, TypeError):
    """
    Used specifically for errors in setting to units not allowed by a class.

    E.g., would be raised if the unit of an `~astropy.coordinates.Angle`
    instances were set to a non-angular unit.
    """


class UnitsWarning(AstropyWarning):
    """
    The base class for unit-specific warnings.
    """


class UnitParserWarning(AstropyWarning):
    """Unit parser warnings"""


class OGIPParserWarning(UnitParserWarning):
    """Base class for OGIP parser warnings."""


class OGIPInvalidMultiplicationWarning(OGIPParserWarning):
    """Warning for OGIP multiplications missing a space between the factors."""

    def __init__(self, left: str, right: str) -> None:
        self.message = (
            f"if '{left}{right}' was meant to be a multiplication, "
            f"it should have been written as '{left} {right}'."
        )

    def __str__(self) -> str:
        return self.message

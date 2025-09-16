# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Custom errors and exceptions for astropy.units."""

__all__ = [
    "UnitConversionError",
    "UnitParserWarning",
    "UnitScaleError",
    "UnitTypeError",
    "UnitsError",
    "UnitsWarning",
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


class UnitParserWarning(UnitsWarning):
    """Unit parser warnings"""

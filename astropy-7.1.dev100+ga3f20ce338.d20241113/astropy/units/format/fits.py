# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "FITS" unit format.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from astropy.units.errors import UnitScaleError
from astropy.utils import classproperty

from . import core, generic, utils

if TYPE_CHECKING:
    from typing import Literal

    from astropy.units import UnitBase


class FITS(generic.Generic):
    """
    The FITS standard unit format.

    This supports the format defined in the Units section of the `FITS
    Standard <https://fits.gsfc.nasa.gov/fits_standard.html>`_.
    """

    @classproperty(lazy=True)
    def _units(cls) -> dict[str, UnitBase]:
        from astropy import units as u

        # add some units up-front for which we don't want to use prefixes
        # and that have different names from the astropy default.
        names = {"Celsius": u.deg_C, "deg C": u.deg_C}
        bases = [
            "m", "g", "s", "rad", "sr", "K", "A", "mol", "cd",
            "Hz", "J", "W", "V", "N", "Pa", "C", "Ohm", "S",
            "F", "Wb", "T", "H", "lm", "lx", "a", "yr", "eV",
            "pc", "Jy", "mag", "R", "bit", "byte", "G", "barn",
        ]  # fmt: skip
        prefixes = [
            "y", "z", "a", "f", "p", "n", "u", "m", "c", "d",
            "", "da", "h", "k", "M", "G", "T", "P", "E", "Z", "Y",
        ]  # fmt: skip

        special_cases = {"dbyte": u.Unit("dbyte", 0.1 * u.byte)}

        for key, _ in utils.get_non_keyword_units(bases, prefixes):
            names[key] = special_cases[key] if key in special_cases else getattr(u, key)
        simple_units = [
            "deg", "arcmin", "arcsec", "mas", "min", "h", "d", "Ry",
            "solMass", "u", "solLum", "solRad", "AU", "lyr", "count",
            "ct", "photon", "ph", "pixel", "pix", "D", "Sun", "chan",
            "bin", "voxel", "adu", "beam", "erg", "Angstrom", "angstrom",
        ]  # fmt: skip
        names.update((unit, getattr(u, unit)) for unit in simple_units)

        return names

    @classmethod
    def _parse_unit(cls, unit: str, detailed_exception: bool = True) -> UnitBase:
        cls._validate_unit(unit, detailed_exception=detailed_exception)
        return cls._units[unit]

    @classmethod
    def to_string(
        cls, unit: UnitBase, fraction: bool | Literal["inline"] = False
    ) -> str:
        # Remove units that aren't known to the format
        unit = cls._decompose_to_known_units(unit)

        parts = []

        base = np.log10(unit.scale)

        if base % 1.0 != 0.0:
            raise UnitScaleError(
                "The FITS unit format is not able to represent scales "
                "that are not powers of 10.  Multiply your data by "
                f"{unit.scale:e}."
            )
        elif unit.scale != 1.0:
            # We could override format_exponential_notation to set the
            # scale factor but that would give the wrong impression that
            # all values in FITS are set that way.  So, instead do it
            # here, and use a unity-scale unit for the rest.
            parts.append(f"10**{int(base)}")
            unit = core.CompositeUnit(1, unit.bases, unit.powers)

        if unit.bases:
            parts.append(super().to_string(unit, fraction=fraction))

        return cls._scale_unit_separator.join(parts)

    @classmethod
    def parse(cls, s: str, debug: bool = False) -> UnitBase:
        result = super().parse(s, debug)
        if hasattr(result, "function_unit"):
            raise ValueError("Function units are not yet supported for FITS units.")
        return result

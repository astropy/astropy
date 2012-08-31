# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles the "FITS" unit format.
"""
import warnings

from . import generic


class Fits(generic.Generic):
    """
    The FITS standard unit format.

    This supports the format defined in the Units section of the `FITS
    Standard <http://fits.gsfc.nasa.gov/fits_standard.html>`_.
    """
    name = 'fits'

    def __init__(self):
        if not '_units' in Fits.__dict__:
            Fits._units, Fits._deprecated_units = self._generate_unit_names()

        if not '_parser' in Fits.__dict__:
            Fits._parser = self._make_parser()

    @staticmethod
    def _generate_unit_names():
        from astropy.units import standard_units as u
        names = {}
        deprecated_names = set()

        bases = [
            'm', 'g', 's', 'rad', 'sr', 'K', 'A', 'mol', 'cd',
            'Hz', 'J', 'W', 'V', 'N', 'Pa', 'C', 'Ohm', 'S',
            'F', 'Wb', 'T', 'H', 'lm', 'lx', 'a', 'yr', 'eV',
            'pc', 'Jy', 'mag', 'R', 'bit', 'byte'
            ]
        deprecated_bases = [
            'G', 'barn'
            ]
        prefixes = [
            'y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', 'c', 'd',
            '', 'da', 'h', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']

        for base in bases + deprecated_bases:
            for prefix in prefixes:
                key = prefix + base
                names[key] = getattr(u, key)
        for base in deprecated_bases:
            for prefix in prefixes:
                deprecated_names.add(prefix + base)

        simple_units = [
            'deg', 'arcmin', 'arcsec', 'mas', 'min', 'h', 'd', 'Ry',
            'solMass', 'u', 'solLum', 'solRad', 'AU', 'lyr', 'count',
            'photon', 'ph', 'pixel', 'pix', 'D', 'Sun', 'chan', 'bin',
            'voxel', 'adu', 'beam'
            ]
        deprecated_units = [
            'erg', 'Angstrom'
        ]

        for unit in simple_units + deprecated_units:
            names[unit] = getattr(u, unit)
        for unit in deprecated_units:
            deprecated_names.add(unit)

        return names, deprecated_names

    @classmethod
    def _parse_unit(cls, s, loc, toks):
        unit = toks[0]
        if unit not in cls._units:
            raise ValueError(
                "Unit {0!r} not supported by the FITS standard.".format(unit))

        if unit in cls._deprecated_units:
            warnings.warn(
                "The unit {0!r} has been deprecated in the FITS "
                "standard.".format(unit),
                DeprecationWarning)

        return cls._units[unit]

    def _get_unit_name(self, unit):
        name = unit.get_format_name('fits')

        if name not in self._units:
            raise ValueError(
                "Unit {0!r} is not part of the FITS standard".format(name))

        if name in self._deprecated_units:
            warnings.warn(
                "The unit {0!r} has been deprecated in the FITS "
                "standard.".format(name),
                DeprecationWarning)

        return name

    def to_string(self, unit):
        from .. import core

        # The FITS standard only allows powers of 10 as a multiplier.

        if isinstance(unit, core.CompositeUnit):
            pairs = zip(unit.bases, unit.powers)
            pairs.sort(key=lambda x: x[1])

            s = self._format_unit_list(pairs)
        elif isinstance(unit, core.NamedUnit):
            s = self._get_unit_name(unit)

        return s

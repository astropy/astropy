# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles the "VOUnit" unit format.
"""
import warnings

from . import generic


class VOUnit(generic.Generic):
    """
    The proposed IVOA standard for units used by the VO.

    This is an implementation of `proposed IVOA standard for units
    <http://www.ivoa.net/Documents/VOUnits/>`_.
    """
    def __init__(self):
        if not '_units' in VOUnit.__dict__:
            VOUnit._units, VOUnit._deprecated_units = \
              VOUnit._generate_unit_names()

        if not '_parser' in VOUnit.__dict__:
            VOUnit._parser = self._make_parser()

    @staticmethod
    def _generate_unit_names():
        from astropy.units import standard_units as u
        names = {}
        deprecated_names = set()

        bases = [
            'm', 's', 'A', 'K', 'mol', 'cd', 'g', 'rad', 'sr',
            'Hz', 'N', 'Pa', 'J', 'W', 'C', 'V', 'S', 'F',
            'Wb', 'T', 'H', 'lm', 'lx', 'Ohm']
        prefixes = [
            'y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', 'c', 'd',
            '', 'da', 'h', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']

        for base in bases:
            for prefix in prefixes:
                key = prefix + base
                names[key] = getattr(u, key)

        simple_units = ['min', 'h', 'd',
            'a', 'yr', 'deg', 'arcsec', 'arcmin', 'deg', 'mas', 'AU',
            'pc', 'u', 'eV', 'Jy']
        deprecated_units = ['angstrom', 'Angstrom',
            'barn', 'erg', 'G', 'mag', 'solMass', 'solLum', 'solRad',
            'lyr', 'ct', 'count', 'photon', 'ph', 'R', 'pix',
            'pixel', 'D', 'Sun', 'chan', 'bin', 'voxel', 'bit',
            'byte', 'adu', 'beam']

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
                "Unit {0!r} not supported by the VOUnit standard.".format(unit))

        if unit in cls._deprecated_units:
            warnings.warn(
                "The use of unit {0!r} is discouraged by the "
                "VOUnit standard.".format(unit),
                DeprecationWarning)

        return cls._units[unit]

    def _get_unit_name(self, unit):
        name = unit.get_format_name('vounit')

        if name not in self._units:
            raise ValueError(
                "Unit {0!r} is not part of the VOUnit standard".format(name))

        if name in self._deprecated_units:
            warnings.warn(
                "The use of unit {0!r} is discouraged by the "
                "VOUnit standard.".format(name),
                DeprecationWarning)

        return name

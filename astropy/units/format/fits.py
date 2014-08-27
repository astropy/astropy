# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "FITS" unit format.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern.six.moves import zip

import keyword
import warnings

from ...utils.exceptions import AstropyDeprecationWarning

from . import generic
from . import utils
from ...utils.misc import did_you_mean

class UnitScaleError(ValueError):
    """
    Used to catch the errors involving scaled units,
    which are not recognized by FITS format.
    """

class Fits(generic.Generic):
    """
    The FITS standard unit format.

    This supports the format defined in the Units section of the `FITS
    Standard <http://fits.gsfc.nasa.gov/fits_standard.html>`_.
    """
    name = 'fits'

    def __init__(self):
        # Build this on the class, so it only gets generated once.
        if '_parser' not in Fits.__dict__:
            Fits._parser, Fits._lexer = self._make_parser()

        if not '_units' in Fits.__dict__:
            Fits._units, Fits._deprecated_units = self._generate_unit_names()

    @staticmethod
    def _generate_unit_names():
        from ... import units as u
        names = {}
        deprecated_names = set()

        bases = [
            'm', 'g', 's', 'rad', 'sr', 'K', 'A', 'mol', 'cd',
            'Hz', 'J', 'W', 'V', 'N', 'Pa', 'C', 'Ohm', 'S',
            'F', 'Wb', 'T', 'H', 'lm', 'lx', 'a', 'yr', 'eV',
            'pc', 'Jy', 'mag', 'R', 'bit', 'byte'
        ]
        deprecated_bases = ['G', 'barn']
        prefixes = [
            'y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', 'c', 'd',
            '', 'da', 'h', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']

        special_cases = {'dbyte': u.Unit('dbyte', 0.1*u.byte)}

        for base in bases + deprecated_bases:
            for prefix in prefixes:
                key = prefix + base
                if keyword.iskeyword(key):
                    continue
                elif key in special_cases:
                    names[key] = special_cases[key]
                else:
                    names[key] = getattr(u, key)
        for base in deprecated_bases:
            for prefix in prefixes:
                deprecated_names.add(prefix + base)

        simple_units = [
            'deg', 'arcmin', 'arcsec', 'mas', 'min', 'h', 'd', 'Ry',
            'solMass', 'u', 'solLum', 'solRad', 'AU', 'lyr', 'count',
            'ct', 'photon', 'ph', 'pixel', 'pix', 'D', 'Sun', 'chan',
            'bin', 'voxel', 'adu', 'beam'
        ]
        deprecated_units = ['erg', 'Angstrom', 'angstrom']

        for unit in simple_units + deprecated_units:
            names[unit] = getattr(u, unit)
        for unit in deprecated_units:
            deprecated_names.add(unit)

        return names, deprecated_names

    @classmethod
    def _parse_unit(cls, unit, detailed_exception=True):
        if unit not in cls._units:
            if detailed_exception:
                raise ValueError(
                    "Unit {0!r} not supported by the FITS standard. {1}".format(
                        unit, did_you_mean(
                            unit, cls._units)))
            else:
                raise ValueError()

        if unit in cls._deprecated_units:
            warnings.warn(
                "The unit {0!r} has been deprecated in the FITS "
                "standard.".format(unit),
                AstropyDeprecationWarning)

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
                AstropyDeprecationWarning)

        return name

    def to_string(self, unit):
        from .. import core

        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, self._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            if unit.scale != 1:
                raise UnitScaleError(
                    "The FITS unit format is not able to represent scale. "
                    "Multiply your data by {0:e}.".format(unit.scale))

            pairs = list(zip(unit.bases, unit.powers))
            pairs.sort(key=lambda x: x[1], reverse=True)

            s = self._format_unit_list(pairs)
        elif isinstance(unit, core.NamedUnit):
            s = self._get_unit_name(unit)

        return s

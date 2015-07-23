# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "FITS" unit format.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern.six.moves import zip

import numpy as np

import copy
import keyword

from . import generic
from . import utils


class Fits(generic.Generic):
    """
    The FITS standard unit format.

    This supports the format defined in the Units section of the `FITS
    Standard <http://fits.gsfc.nasa.gov/fits_standard.html>`_.
    """

    name = 'fits'

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

        return names, deprecated_names, []

    @classmethod
    def _validate_unit(cls, unit, detailed_exception=True):
        if unit not in cls._units:
            if detailed_exception:
                raise ValueError(
                    "Unit '{0}' not supported by the FITS standard. {1}".format(
                        unit, utils.did_you_mean_units(
                            unit, cls._units, cls._deprecated_units,
                            cls._to_decomposed_alternative)))
            else:
                raise ValueError()

        if unit in cls._deprecated_units:
            utils.unit_deprecation_warning(
                unit, cls._units[unit], 'FITS',
                cls._to_decomposed_alternative)

    @classmethod
    def _parse_unit(cls, unit, detailed_exception=True):
        cls._validate_unit(unit)
        return cls._units[unit]

    @classmethod
    def _get_unit_name(cls, unit):
        name = unit.get_format_name('fits')
        cls._validate_unit(name)
        return name

    @classmethod
    def to_string(cls, unit):
        from .. import core

        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, cls._get_unit_name)

        parts = []

        if isinstance(unit, core.CompositeUnit):
            base = np.log10(unit.scale)

            if base % 1.0 != 0.0:
                raise core.UnitScaleError(
                    "The FITS unit format is not able to represent scales "
                    "that are not powers of 10.  Multiply your data by "
                    "{0:e}.".format(unit.scale))
            elif unit.scale != 1.0:
                parts.append('10**{0}'.format(int(base)))

            pairs = list(zip(unit.bases, unit.powers))
            if len(pairs):
                pairs.sort(key=lambda x: x[1], reverse=True)
                parts.append(cls._format_unit_list(pairs))

            s = ' '.join(parts)
        elif isinstance(unit, core.NamedUnit):
            s = cls._get_unit_name(unit)

        return s

    @classmethod
    def _to_decomposed_alternative(cls, unit):
        from .. import core

        try:
            s = cls.to_string(unit)
        except core.UnitScaleError:
            scale = unit.scale
            unit = copy.copy(unit)
            unit._scale = 1.0
            return '{0} (with data multiplied by {1})'.format(
                cls.to_string(unit), scale)
        return s

    @classmethod
    def parse(cls, s, debug=False):
        result = super(Fits, cls).parse(s, debug)
        if hasattr(result, 'function_unit'):
            raise ValueError("Function units are not yet supported for "
                             "FITS units.")
        return result

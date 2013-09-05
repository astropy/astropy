# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles the "VOUnit" unit format.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import warnings

from . import generic
from . import utils


class VOUnit(generic.Generic):
    """
    The proposed IVOA standard for units used by the VO.

    This is an implementation of `proposed IVOA standard for units
    <http://www.ivoa.net/Documents/VOUnits/>`_.
    """
    def __init__(self):
        super(VOUnit, self).__init__()

        if not '_units' in VOUnit.__dict__:
            unit_names = VOUnit._generate_unit_names()
            VOUnit._units, VOUnit._deprecated_units = unit_names

    @staticmethod
    def _generate_unit_names():
        import keyword
        from ... import units as u
        names = {}
        deprecated_names = set()

        bases = [
            'm', 's', 'A', 'K', 'mol', 'cd', 'g', 'rad', 'sr', 'Hz', 'N', 'Pa',
            'J', 'W', 'C', 'V', 'S', 'F', 'Wb', 'T', 'H', 'lm', 'lx', 'Ohm']
        prefixes = [
            'y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', 'c', 'd',
            '', 'da', 'h', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']

        for base in bases:
            for prefix in prefixes:
                key = prefix + base
                if keyword.iskeyword(key):
                    continue
                names[key] = getattr(u, key)

        simple_units = [
            'min', 'h', 'd', 'a', 'yr', 'deg', 'arcsec', 'arcmin', 'deg',
            'mas', 'AU', 'pc', 'u', 'eV', 'Jy']

        deprecated_units = [
            'angstrom', 'Angstrom', 'barn', 'erg', 'G', 'mag', 'solMass',
            'solLum', 'solRad', 'lyr', 'ct', 'count', 'photon', 'ph', 'R',
            'pix', 'pixel', 'D', 'Sun', 'chan', 'bin', 'voxel', 'bit', 'byte',
            'adu', 'beam']

        for unit in simple_units + deprecated_units:
            names[unit] = getattr(u, unit)
        for unit in deprecated_units:
            deprecated_names.add(unit)

        return names, deprecated_names

    @classmethod
    def _parse_unit(cls, unit):
        if unit not in cls._units:
            raise ValueError(
                "Unit {0!r} not supported by the VOUnit "
                "standard.".format(unit))

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

    def to_string(self, unit):
        from .. import core

        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, self._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            s = ''
            if unit.scale != 1:
                m, ex = utils.split_mantissa_exponent(unit.scale)
                if m:
                    s += m + ' '
                if ex:
                    s += ' 10'
                    if not ex.startswith('-'):
                        s += '+'
                    s += ex

            pairs = zip(unit.bases, unit.powers)
            pairs.sort(key=lambda x: x[1], reverse=True)

            s += self._format_unit_list(pairs)
        elif isinstance(unit, core.NamedUnit):
            s = self._get_unit_name(unit)

        return s

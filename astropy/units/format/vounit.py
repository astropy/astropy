# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles the "VOUnit" unit format.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six
from ...extern.six.moves import zip

import keyword
import warnings
from ...utils.exceptions import AstropyDeprecationWarning

from . import generic
from . import utils
from ...utils.misc import did_you_mean


class VOUnit(generic.Generic):
    """
    The proposed IVOA standard for units used by the VO.

    This is an implementation of `proposed IVOA standard for units
    <http://www.ivoa.net/Documents/VOUnits/>`_.
    """
    def __init__(self):
        if '_parser' not in VOUnit.__dict__:
            VOUnit._parser, VOUnit._lexer = self._make_parser()

        if not '_units' in VOUnit.__dict__:
            unit_names = VOUnit._generate_unit_names()
            VOUnit._units, VOUnit._deprecated_units = unit_names

    @staticmethod
    def _generate_unit_names():
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
            'min', 'h', 'd', 'a', 'yr', 'arcsec', 'arcmin', 'deg',
            'mas', 'AU', 'pc', 'u', 'eV', 'Jy']

        deprecated_units = [
            'angstrom', 'Angstrom', 'barn', 'erg', 'G', 'mag', 'dB', 'solMass',
            'solLum', 'solRad', 'lyr', 'ct', 'count', 'photon', 'ph', 'R',
            'pix', 'pixel', 'D', 'Sun', 'chan', 'bin', 'voxel', 'bit', 'byte',
            'adu', 'beam']

        for unit in simple_units + deprecated_units:
            names[unit] = getattr(u, unit)
        for unit in deprecated_units:
            deprecated_names.add(unit)

        return names, deprecated_names

    def parse(self, s, debug=False):
        result = self._do_parse(s, debug=debug)
        if s.count('/') > 1:
            from ..core import UnitsError
            raise UnitsError(
                "'{0}' contains multiple slashes, which is "
                "disallowed by the VOUnit standard".format(s))
        return result

    @classmethod
    def _parse_unit(cls, unit, detailed_exception=True):
        if unit not in cls._units:
            if detailed_exception:
                raise ValueError(
                    "Unit {0!r} not supported by the VOUnit "
                    "standard. {1}".format(
                        unit, did_you_mean(
                            unit, cls._units)))
            else:
                raise ValueError()

        if unit in cls._deprecated_units:
            warnings.warn(
                "The use of unit {0!r} is discouraged by the "
                "VOUnit standard.".format(unit),
                AstropyDeprecationWarning)

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
                AstropyDeprecationWarning)

        return name

    def to_string(self, unit):
        from .. import core

        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, self._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            if unit.physical_type == 'dimensionless' and unit.scale != 1:
                raise ValueError("The VOUnit format is not able to "
                                 "represent scale for dimensionless units. "
                                 "Multiply your data by {0:e}."
                                 .format(unit.scale))
            s = ''
            if unit.scale != 1:
                m, ex = utils.split_mantissa_exponent(unit.scale)
                parts = []
                if m:
                    parts.append(m)
                if ex:
                    fex = '10'
                    if not ex.startswith('-'):
                        fex += '+'
                    fex += ex
                    parts.append(fex)
                s += ' '.join(parts)

            pairs = list(zip(unit.bases, unit.powers))
            pairs.sort(key=lambda x: x[1], reverse=True)

            s += self._format_unit_list(pairs)
        elif isinstance(unit, core.NamedUnit):
            s = self._get_unit_name(unit)

        return s

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles the "VOUnit" unit format.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six
from ...extern.six.moves import zip

import copy
import keyword
import operator
import re
import warnings

from . import core, generic, utils


class VOUnit(generic.Generic):
    """
    The IVOA standard for units used by the VO.

    This is an implementation of `Units in the VO 1.0
    <http://www.ivoa.net/Documents/VOUnits/>`_.
    """
    _explicit_custom_unit_regex = re.compile(
        r"^[YZEPTGMkhdcmunpfazy]?'((?!\d)\w)+'$")
    _custom_unit_regex = re.compile(r"^((?!\d)\w)+$")
    _custom_units = {}

    @staticmethod
    def _generate_unit_names():
        from ... import units as u
        from ...units import required_by_vounit as uvo

        names = {}
        deprecated_names = set()

        bases = [
            'A', 'C', 'D', 'F', 'G', 'H', 'Hz', 'J', 'Jy', 'K', 'N',
            'Ohm', 'Pa', 'R', 'Ry', 'S', 'T', 'V', 'W', 'Wb', 'a',
            'adu', 'arcmin', 'arcsec', 'barn', 'beam', 'bin', 'cd',
            'chan', 'count', 'ct', 'd', 'deg', 'eV', 'erg', 'g', 'h',
            'lm', 'lx', 'lyr', 'm', 'mag', 'min', 'mol', 'pc', 'ph',
            'photon', 'pix', 'pixel', 'rad', 'rad', 's', 'solLum',
            'solMass', 'solRad', 'sr', 'u', 'voxel', 'yr'
        ]
        binary_bases = [
            'bit', 'byte', 'B'
        ]
        simple_units = [
            'Angstrom', 'angstrom', 'AU', 'au', 'Ba', 'dB', 'mas'
        ]
        si_prefixes = [
            'y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', 'c', 'd',
            '', 'da', 'h', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'
        ]
        binary_prefixes = [
            'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei'
        ]
        deprecated_units = set([
            'a', 'angstrom', 'Angstrom', 'au', 'Ba', 'barn', 'ct',
            'erg', 'G', 'ph', 'pix'
        ])

        def do_defines(bases, prefixes, skips=[]):
            for base in bases:
                for prefix in prefixes:
                    key = prefix + base
                    if key in skips:
                        continue
                    if keyword.iskeyword(key):
                        continue

                    names[key] = getattr(u if hasattr(u, key) else uvo, key)
                    if base in deprecated_units:
                        deprecated_names.add(key)

        do_defines(bases, si_prefixes, ['pct', 'pcount', 'yd'])
        do_defines(binary_bases, si_prefixes + binary_prefixes, ['dB', 'dbyte'])
        do_defines(simple_units, [''])

        return names, deprecated_names, []

    @classmethod
    def parse(cls, s, debug=False):
        if s in ('unknown', 'UNKNOWN'):
            return None
        if s == '':
            return core.dimensionless_unscaled
        if s.count('/') > 1:
            raise core.UnitsError(
                "'{0}' contains multiple slashes, which is "
                "disallowed by the VOUnit standard".format(s))
        result = cls._do_parse(s, debug=debug)
        if hasattr(result, 'function_unit'):
            raise ValueError("Function units are not yet supported in "
                             "VOUnit.")
        return result

    @classmethod
    def _parse_unit(cls, unit, detailed_exception=True):
        if unit not in cls._units:
            if cls._explicit_custom_unit_regex.match(unit):
                return cls._def_custom_unit(unit)

            if not cls._custom_unit_regex.match(unit):
                raise ValueError()

            warnings.warn(
                "Unit {0!r} not supported by the VOUnit "
                "standard. {1}".format(
                    unit, utils.did_you_mean_units(
                        unit, cls._units, cls._deprecated_units,
                        cls._to_decomposed_alternative)),
                core.UnitsWarning)

            return cls._def_custom_unit(unit)

        if unit in cls._deprecated_units:
            utils.unit_deprecation_warning(
                unit, cls._units[unit], 'VOUnit',
                cls._to_decomposed_alternative)

        return cls._units[unit]

    @classmethod
    def _get_unit_name(cls, unit):
        # The da- and d- prefixes are discouraged.  This has the
        # effect of adding a scale to value in the result.
        if isinstance(unit, core.PrefixUnit):
            if unit._represents.scale == 10.0:
                raise ValueError(
                    "In '{0}': VOUnit can not represent units with the 'da' "
                    "(deka) prefix".format(unit))
            elif unit._represents.scale == 0.1:
                raise ValueError(
                    "In '{0}': VOUnit can not represent units with the 'd' "
                    "(deci) prefix".format(unit))

        name = unit.get_format_name('vounit')

        if unit in six.itervalues(cls._custom_units):
            return name

        if name not in cls._units:
            raise ValueError(
                "Unit {0!r} is not part of the VOUnit standard".format(name))

        if name in cls._deprecated_units:
            utils.unit_deprecation_warning(
                name, unit, 'VOUnit',
                cls._to_decomposed_alternative)

        return name

    @classmethod
    def _def_custom_unit(cls, unit):
        def def_base(name):
            if name in cls._custom_units:
                return cls._custom_units[name]

            if name.startswith("'"):
                return core.def_unit(
                    [name[1:-1], name],
                    format={'vounit': name},
                    namespace=cls._custom_units)
            else:
                return core.def_unit(
                    name, namespace=cls._custom_units)

        if unit in cls._custom_units:
            return cls._custom_units[unit]

        for short, full, factor in core.si_prefixes:
            for prefix in short:
                if unit.startswith(prefix):
                    base_name = unit[len(prefix):]
                    base_unit = def_base(base_name)
                    return core.PrefixUnit(
                        [prefix + x for x in base_unit.names],
                        core.CompositeUnit(factor, [base_unit], [1],
                                        _error_check=False),
                        format={'vounit': prefix + base_unit.names[-1]},
                        namespace=cls._custom_units)

        return def_base(unit)

    @classmethod
    def to_string(cls, unit):
        from .. import core

        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, cls._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            if unit.physical_type == 'dimensionless' and unit.scale != 1:
                raise core.UnitScaleError(
                    "The VOUnit format is not able to "
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
            pairs.sort(key=operator.itemgetter(1), reverse=True)

            s += cls._format_unit_list(pairs)
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

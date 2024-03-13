# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles the "VOUnit" unit format.
"""

import copy
import keyword
import re
import warnings

from . import core, generic, utils


class VOUnit(generic.Generic):
    """
    The IVOA standard for units used by the VO.

    This is an implementation of `Units in the VO 1.0
    <http://www.ivoa.net/documents/VOUnits/>`_.
    """

    _explicit_custom_unit_regex = re.compile(r"^[YZEPTGMkhdcmunpfazy]?'((?!\d)\w)+'$")
    _custom_unit_regex = re.compile(r"^((?!\d)\w)+$")
    _custom_units = {}
    _space = "."
    _scale_unit_separator = ""

    @staticmethod
    def _generate_unit_names():
        from astropy import units as u
        from astropy.units import required_by_vounit as uvo

        names = {}
        deprecated_names = set()
        # The tropical year is missing here compared to the standard
        bases = [
            "A", "a", "adu", "arcmin", "arcsec", "barn", "beam", "bin",
            "C", "cd", "chan", "count", "ct", "d", "D", "deg", "erg", "eV",
            "F", "g", "G", "H", "h", "Hz", "J", "Jy", "K", "lm", "lx", "lyr",
            "m", "mag", "min", "mol", "N", "Ohm", "Pa", "pc", "ph", "photon",
            "pix", "pixel", "R", "rad", "Ry", "s", "S", "solLum", "solMass",
            "solRad", "sr", "T", "u", "V", "voxel", "W", "Wb", "yr",
        ]  # fmt: skip
        binary_bases = ["bit", "byte", "B"]
        simple_units = ["Angstrom", "angstrom", "AU", "au", "Ba", "dB", "mas", "Sun"]
        si_prefixes = [
            "y", "z", "a", "f", "p", "n", "u", "m", "c", "d",
            "", "da", "h", "k", "M", "G", "T", "P", "E", "Z", "Y"
        ]  # fmt: skip
        # While zebi and yobi are part of the standard for binary prefixes,
        # they are not implemented here due to computation limitations
        binary_prefixes = ["Ki", "Mi", "Gi", "Ti", "Pi", "Ei"]
        deprecated_units = {
            "angstrom", "Angstrom", "Ba", "barn", "erg", "G", "ta",
        }  # fmt: skip

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

        do_defines(bases, si_prefixes, ["pct", "pcount", "yd"])
        do_defines(binary_bases, si_prefixes + binary_prefixes, ["dB", "dbyte"])
        do_defines(simple_units, [""])

        return names, deprecated_names, []

    @classmethod
    def parse(cls, s, debug=False):
        if s in ("unknown", "UNKNOWN"):
            return None
        if s == "":
            return core.dimensionless_unscaled
        # Check for excess solidi, but exclude fractional exponents (allowed)
        if s.count("/") > 1 and s.count("/") - len(re.findall(r"\(\d+/\d+\)", s)) > 1:
            raise core.UnitsError(
                f"'{s}' contains multiple slashes, which is "
                "disallowed by the VOUnit standard."
            )
        result = cls._do_parse(s, debug=debug)
        if hasattr(result, "function_unit"):
            raise ValueError("Function units are not yet supported in VOUnit.")
        return result

    @classmethod
    def _get_unit(cls, t):
        try:
            return super()._get_unit(t)
        except ValueError:
            if cls._explicit_custom_unit_regex.match(t.value):
                return cls._def_custom_unit(t.value)

            if cls._custom_unit_regex.match(t.value):
                warnings.warn(
                    f"Unit {t.value!r} not supported by the VOUnit standard. "
                    + utils.did_you_mean_units(
                        t.value,
                        cls._units,
                        cls._deprecated_units,
                        cls._to_decomposed_alternative,
                    ),
                    core.UnitsWarning,
                )

                return cls._def_custom_unit(t.value)

            raise

    @classmethod
    def _parse_unit(cls, unit, detailed_exception=True):
        if unit not in cls._units:
            raise ValueError()

        if unit in cls._deprecated_units:
            utils.unit_deprecation_warning(
                unit, cls._units[unit], "VOUnit", cls._to_decomposed_alternative
            )

        return cls._units[unit]

    @classmethod
    def _get_unit_name(cls, unit):
        # The da- and d- prefixes are discouraged.  This has the
        # effect of adding a scale to value in the result.
        if isinstance(unit, core.PrefixUnit):
            if unit._represents.scale == 10.0:
                raise ValueError(
                    f"In '{unit}': VOUnit can not represent units with the 'da' "
                    "(deka) prefix"
                )
            elif unit._represents.scale == 0.1:
                raise ValueError(
                    f"In '{unit}': VOUnit can not represent units with the 'd' "
                    "(deci) prefix"
                )

        name = super()._get_unit_name(unit)

        if unit in cls._custom_units.values():
            return name

        if name not in cls._units:
            raise ValueError(f"Unit {name!r} is not part of the VOUnit standard")

        if name in cls._deprecated_units:
            utils.unit_deprecation_warning(
                name, unit, "VOUnit", cls._to_decomposed_alternative
            )

        return name

    @classmethod
    def _def_custom_unit(cls, unit):
        def def_base(name):
            if name in cls._custom_units:
                return cls._custom_units[name]

            if name.startswith("'"):
                return core.def_unit(
                    [name[1:-1], name],
                    format={"vounit": name},
                    namespace=cls._custom_units,
                )
            else:
                return core.def_unit(name, namespace=cls._custom_units)

        if unit in cls._custom_units:
            return cls._custom_units[unit]

        for short, full, factor in core.si_prefixes:
            for prefix in short:
                if unit.startswith(prefix):
                    base_name = unit[len(prefix) :]
                    base_unit = def_base(base_name)
                    return core.PrefixUnit(
                        [prefix + x for x in base_unit.names],
                        core.CompositeUnit(
                            factor, [base_unit], [1], _error_check=False
                        ),
                        format={"vounit": prefix + base_unit.names[-1]},
                        namespace=cls._custom_units,
                    )

        return def_base(unit)

    @classmethod
    def _format_superscript(cls, number):
        return f"**({number})" if "/" in number or "." in number else f"**{number}"

    @classmethod
    def format_exponential_notation(cls, val, format_spec=".8g"):
        return super().format_exponential_notation(val, format_spec)

    @classmethod
    def _format_fraction(cls, scale, numerator, denominator, *, fraction="inline"):
        if not (fraction is True or fraction == "inline"):
            raise ValueError(
                "format {cls.name!r} only supports inline fractions,"
                f"not fraction={fraction!r}."
            )

        if cls._space in denominator:
            denominator = f"({denominator})"
        if scale and numerator == "1":
            return f"{scale}/{denominator}"
        return f"{scale}{numerator}/{denominator}"

    @classmethod
    def to_string(cls, unit, fraction=False):
        from astropy.units import core

        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, cls._get_unit_name)

        if unit.physical_type == "dimensionless" and unit.scale != 1:
            raise core.UnitScaleError(
                "The VOUnit format is not able to "
                "represent scale for dimensionless units. "
                f"Multiply your data by {unit.scale:e}."
            )

        return super().to_string(unit, fraction=fraction)

    @classmethod
    def _to_decomposed_alternative(cls, unit):
        from astropy.units import core

        try:
            s = cls.to_string(unit)
        except core.UnitScaleError:
            scale = unit.scale
            unit = copy.copy(unit)
            unit._scale = 1.0
            return f"{cls.to_string(unit)} (with data multiplied by {scale})"
        return s

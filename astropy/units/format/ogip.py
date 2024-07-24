# Licensed under a 3-clause BSD style license - see LICNSE.rst

# This module includes files automatically generated from ply (these end in
# _lextab.py and _parsetab.py). To generate these files, remove them from this
# folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/units
#
# You can then commit the changes to the re-generated _lextab.py and
# _parsetab.py files.

"""
Handles units in `Office of Guest Investigator Programs (OGIP)
FITS files
<https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
"""

from __future__ import annotations

import copy
import math
import warnings
from fractions import Fraction
from typing import TYPE_CHECKING

from astropy.utils import classproperty, parsing

from . import core, generic, utils

if TYPE_CHECKING:
    from typing import ClassVar


class OGIP(generic.Generic):
    """
    Support the units in `Office of Guest Investigator Programs (OGIP)
    FITS files
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
    """

    _tokens = (
        "DIVISION",
        "OPEN_PAREN",
        "CLOSE_PAREN",
        "WHITESPACE",
        "STARSTAR",
        "STAR",
        "SIGN",
        "UFLOAT",
        "LIT10",
        "UINT",
        "UNKNOWN",
        "UNIT",
    )

    _functions: ClassVar[frozenset[str]] = frozenset((
        "log", "ln", "exp", "sqrt", "sin", "cos", "tan",
        "asin", "acos", "atan", "sinh", "cosh", "tanh",
    ))  # fmt: skip

    @classmethod
    def _generate_unit_names(cls):
        from astropy import units as u

        bases = [
            "A", "C", "cd", "eV", "F", "g", "H", "Hz", "J",
            "Jy", "K", "lm", "lx", "m", "mol", "N", "ohm", "Pa",
            "pc", "rad", "s", "S", "sr", "T", "V", "W", "Wb",
        ]  # fmt: skip
        prefixes = [
            "y", "z", "a", "f", "p", "n", "u", "m", "c", "d",
            "", "da", "h", "k", "M", "G", "T", "P", "E", "Z", "Y",
        ]  # fmt: skip

        names = {
            unit: getattr(u, unit)
            for unit, _ in utils.get_non_keyword_units(bases, prefixes)
        }
        simple_units = [
            "angstrom", "arcmin", "arcsec", "AU", "barn", "bin",
            "byte", "chan", "count", "day", "deg", "erg", "G",
            "h", "lyr", "mag", "min", "photon", "pixel",
            "voxel", "yr",
        ]  # fmt: skip
        names.update((unit, getattr(u, unit)) for unit in simple_units)

        # Create a separate, disconnected unit for the special case of
        # Crab and mCrab, since OGIP doesn't define their quantities.
        names["Crab"] = u.def_unit(["Crab"], prefixes=False, doc="Crab (X-ray flux)")
        names["mCrab"] = u.Unit(10**-3 * names["Crab"])

        names.update((name, name) for name in cls._functions)

        return names, {"Crab", "mCrab"}

    @classproperty(lazy=True)
    def _lexer(cls):
        tokens = cls._tokens

        t_DIVISION = r"/"
        t_OPEN_PAREN = r"\("
        t_CLOSE_PAREN = r"\)"
        t_WHITESPACE = "[ \t]+"
        t_STARSTAR = r"\*\*"
        t_STAR = r"\*"

        # NOTE THE ORDERING OF THESE RULES IS IMPORTANT!!
        # Regular expression rules for simple tokens
        def t_UFLOAT(t):
            r"(((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+))|(((\d+\.\d*)|(\.\d+))([eE][+-]?\d+)?)"
            t.value = float(t.value)
            return t

        def t_UINT(t):
            r"\d+"
            t.value = int(t.value)
            return t

        def t_SIGN(t):
            r"[+-](?=\d)"
            t.value = float(t.value + "1")
            return t

        def t_X(t):  # multiplication for factor in front of unit
            r"[x×]"
            return t

        def t_LIT10(t):
            r"10"
            return 10

        def t_UNKNOWN(t):
            r"[Uu][Nn][Kk][Nn][Oo][Ww][Nn]"
            return None

        def t_UNIT(t):
            r"[a-zA-Z][a-zA-Z_]*"
            t.value = cls._get_unit(t)
            return t

        # Don't ignore whitespace
        t_ignore = ""

        # Error handling rule
        def t_error(t):
            raise ValueError(f"Invalid character at col {t.lexpos}")

        return parsing.lex(lextab="ogip_lextab", package="astropy/units")

    @classproperty(lazy=True)
    def _parser(cls):
        """
        The grammar here is based on the description in the
        `Specification of Physical Units within OGIP FITS files
        <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__,
        which is not terribly precise.  The exact grammar is here is
        based on the YACC grammar in the `unity library
        <https://bitbucket.org/nxg/unity/>`_.
        """
        tokens = cls._tokens

        def p_main(p):
            """
            main : UNKNOWN
                 | complete_expression
                 | scale_factor complete_expression
                 | scale_factor WHITESPACE complete_expression
            """
            if len(p) == 4:
                p[0] = p[1] * p[3]
            elif len(p) == 3:
                p[0] = p[1] * p[2]
            else:
                p[0] = p[1]

        def p_complete_expression(p):
            """
            complete_expression : product_of_units
            """
            p[0] = p[1]

        def p_product_of_units(p):
            """
            product_of_units : unit_expression
                             | division unit_expression
                             | product_of_units product unit_expression
                             | product_of_units division unit_expression
            """
            if len(p) == 4:
                if p[2] == "DIVISION":
                    p[0] = p[1] / p[3]
                else:
                    p[0] = p[1] * p[3]
            elif len(p) == 3:
                p[0] = p[2] ** -1
            else:
                p[0] = p[1]

        def p_unit_expression(p):
            """
            unit_expression : unit
                            | UNIT OPEN_PAREN complete_expression CLOSE_PAREN
                            | OPEN_PAREN complete_expression CLOSE_PAREN
                            | UNIT OPEN_PAREN complete_expression CLOSE_PAREN power numeric_power
                            | OPEN_PAREN complete_expression CLOSE_PAREN power numeric_power
            """
            bad_function_message = (
                "The function '{}' is valid in OGIP, but not understood "
                "by astropy.units."
            )
            bad_multiplication_message = (
                "if '{0}{1}' was meant to be a multiplication, "
                "it should have been written as '{0} {1}'."
            )

            if len(p) == 7:
                if p[1] == "sqrt":
                    p[0] = p[3] ** (0.5 * p[6])
                elif p[1] in cls._functions:
                    raise ValueError(bad_function_message.format(p[1]))
                else:
                    raise ValueError(
                        bad_multiplication_message.format(p[1], f"({p[3]})**{p[6]}")
                    )
            elif len(p) == 6:
                p[0] = p[2] ** p[5]
            elif len(p) == 5:
                if p[1] == "sqrt":
                    p[0] = p[3] ** 0.5
                elif p[1] in cls._functions:
                    raise ValueError(bad_function_message.format(p[1]))
                else:
                    raise ValueError(
                        bad_multiplication_message.format(p[1], f"({p[3]})")
                    )
            elif len(p) == 4:
                p[0] = p[2]
            else:
                p[0] = p[1]

        def p_scale_factor(p):
            """
            scale_factor : LIT10 power numeric_power
                         | LIT10
                         | signed_float
                         | signed_float power numeric_power
                         | signed_int power numeric_power
            """
            if len(p) == 4:
                p[0] = 10 ** p[3]
            else:
                p[0] = p[1]
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(p[0]) % 1.0 != 0.0:
                from astropy.units.core import UnitsWarning

                warnings.warn(
                    f"'{p[0]}' scale should be a power of 10 in OGIP format",
                    UnitsWarning,
                )

        def p_division(p):
            """
            division : DIVISION
                     | WHITESPACE DIVISION
                     | WHITESPACE DIVISION WHITESPACE
                     | DIVISION WHITESPACE
            """
            p[0] = "DIVISION"

        def p_product(p):
            """
            product : WHITESPACE
                    | STAR
                    | WHITESPACE STAR
                    | WHITESPACE STAR WHITESPACE
                    | STAR WHITESPACE
            """
            p[0] = "PRODUCT"

        def p_power(p):
            """
            power : STARSTAR
            """
            p[0] = "POWER"

        def p_unit(p):
            """
            unit : UNIT
                 | UNIT power numeric_power
            """
            if len(p) == 4:
                p[0] = p[1] ** p[3]
            else:
                p[0] = p[1]

        def p_numeric_power(p):
            """
            numeric_power : UINT
                          | signed_float
                          | OPEN_PAREN signed_int CLOSE_PAREN
                          | OPEN_PAREN signed_float CLOSE_PAREN
                          | OPEN_PAREN signed_float division UINT CLOSE_PAREN
            """
            if len(p) == 6:
                p[0] = Fraction(int(p[2]), int(p[4]))
            elif len(p) == 4:
                p[0] = p[2]
            else:
                p[0] = p[1]

        def p_sign(p):
            """
            sign : SIGN
                 |
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1.0

        def p_signed_int(p):
            """
            signed_int : SIGN UINT
            """
            p[0] = p[1] * p[2]

        def p_signed_float(p):
            """
            signed_float : sign UINT
                         | sign UFLOAT
            """
            p[0] = p[1] * p[2]

        def p_error(p):
            raise ValueError()

        return parsing.yacc(tabmodule="ogip_parsetab", package="astropy/units")

    @classmethod
    def _parse_unit(cls, unit, detailed_exception=True):
        cls._validate_unit(unit, detailed_exception=detailed_exception)
        return cls._units[unit]

    @classmethod
    def parse(cls, s, debug=False):
        s = s.strip()
        try:
            # This is a short circuit for the case where the string is
            # just a single unit name
            return cls._parse_unit(s, detailed_exception=False)
        except ValueError:
            try:
                return core.Unit(cls._parser.parse(s, lexer=cls._lexer, debug=debug))
            except ValueError as e:
                if str(e):
                    raise
                else:
                    raise ValueError(f"Syntax error parsing unit '{s}'")

    @classmethod
    def _format_superscript(cls, number):
        return f"**({number})" if "/" in number else f"**{number}"

    @classmethod
    def to_string(cls, unit, fraction="inline"):
        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, cls._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(unit.scale) % 1.0 != 0.0:
                warnings.warn(
                    f"'{unit.scale}' scale should be a power of 10 in OGIP format",
                    core.UnitsWarning,
                )

        return super().to_string(unit, fraction=fraction)

    @classmethod
    def _to_decomposed_alternative(cls, unit):
        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, cls._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(unit.scale) % 1.0 != 0.0:
                scale = unit.scale
                unit = copy.copy(unit)
                unit._scale = 1.0
                return (
                    f"{generic._to_string(cls, unit)} (with data multiplied by {scale})"
                )

        return super().to_string(unit)

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

import math
import warnings
from fractions import Fraction
from typing import TYPE_CHECKING

from astropy.units.core import CompositeUnit
from astropy.units.errors import UnitParserWarning, UnitsWarning
from astropy.utils import classproperty, parsing

from . import utils
from .base import Base, _ParsingFormatMixin

if TYPE_CHECKING:
    from typing import ClassVar, Literal

    import numpy as np

    from astropy.extern.ply.lex import Lexer
    from astropy.units import UnitBase
    from astropy.units.typing import UnitScale
    from astropy.utils.parsing import ThreadSafeParser


class OGIP(Base, _ParsingFormatMixin):
    """
    Support the units in `Office of Guest Investigator Programs (OGIP)
    FITS files
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
    """

    _tokens: ClassVar[tuple[str, ...]] = (
        "DIVISION",
        "OPEN_PAREN",
        "CLOSE_PAREN",
        "WHITESPACE",
        "POWER",
        "STAR",
        "SIGN",
        "UFLOAT",
        "LIT10",
        "UINT",
        "UNKNOWN",
        "FUNCNAME",
        "UNIT",
    )

    _deprecated_units: ClassVar[frozenset[str]] = frozenset(("Crab", "mCrab"))

    @classproperty(lazy=True)
    def _units(cls) -> dict[str, UnitBase]:
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
            "byte", "chan", "count", "d", "deg", "erg", "G",
            "h", "lyr", "mag", "min", "photon", "pixel",
            "voxel", "yr",
        ]  # fmt: skip
        names.update((unit, getattr(u, unit)) for unit in simple_units)

        # Create a separate, disconnected unit for the special case of
        # Crab and mCrab, since OGIP doesn't define their quantities.
        names["Crab"] = u.def_unit(["Crab"], prefixes=False, doc="Crab (X-ray flux)")
        names["mCrab"] = u.Unit(10**-3 * names["Crab"])

        return names

    @classproperty(lazy=True)
    def _lexer(cls) -> Lexer:
        tokens = cls._tokens

        t_DIVISION = "[ \t]*/[ \t]*"
        t_OPEN_PAREN = r"\("
        t_CLOSE_PAREN = r"\)"
        t_WHITESPACE = "[ \t]+"
        t_POWER = r"\*\*"
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
            t.value = 1 if t.value == "+" else -1
            return t

        def t_LIT10(t):
            r"10"
            return 10

        def t_UNKNOWN(t):
            r"[Uu][Nn][Kk][Nn][Oo][Ww][Nn]"
            return None

        def t_FUNCNAME(t):
            r"((sqrt)|(ln)|(exp)|(log)|(sin)|(cos)|(tan)|(asin)|(acos)|(atan)|(sinh)|(cosh)|(tanh))(?=\ *\()"
            return t

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
    def _parser(cls) -> ThreadSafeParser:
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
            match p[1:]:
                case (factor, unit) | (factor, _, unit):
                    p[0] = CompositeUnit(factor * unit.scale, unit.bases, unit.powers)
                case _:
                    p[0] = p[1]

        def p_complete_expression(p):
            """
            complete_expression : unit_expression
                                | product_of_units
                                | division_of_units
            """
            # product_of_units is not in unit_expression for performance
            # division_of_units is separate to enforce the correct order of operations
            p[0] = p[1]

        def p_product_of_units(p):
            """
            product_of_units : complete_expression product unit_expression
            """
            p[0] = p[1] * p[3]

        def p_division_of_units(p):
            """
            division_of_units : DIVISION unit_expression
                              | complete_expression DIVISION unit_expression
            """
            match p[1:]:
                case _, unit:
                    p[0] = unit**-1
                case num, _, denom:
                    p[0] = num / denom

        def p_unit_expression(p):
            """
            unit_expression : UNIT
                            | function
                            | UNIT POWER numeric_power
                            | UNIT OPEN_PAREN complete_expression CLOSE_PAREN
                            | OPEN_PAREN complete_expression CLOSE_PAREN
                            | UNIT OPEN_PAREN complete_expression CLOSE_PAREN POWER numeric_power
                            | OPEN_PAREN complete_expression CLOSE_PAREN POWER numeric_power
            """
            bad_multiplication_message = (
                "if '{0}{1}' was meant to be a multiplication, "
                "it should have been written as '{0} {1}'."
            )

            match p[1:]:
                case factor, _, unit, _, _, power:
                    warnings.warn(
                        bad_multiplication_message.format(factor, f"({unit})**{power}"),
                        UnitParserWarning,
                    )
                    p[0] = factor * unit**power
                case (_, unit, _, _, power) | (unit, "**", power):
                    p[0] = unit**power
                case left, _, right, _:
                    warnings.warn(
                        bad_multiplication_message.format(left, f"({right})"),
                        UnitParserWarning,
                    )
                    p[0] = left * right
                case _, unit, _:
                    p[0] = unit
                case _:
                    p[0] = p[1]

        def p_function(p):
            """
            function : FUNCNAME OPEN_PAREN complete_expression CLOSE_PAREN
                     | FUNCNAME OPEN_PAREN complete_expression CLOSE_PAREN POWER numeric_power
            """
            match p[1:]:
                case "sqrt", _, unit, _:
                    p[0] = unit**0.5
                case "sqrt", _, unit, _, _, numeric_power:
                    p[0] = unit ** (0.5 * numeric_power)
                case func, *_:
                    raise ValueError(
                        f"The function '{func}' is valid in OGIP, but not understood "
                        "by astropy.units."
                    )

        def p_scale_factor(p):
            """
            scale_factor : LIT10 POWER numeric_power
                         | LIT10
                         | signed_float
                         | signed_float POWER numeric_power
                         | signed_int POWER numeric_power
            """
            if len(p) == 4:
                p[0] = 10 ** p[3]
            else:
                p[0] = p[1]
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(p[0]) % 1.0 != 0.0:
                warnings.warn(
                    f"'{p[0]}' scale should be a power of 10 in OGIP format",
                    UnitsWarning,
                )

        def p_product(p):
            """
            product : WHITESPACE
                    | STAR
                    | WHITESPACE STAR
                    | WHITESPACE STAR WHITESPACE
                    | STAR WHITESPACE
            """

        def p_numeric_power(p):
            """
            numeric_power : UINT
                          | signed_float
                          | OPEN_PAREN signed_int CLOSE_PAREN
                          | OPEN_PAREN signed_float CLOSE_PAREN
                          | OPEN_PAREN signed_float DIVISION UINT CLOSE_PAREN
            """
            if len(p) == 6:
                p[0] = Fraction(int(p[2]), int(p[4]))
            elif len(p) == 4:
                p[0] = p[2]
            else:
                p[0] = p[1]
                if p[1] < 0:
                    warnings.warn(
                        UnitParserWarning(
                            "negative exponents must be enclosed in parenthesis. "
                            f"Expected '**({p[1]})' instead of '**{p[1]}'."
                        )
                    )

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
    def parse(cls, s: str, debug: bool = False) -> UnitBase:
        return cls._do_parse(s.strip(), debug)

    @classmethod
    def _format_superscript(cls, number: str) -> str:
        return f"**({number})" if "/" in number else f"**{number}"

    @classmethod
    def to_string(
        cls, unit: UnitBase, fraction: bool | Literal["inline", "multiline"] = "inline"
    ) -> str:
        # Remove units that aren't known to the format
        unit = cls._decompose_to_known_units(unit)

        if isinstance(unit, CompositeUnit):
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(unit.scale) % 1.0 != 0.0:
                warnings.warn(
                    f"'{unit.scale}' scale should be a power of 10 in OGIP format",
                    UnitsWarning,
                )

        return super().to_string(unit, fraction=fraction)

    @classmethod
    def format_exponential_notation(
        cls, val: UnitScale | np.number, format_spec: str = "g"
    ) -> str:
        return format(val, format_spec)

    @classmethod
    def _validate_unit(cls, unit: str, detailed_exception: bool = True) -> UnitBase:
        if unit in cls._deprecated_units:
            warnings.warn(
                f"The unit '{unit}' has been deprecated in the OGIP standard.",
                UnitsWarning,
            )
        return super()._validate_unit(unit, detailed_exception)

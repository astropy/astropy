# Licensed under a 3-clause BSD style license - see LICENSE.rst

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
Handles a "generic" string format for units
"""

from __future__ import annotations

import re
import unicodedata
import warnings
from fractions import Fraction
from typing import TYPE_CHECKING

from astropy.units.errors import UnitsWarning
from astropy.utils import classproperty, parsing
from astropy.utils.misc import did_you_mean

from . import core
from .base import Base

if TYPE_CHECKING:
    from re import Match, Pattern
    from typing import ClassVar, Final

    import numpy as np

    from astropy.extern.ply.lex import Lexer, LexToken
    from astropy.units import CompositeUnit, NamedUnit, UnitBase
    from astropy.utils.parsing import ThreadSafeParser


class Generic(Base):
    """
    A "generic" format.

    The syntax of the format is based directly on the FITS standard,
    but instead of only supporting the units that FITS knows about, it
    supports any unit available in the `astropy.units` namespace.
    """

    _tokens: ClassVar[tuple[str, ...]] = (
        "COMMA",
        "DOUBLE_STAR",
        "STAR",
        "PERIOD",
        "SOLIDUS",
        "CARET",
        "OPEN_PAREN",
        "CLOSE_PAREN",
        "FUNCNAME",
        "UNIT",
        "SIGN",
        "UINT",
        "UFLOAT",
    )

    _deprecated_units: ClassVar[frozenset[str]] = frozenset()

    @classproperty(lazy=True)
    def _lexer(cls) -> Lexer:
        tokens = cls._tokens

        t_COMMA = r"\,"
        t_STAR = r"\*"
        t_PERIOD = r"\."
        t_SOLIDUS = r"/"
        t_DOUBLE_STAR = r"\*\*"
        t_CARET = r"\^"
        t_OPEN_PAREN = r"\("
        t_CLOSE_PAREN = r"\)"

        # NOTE THE ORDERING OF THESE RULES IS IMPORTANT!!
        # Regular expression rules for simple tokens
        def t_UFLOAT(t):
            r"((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+)?"
            if not re.search(r"[eE\.]", t.value):
                t.type = "UINT"
                t.value = int(t.value)
            elif t.value.endswith("."):
                t.type = "UINT"
                t.value = int(t.value[:-1])
            else:
                t.value = float(t.value)
            return t

        def t_UINT(t):
            r"\d+"
            t.value = int(t.value)
            return t

        def t_SIGN(t):
            r"[+-](?=\d)"
            t.value = int(t.value + "1")
            return t

        # This needs to be a function so we can force it to happen
        # before t_UNIT
        def t_FUNCNAME(t):
            r"((sqrt)|(ln)|(exp)|(log)|(mag)|(dB)|(dex))(?=\ *\()"
            return t

        # A possible unit is something that consists of characters not used
        # for anything else: no spaces, no digits, signs, periods, stars,
        # carets, parentheses or commas.
        def t_UNIT(t):
            r"[^\s\d+\-\./\*\^\(\)\,]+"
            t.value = cls._get_unit(t)
            return t

        t_ignore = " "

        # Error handling rule
        def t_error(t):
            raise ValueError(f"Invalid character at col {t.lexpos}")

        return parsing.lex(
            lextab="generic_lextab", package="astropy/units", reflags=int(re.UNICODE)
        )

    @classproperty(lazy=True)
    def _parser(cls) -> ThreadSafeParser:
        """
        The grammar here is based on the description in the `FITS
        standard
        <http://fits.gsfc.nasa.gov/standard30/fits_standard30aa.pdf>`_,
        Section 4.3, which is not terribly precise.  The exact grammar
        is here is based on the YACC grammar in the `unity library
        <https://bitbucket.org/nxg/unity/>`_.

        This same grammar is used by the `"fits"` and `"vounit"`
        formats, the only difference being the set of available unit
        strings.
        """
        tokens = cls._tokens

        def p_main(p):
            """
            main : unit
                 | structured_unit
                 | structured_subunit
            """
            if isinstance(p[1], tuple):
                # Unpack possible StructuredUnit inside a tuple, ie.,
                # ignore any set of very outer parentheses.
                p[0] = p[1][0]
            else:
                p[0] = p[1]

        def p_structured_subunit(p):
            """
            structured_subunit : OPEN_PAREN structured_unit CLOSE_PAREN
            """
            # We hide a structured unit enclosed by parentheses inside
            # a tuple, so that we can easily distinguish units like
            # "(au, au/day), yr" from "au, au/day, yr".
            p[0] = (p[2],)

        def p_structured_unit(p):
            """
            structured_unit : subunit COMMA
                            | subunit COMMA subunit
            """
            from astropy.units.structured import StructuredUnit

            inputs = (p[1],) if len(p) == 3 else (p[1], p[3])
            units = ()
            for subunit in inputs:
                if isinstance(subunit, tuple):
                    # Structured unit that should be its own entry in the
                    # new StructuredUnit (was enclosed in parentheses).
                    units += subunit
                elif isinstance(subunit, StructuredUnit):
                    # Structured unit whose entries should be
                    # individually added to the new StructuredUnit.
                    units += subunit.values()
                else:
                    # Regular unit to be added to the StructuredUnit.
                    units += (subunit,)

            p[0] = StructuredUnit(units)

        def p_subunit(p):
            """
            subunit : unit
                    | structured_unit
                    | structured_subunit
            """
            p[0] = p[1]

        def p_unit(p):
            """
            unit : product_of_units
                 | factor product_of_units
                 | factor product product_of_units
                 | division_product_of_units
                 | factor division_product_of_units
                 | factor product division_product_of_units
                 | inverse_unit
                 | factor inverse_unit
                 | factor product inverse_unit
                 | factor
            """
            from astropy.units.core import CompositeUnit, Unit

            if len(p) == 2:
                p[0] = Unit(p[1])
            elif len(p) == 3:
                p[0] = CompositeUnit(p[1] * p[2].scale, p[2].bases, p[2].powers)
            elif len(p) == 4:
                p[0] = CompositeUnit(p[1] * p[3].scale, p[3].bases, p[3].powers)

        def p_division_product_of_units(p):
            """
            division_product_of_units : division_product_of_units division product_of_units
                                      | product_of_units
            """
            from astropy.units.core import Unit

            if len(p) == 4:
                p[0] = Unit(p[1] / p[3])
            else:
                p[0] = p[1]

        def p_inverse_unit(p):
            """
            inverse_unit : division unit_expression
            """
            p[0] = p[2] ** -1

        def p_factor(p):
            """
            factor : factor_fits
                   | factor_float
                   | factor_int
            """
            p[0] = p[1]

        def p_factor_float(p):
            """
            factor_float : signed_float
                         | signed_float UINT signed_int
                         | signed_float UINT power numeric_power
            """
            if cls.name == "fits":
                raise ValueError("Numeric factor not supported by FITS")
            if len(p) == 4:
                p[0] = p[1] * p[2] ** float(p[3])
            elif len(p) == 5:
                p[0] = p[1] * p[2] ** float(p[4])
            elif len(p) == 2:
                p[0] = p[1]

        def p_factor_int(p):
            """
            factor_int : UINT
                       | UINT signed_int
                       | UINT power numeric_power
                       | UINT UINT signed_int
                       | UINT UINT power numeric_power
            """
            if cls.name == "fits":
                raise ValueError("Numeric factor not supported by FITS")
            if len(p) == 2:
                p[0] = p[1]
            elif len(p) == 3:
                p[0] = p[1] ** float(p[2])
            elif len(p) == 4:
                if isinstance(p[2], int):
                    p[0] = p[1] * p[2] ** float(p[3])
                else:
                    p[0] = p[1] ** float(p[3])
            elif len(p) == 5:
                p[0] = p[1] * p[2] ** p[4]

        def p_factor_fits(p):
            """
            factor_fits : UINT power OPEN_PAREN signed_int CLOSE_PAREN
                        | UINT power OPEN_PAREN UINT CLOSE_PAREN
                        | UINT power signed_int
                        | UINT power UINT
                        | UINT SIGN UINT
                        | UINT OPEN_PAREN signed_int CLOSE_PAREN
            """
            if p[1] != 10:
                if cls.name == "fits":
                    raise ValueError("Base must be 10")
                else:
                    return
            if len(p) == 4:
                if p[2] in ("**", "^"):
                    p[0] = 10 ** p[3]
                else:
                    p[0] = 10 ** (p[2] * p[3])
            elif len(p) == 5:
                p[0] = 10 ** p[3]
            elif len(p) == 6:
                p[0] = 10 ** p[4]

        def p_product_of_units(p):
            """
            product_of_units : unit_expression product product_of_units
                             | unit_expression product_of_units
                             | unit_expression
            """
            if len(p) == 2:
                p[0] = p[1]
            elif len(p) == 3:
                p[0] = p[1] * p[2]
            else:
                p[0] = p[1] * p[3]

        def p_unit_expression(p):
            """
            unit_expression : function
                            | unit_with_power
                            | OPEN_PAREN product_of_units CLOSE_PAREN
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = p[2]

        def p_unit_with_power(p):
            """
            unit_with_power : UNIT power numeric_power
                            | UNIT numeric_power
                            | UNIT
            """
            if len(p) == 2:
                p[0] = p[1]
            elif len(p) == 3:
                p[0] = p[1] ** p[2]
            else:
                p[0] = p[1] ** p[3]

        def p_numeric_power(p):
            """
            numeric_power : sign UINT
                          | OPEN_PAREN paren_expr CLOSE_PAREN
            """
            if len(p) == 3:
                p[0] = p[1] * p[2]
            elif len(p) == 4:
                p[0] = p[2]

        def p_paren_expr(p):
            """
            paren_expr : sign UINT
                       | signed_float
                       | frac
            """
            if len(p) == 3:
                p[0] = p[1] * p[2]
            else:
                p[0] = p[1]

        def p_frac(p):
            """
            frac : sign UINT division sign UINT
            """
            p[0] = Fraction(p[1] * p[2], p[4] * p[5])

        def p_sign(p):
            """
            sign : SIGN
                 |
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1

        def p_product(p):
            """
            product : STAR
                    | PERIOD
            """

        def p_division(p):
            """
            division : SOLIDUS
            """

        def p_power(p):
            """
            power : DOUBLE_STAR
                  | CARET
            """
            p[0] = p[1]

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

        def p_function_name(p):
            """
            function_name : FUNCNAME
            """
            p[0] = p[1]

        def p_function(p):
            """
            function : function_name OPEN_PAREN main CLOSE_PAREN
            """
            if p[1] == "sqrt":
                p[0] = p[3] ** 0.5
                return
            elif p[1] in ("mag", "dB", "dex"):
                function_unit = cls._parse_unit(p[1])
                # In Generic, this is callable, but that does not have to
                # be the case in subclasses (e.g., in VOUnit it is not).
                if callable(function_unit):
                    p[0] = function_unit(p[3])
                    return

            raise ValueError(f"'{p[1]}' is not a recognized function")

        def p_error(p):
            raise ValueError()

        return parsing.yacc(tabmodule="generic_parsetab", package="astropy/units")

    @classmethod
    def _get_unit(cls, t: LexToken) -> UnitBase:
        try:
            return cls._parse_unit(t.value)
        except ValueError as e:
            registry = core.get_current_unit_registry()
            if t.value in registry.aliases:
                return registry.aliases[t.value]

            raise ValueError(f"At col {t.lexpos}, {str(e)}")

    @classmethod
    def _parse_unit(cls, s: str, detailed_exception: bool = True) -> UnitBase:
        registry = core.get_current_unit_registry().registry
        if s in cls._unit_symbols:
            s = cls._unit_symbols[s]

        elif not s.isascii():
            if s[0] == "\N{MICRO SIGN}":
                s = "u" + s[1:]
            elif s[0] == "°":
                s = "deg" if len(s) == 1 else "deg_" + s[1:]
            if s[-1] in cls._prefixable_unit_symbols:
                s = s[:-1] + cls._prefixable_unit_symbols[s[-1]]
            elif len(s) > 1 and s[-1] in cls._unit_suffix_symbols:
                s = s[:-1] + cls._unit_suffix_symbols[s[-1]]
            elif s.endswith("R\N{INFINITY}"):
                s = s[:-2] + "Ry"

        if s in registry:
            return registry[s]

        if detailed_exception:
            raise ValueError(f"{s} is not a valid unit. {did_you_mean(s, registry)}")
        else:
            raise ValueError()

    _unit_symbols: ClassVar[dict[str, str]] = {
        "%": "percent",
        "\N{PRIME}": "arcmin",
        "\N{DOUBLE PRIME}": "arcsec",
        "\N{MODIFIER LETTER SMALL H}": "hourangle",
        "e\N{SUPERSCRIPT MINUS}": "electron",
    }

    _prefixable_unit_symbols: ClassVar[dict[str, str]] = {
        "\N{GREEK CAPITAL LETTER OMEGA}": "Ohm",
        "\N{LATIN CAPITAL LETTER A WITH RING ABOVE}": "Angstrom",
        "\N{SCRIPT SMALL L}": "l",
    }

    _unit_suffix_symbols: ClassVar[dict[str, str]] = {
        "\N{CIRCLED DOT OPERATOR}": "sun",
        "\N{SUN}": "sun",
        "\N{CIRCLED PLUS}": "earth",
        "\N{EARTH}": "earth",
        "\N{JUPITER}": "jupiter",
        "\N{LATIN SUBSCRIPT SMALL LETTER E}": "_e",
        "\N{LATIN SUBSCRIPT SMALL LETTER P}": "_p",
    }

    _translations: ClassVar[dict[int, str]] = str.maketrans(
        {
            "\N{GREEK SMALL LETTER MU}": "\N{MICRO SIGN}",
            "\N{MINUS SIGN}": "-",
        }
    )
    """Character translations that should be applied before parsing a string.

    Note that this does explicitly *not* generally translate MICRO SIGN to u,
    since then a string like 'µ' would be interpreted as unit mass.
    """

    _superscripts: Final[str] = (
        "\N{SUPERSCRIPT MINUS}"
        "\N{SUPERSCRIPT PLUS SIGN}"
        "\N{SUPERSCRIPT ZERO}"
        "\N{SUPERSCRIPT ONE}"
        "\N{SUPERSCRIPT TWO}"
        "\N{SUPERSCRIPT THREE}"
        "\N{SUPERSCRIPT FOUR}"
        "\N{SUPERSCRIPT FIVE}"
        "\N{SUPERSCRIPT SIX}"
        "\N{SUPERSCRIPT SEVEN}"
        "\N{SUPERSCRIPT EIGHT}"
        "\N{SUPERSCRIPT NINE}"
    )

    _superscript_translations: ClassVar[dict[int, int]] = str.maketrans(
        _superscripts, "-+0123456789"
    )
    _regex_superscript: ClassVar[Pattern[str]] = re.compile(
        f"[{_superscripts}]?[{_superscripts[2:]}]+"
    )

    @classmethod
    def _convert_superscript(cls, m: Match[str]) -> str:
        return f"({m.group().translate(cls._superscript_translations)})"

    @classmethod
    def parse(cls, s: str, debug: bool = False) -> UnitBase:
        if not isinstance(s, str):
            s = s.decode("ascii")
        elif not s.isascii():
            # common normalization of unicode strings to avoid
            # having to deal with multiple representations of
            # the same character. This normalizes to "composed" form
            # and will e.g. convert OHM SIGN to GREEK CAPITAL LETTER OMEGA
            s = unicodedata.normalize("NFC", s)
            # Translate some basic unicode items that we'd like to support on
            # input but are not standard.
            s = s.translate(cls._translations)

            # TODO: might the below be better done in the parser/lexer?
            # Translate superscripts to parenthesized numbers; this ensures
            # that mixes of superscripts and regular numbers fail.
            s = cls._regex_superscript.sub(cls._convert_superscript, s)

        result = cls._do_parse(s, debug=debug)
        # Check for excess solidi, but exclude fractional exponents (accepted)
        n_slashes = s.count("/")
        if n_slashes > 1 and (n_slashes - len(re.findall(r"\(\d+/\d+\)", s))) > 1:
            warnings.warn(
                f"'{s}' contains multiple slashes, which is "
                "discouraged by the FITS standard",
                UnitsWarning,
            )
        return result

    @classmethod
    def _do_parse(cls, s: str, debug: bool = False) -> UnitBase:
        try:
            return cls._parser.parse(s, lexer=cls._lexer, debug=debug)
        except ValueError as e:
            if str(e):
                raise
            else:
                raise ValueError(f"Syntax error parsing unit '{s}'")

    @classmethod
    def _get_unit_name(cls, unit: NamedUnit) -> str:
        name = unit._get_format_name(cls.name)
        cls._validate_unit(name)
        return name

    @classmethod
    def _validate_unit(cls, unit: str, detailed_exception: bool = True) -> None:
        if unit not in cls._units:
            if detailed_exception:
                raise ValueError(
                    f"Unit '{unit}' not supported by the {cls.__name__} standard. "
                    + cls._did_you_mean_units(unit)
                )
            raise ValueError()
        if unit in cls._deprecated_units:
            message = (
                f"The unit '{unit}' has been deprecated in the {cls.__name__} standard."
            )
            if (decomposed := cls._try_decomposed(cls._units[unit])) is not None:
                message += f" Suggested: {decomposed}."
            warnings.warn(message, UnitsWarning)

    @classmethod
    def _did_you_mean_units(cls, unit: str) -> str:
        """
        A wrapper around `astropy.utils.misc.did_you_mean` that deals with
        the display of deprecated units.

        Parameters
        ----------
        unit : str
            The invalid unit string

        Returns
        -------
        msg : str
            A message with alternatives, or the empty string.
        """
        return did_you_mean(unit, cls._units, fix=cls._fix_deprecated)

    @classmethod
    def _fix_deprecated(cls, x: str) -> list[str]:
        return [x + " (deprecated)" if x in cls._deprecated_units else x]

    @classmethod
    def _try_decomposed(cls, unit: UnitBase) -> str | None:
        return None

    @classmethod
    def _decompose_to_known_units(cls, unit: CompositeUnit | NamedUnit) -> UnitBase:
        """
        Partially decomposes a unit so it is only composed of units that
        are "known" to a given format.
        """
        if isinstance(unit, core.CompositeUnit):
            return core.CompositeUnit(
                unit.scale,
                [cls._decompose_to_known_units(base) for base in unit.bases],
                unit.powers,
                _error_check=False,
            )
        if isinstance(unit, core.NamedUnit):
            try:
                cls._get_unit_name(unit)
            except ValueError:
                if isinstance(unit, core.Unit):
                    return cls._decompose_to_known_units(unit._represents)
                raise
            return unit
        raise TypeError(
            f"unit argument must be a 'NamedUnit' or 'CompositeUnit', not {type(unit)}"
        )

    @classmethod
    def format_exponential_notation(
        cls, val: float | np.number, format_spec: str = "g"
    ) -> str:
        return format(val, format_spec)

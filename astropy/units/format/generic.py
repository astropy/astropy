# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles a "generic" string format for units
"""

import functools
import operator
import re
import unicodedata
import warnings
from fractions import Fraction
from re import Match, Pattern
from typing import ClassVar, Final

import numpy as np
from lark import Token, Transformer, v_args

from astropy.units.core import CompositeUnit, Unit, UnitBase, get_current_unit_registry
from astropy.units.enums import DeprecatedUnitAction
from astropy.units.errors import UnitsWarning
from astropy.units.typing import UnitScale
from astropy.utils import classproperty
from astropy.utils.parsing import token_values

from .base import Base, _ParsingFormatMixin

_GRAMMAR: Final[str] = r"""
main: unit
    | structured_unit
    | structured_subunit

structured_subunit: "(" structured_unit ")"

structured_unit: subunit ","
               | subunit "," subunit

subunit: unit
       | structured_unit
       | structured_subunit

unit: factor
    | factor _PRODUCT? units
    | units

?units: product_of_units
      | division_product_of_units
      | inverse_unit

division_product_of_units: product_of_units ("/" product_of_units)+

inverse_unit: "/" unit_expression

?factor: factor_fits
       | factor_float
       | factor_int

factor_float: signed_float
            | signed_float UINT signed_int
            | signed_float UINT _POWER numeric_power

factor_int: UINT
          | UINT _POWER numeric_power         -> factor_int_power
          | UINT UINT signed_int              -> factor_int_scaled_sign
          | UINT UINT _POWER numeric_power    -> factor_int_scaled_power

factor_fits: UINT _POWER "(" signed_int ")"
           | UINT _POWER "(" UINT ")"
           | UINT _POWER signed_int
           | UINT _POWER UINT
           | UINT "(" signed_int ")"
           | UINT SIGN UINT                   -> factor_fits_sign

product_of_units: unit_expression (_PRODUCT? unit_expression)*

?unit_expression: function
                | unit_with_power
                | "(" product_of_units ")"

unit_with_power: UNIT
               | UNIT _POWER? numeric_power

numeric_power: sign UINT
             | "(" paren_expr ")"

?paren_expr: signed_float
           | frac

frac: sign UINT "/" sign UINT

sign: SIGN?

signed_int: SIGN UINT

signed_float: sign UINT
            | sign UFLOAT

function: FUNCNAME "(" main ")"

// The relative priorities of the terminals matter for the ones that can
// match the same text, since the first matching terminal is used.
_POWER.2: "**" | "^"
_PRODUCT: "*" | "."
FUNCNAME.2: /(sqrt|ln|exp|log|mag|dB|dex)(?= *\()/
UNIT: /[^\s\d+\-.\/*^(),]+/
SIGN: /[+-](?=\d)/
UFLOAT.2: /(\d+\.?\d*|\.\d+)[eE][+-]?\d+|\d+\.\d+|\.\d+/
UINT: /\d+\.?/
%ignore " "
"""


@v_args(wrapper=token_values)
class _GenericTransformer(Transformer):
    """Turn a parsed generic unit string into a unit.

    Terminal methods are applied as soon as a token is created and set its
    value, rule methods are applied when the rule is reduced.
    """

    def __init__(self, format_cls: type[_ParsingFormatMixin]) -> None:
        super().__init__()
        self._format = format_cls

    def t_UINT(self, t: Token) -> Token:
        t.value = int(t.rstrip("."))
        return t

    def t_UFLOAT(self, t: Token) -> Token:
        t.value = float(t)
        return t

    def t_SIGN(self, t: Token) -> Token:
        t.value = int(t + "1")
        return t

    def t_UNIT(self, t: Token) -> Token:
        t.value = self._format._get_unit(t)
        return t

    def main(self, unit):
        # Unpack possible StructuredUnit inside a tuple, ie., ignore any set
        # of very outer parentheses.
        return unit[0] if isinstance(unit, tuple) else unit

    def structured_subunit(self, unit):
        # We hide a structured unit enclosed by parentheses inside a tuple,
        # so that we can easily distinguish units like "(au, au/day), yr"
        # from "au, au/day, yr".
        return (unit,)

    def structured_unit(self, *inputs):
        from astropy.units.structured import StructuredUnit

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

        return StructuredUnit(units)

    def subunit(self, unit):
        return unit

    def unit(self, *items):
        match items:
            case [unit]:
                return Unit(unit)
            case [factor, unit]:
                return CompositeUnit(factor * unit.scale, unit.bases, unit.powers)

    def division_product_of_units(self, *units):
        return functools.reduce(operator.truediv, units)

    def inverse_unit(self, unit):
        return unit**-1

    def factor_float(self, *items):
        if self._format.name == "fits":
            raise ValueError("Numeric factor not supported by FITS")
        match items:
            case [factor]:
                return factor
            case [factor, base, power]:
                return factor * base ** float(power)

    def _check_factor_int(self) -> None:
        if self._format.name == "fits":
            raise ValueError("Numeric factor not supported by FITS")

    def factor_int(self, factor):
        self._check_factor_int()
        return factor

    def factor_int_power(self, base, power):
        self._check_factor_int()
        return base ** float(power)

    def factor_int_scaled_sign(self, factor, base, power):
        self._check_factor_int()
        return factor * base ** float(power)

    def factor_int_scaled_power(self, factor, base, power):
        self._check_factor_int()
        return factor * base**power

    def factor_fits(self, base, power):
        if base != 10:
            if self._format.name == "fits":
                raise ValueError("Base must be 10")
            else:
                return None
        return 10**power

    def factor_fits_sign(self, base, sign, power):
        return self.factor_fits(base, sign * power)

    def product_of_units(self, *units):
        # Multiply from the right, to get the same order as in the past.
        return functools.reduce(lambda product, unit: unit * product, reversed(units))

    def unit_with_power(self, unit, power=None):
        return unit if power is None else unit**power

    def numeric_power(self, *items):
        match items:
            case [sign, uint]:
                return sign * uint
            case [power]:
                return power

    def frac(self, sign1, uint1, sign2, uint2):
        return Fraction(sign1 * uint1, sign2 * uint2)

    def sign(self, sign=1):
        return sign

    def signed_int(self, sign, uint):
        return sign * uint

    def signed_float(self, sign, number):
        return sign * number

    def function(self, name, unit):
        if name == "sqrt":
            return unit**0.5
        elif name in ("mag", "dB", "dex"):
            try:
                function_unit = self._format._validate_unit(name)
            except KeyError:
                raise ValueError(
                    self._format._invalid_unit_error_message(name)
                ) from None
            # In Generic, this is callable, but that does not have to
            # be the case in subclasses (e.g., in VOUnit it is not).
            if callable(function_unit):
                return function_unit(unit)

        raise ValueError(f"'{name}' is not a recognized function")


class _GenericParserMixin(_ParsingFormatMixin):
    """Provide the parser used by Generic, FITS and VOUnit.

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

    _grammar: ClassVar[str] = _GRAMMAR
    _transformer: ClassVar[type[Transformer]] = _GenericTransformer


class Generic(Base, _GenericParserMixin):
    """
    A "generic" format.

    The syntax of the format is based directly on the FITS standard,
    but instead of only supporting the units that FITS knows about, it
    supports any unit available in the `astropy.units` namespace.
    """

    @classproperty
    def _units(cls) -> dict[str, UnitBase]:
        return get_current_unit_registry().registry

    @classmethod
    def _validate_unit(
        cls, s: str, deprecations: DeprecatedUnitAction = DeprecatedUnitAction.WARN
    ) -> UnitBase:
        if s in cls._unit_symbols:
            s = cls._unit_symbols[s]

        elif not s.isascii():
            if s[0].startswith("°"):
                s = "deg" if len(s) == 1 else "deg_" + s[1:]
            if len(s) > 1 and s[-1] in cls._unit_suffix_symbols:
                s = s[:-1] + cls._unit_suffix_symbols[s[-1]]
            elif s.endswith("R\N{INFINITY}"):
                s = s[:-2] + "Ry"

        return cls._units[s]

    @classmethod
    def _invalid_unit_error_message(cls, unit: str) -> str:
        return f"{unit} is not a valid unit. {cls._did_you_mean_units(unit)}"

    _unit_symbols: ClassVar[dict[str, str]] = {
        "%": "percent",
        "\N{PRIME}": "arcmin",
        "\N{DOUBLE PRIME}": "arcsec",
        "\N{MODIFIER LETTER SMALL H}": "hourangle",
        "e\N{SUPERSCRIPT MINUS}": "electron",
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

    _translations: ClassVar[dict[int, str]] = str.maketrans({"\N{MINUS SIGN}": "-"})
    """Character translations that should be applied before parsing a string."""

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

        result = cls._do_parse(s, debug)

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
    def format_exponential_notation(
        cls, val: UnitScale | np.number, format_spec: str = "g"
    ) -> str:
        return format(val, format_spec)

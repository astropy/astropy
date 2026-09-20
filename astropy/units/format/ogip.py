# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles units in `Office of Guest Investigator Programs (OGIP)
FITS files
<https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
"""

import math
import warnings
from fractions import Fraction
from typing import ClassVar, Final, Literal

import numpy as np
from lark import Token, Transformer, v_args

from astropy.units.core import CompositeUnit, UnitBase
from astropy.units.enums import DeprecatedUnitAction
from astropy.units.errors import UnitParserWarning, UnitsWarning
from astropy.units.typing import UnitScale
from astropy.utils import classproperty
from astropy.utils.parsing import token_values

from .base import Base, _ParsingFormatMixin

_GRAMMAR: Final[str] = r"""
main: complete_expression
    | scale_factor complete_expression                  -> scaled_unit
    | scale_factor _WHITESPACE complete_expression      -> scaled_unit

// product_of_units is not in unit_expression for performance
// division_of_units is separate to enforce the correct order of operations
?complete_expression: unit_expression
                    | product_of_units
                    | division_of_units

product_of_units: complete_expression _product unit_expression

division_of_units: _DIVISION unit_expression                   -> inverse_unit
                 | complete_expression _DIVISION unit_expression

unit_expression: UNIT
               | function
               | "(" complete_expression ")"
               | UNIT _POWER numeric_power                              -> unit_power
               | "(" complete_expression ")" _POWER numeric_power       -> unit_power
               | UNIT "(" complete_expression ")"                       -> implicit_product
               | UNIT "(" complete_expression ")" _POWER numeric_power  -> implicit_product_power

function: FUNCNAME "(" complete_expression ")"
        | FUNCNAME "(" complete_expression ")" _POWER numeric_power

scale_factor: number
            | number _POWER numeric_power   -> scale_factor_power

_product: _WHITESPACE
        | _STAR
        | _WHITESPACE _STAR
        | _WHITESPACE _STAR _WHITESPACE
        | _STAR _WHITESPACE

numeric_power: number
             | "(" number ")"              -> numeric_power_in_parentheses
             | "(" INT _DIVISION INT ")"   -> numeric_power_fraction

?number: INT
       | FLOAT

// The relative priorities of the terminals matter for the ones that can
// match the same text, since the first matching terminal is used.
FLOAT.9: /[+-]?((((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+))|(((\d+\.\d*)|(\.\d+))([eE][+-]?\d+)?))/
INT.8: /[+-]?\d+/
FUNCNAME.7: /(sqrt|ln|exp|log|sin|cos|tan|asin|acos|atan|sinh|cosh|tanh)(?= *\()/
UNIT.6: /[a-zA-Z][a-zA-Z_]*/
_DIVISION.5: /[ \t]*\/[ \t]*/
_WHITESPACE.4: /[ \t]+/
_POWER.3: "**"
_STAR.2: "*"
"""


@v_args(wrapper=token_values)
class _OGIPTransformer(Transformer):
    """Turn a parsed OGIP unit string into a unit.

    Terminal methods are applied as soon as a token is created and set its
    value, rule methods are applied when the rule is reduced.
    """

    _bad_multiplication_message: Final[str] = (
        "if '{0}{1}' was meant to be a multiplication, "
        "it should have been written as '{0} {1}'."
    )

    def __init__(self, format_cls: type[_ParsingFormatMixin]) -> None:
        super().__init__()
        self._format = format_cls

    def t_INT(self, t: Token) -> Token:
        t.value = int(t)
        return t

    def t_FLOAT(self, t: Token) -> Token:
        t.value = float(t)
        return t

    def t_UNIT(self, t: Token) -> Token:
        t.value = self._format._get_unit(t)
        return t

    def main(self, unit):
        return unit

    def scaled_unit(self, factor, unit):
        return CompositeUnit(factor * unit.scale, unit.bases, unit.powers)

    def product_of_units(self, unit1, unit2):
        return unit1 * unit2

    def inverse_unit(self, unit):
        return unit**-1

    def division_of_units(self, numerator, denominator):
        return numerator / denominator

    def unit_expression(self, unit):
        return unit

    def unit_power(self, unit, power):
        return unit**power

    def implicit_product(self, left, right):
        warnings.warn(
            self._bad_multiplication_message.format(left, f"({right})"),
            UnitParserWarning,
        )
        return left * right

    def implicit_product_power(self, factor, unit, power):
        warnings.warn(
            self._bad_multiplication_message.format(factor, f"({unit})**{power}"),
            UnitParserWarning,
        )
        return factor * unit**power

    def function(self, name, unit, power=1):
        if name == "sqrt":
            return unit ** (0.5 * power)
        raise ValueError(
            f"The function '{name}' is valid in OGIP, but not understood "
            "by astropy.units."
        )

    def scale_factor(self, factor):
        # Can't use np.log10 here, because the factor may be a Python long.
        if math.log10(factor) % 1.0 != 0.0:
            warnings.warn(
                f"'{factor}' scale should be a power of 10 in OGIP format",
                UnitsWarning,
            )
        return factor

    def scale_factor_power(self, base, power):
        return self.scale_factor(10**power)

    def numeric_power(self, power):
        if power < 0:
            warnings.warn(
                UnitParserWarning(
                    "negative exponents must be enclosed in parenthesis. "
                    f"Expected '**({power})' instead of '**{power}'."
                )
            )
        return power

    def numeric_power_in_parentheses(self, power):
        return power

    def numeric_power_fraction(self, numerator, denominator):
        return Fraction(numerator, denominator)


class OGIP(Base, _ParsingFormatMixin):
    """
    Support the units in `Office of Guest Investigator Programs (OGIP)
    FITS files
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
    """

    _deprecated_units: ClassVar[frozenset[str]] = frozenset(("Crab", "mCrab"))

    @classproperty(lazy=True)
    def _units(cls) -> dict[str, UnitBase]:
        from astropy import units as u

        names = {"as": u.attosecond}
        for non_prefixed_unit in [
            "angstrom", "arcmin", "arcsec", "AU", "barn", "bin",
            "byte", "chan", "count", "d", "deg", "erg", "G",
            "h", "lyr", "mag", "min", "photon", "pixel",
            "voxel", "yr",
        ]:  # fmt: skip
            names[non_prefixed_unit] = getattr(u, non_prefixed_unit)

        bases = [
            "A", "C", "cd", "eV", "F", "g", "H", "Hz", "J",
            "Jy", "K", "lm", "lx", "m", "mol", "N", "ohm", "Pa",
            "pc", "rad", "s", "S", "sr", "T", "V", "W", "Wb",
        ]  # fmt: skip
        prefixes = [
            "y", "z", "a", "f", "p", "n", "u", "m", "c", "d",
            "", "da", "h", "k", "M", "G", "T", "P", "E", "Z", "Y",
        ]  # fmt: skip

        for name in (prefix + base for base in bases for prefix in prefixes):
            if name not in names:
                names[name] = getattr(u, name)

        # Create a separate, disconnected unit for the special case of
        # Crab and mCrab, since OGIP doesn't define their quantities.
        names["Crab"] = u.def_unit(["Crab"], prefixes=False, doc="Crab (X-ray flux)")
        names["mCrab"] = u.Unit(10**-3 * names["Crab"])

        return names

    _grammar: ClassVar[str] = _GRAMMAR
    _transformer: ClassVar[type[Transformer]] = _OGIPTransformer

    @classmethod
    def parse(cls, s: str, debug: bool = False) -> UnitBase:
        return cls._do_parse(s.strip(), debug)

    @classmethod
    def _format_superscript(cls, number: str) -> str:
        return f"**({number})" if "/" in number else f"**{number}"

    @classmethod
    def to_string(
        cls,
        unit: UnitBase,
        fraction: bool | Literal["inline", "multiline"] = "inline",
        deprecations: DeprecatedUnitAction = DeprecatedUnitAction.WARN,
    ) -> str:
        # Remove units that aren't known to the format
        unit = cls._decompose_to_known_units(unit, deprecations)

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

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Handles the CDS string format for units."""

from typing import ClassVar, Final, Literal

from lark import Token, Transformer, v_args

from astropy.units.core import CompositeUnit, Unit, UnitBase
from astropy.units.enums import DeprecatedUnitAction
from astropy.units.utils import is_effectively_unity
from astropy.utils import classproperty
from astropy.utils.parsing import token_values

from .base import Base, _ParsingFormatMixin

_GRAMMAR: Final[str] = r"""
main: factor combined_units             -> scaled_unit
    | combined_units
    | DIMENSIONLESS
    | "[" combined_units "]"            -> dex_unit
    | "[" DIMENSIONLESS "]"             -> dex_unit
    | factor

?combined_units: product_of_units
               | division_of_units

product_of_units: unit_expression _PRODUCT combined_units
                | unit_expression

division_of_units: _DIVISION unit_expression                  -> inverse_unit
                 | combined_units _DIVISION unit_expression

?unit_expression: unit_with_power
                | "(" combined_units ")"

factor: FLOAT _X INT INT    -> scaled_factor
      | INT _X INT INT      -> scaled_factor
      | INT INT             -> power_of_ten
      | INT
      | FLOAT

unit_with_power: UNIT INT
               | UNIT

// The relative priorities of the terminals matter for the ones that can
// match the same text, since the first matching terminal is used.
_PRODUCT: "."
_DIVISION: "/"
FLOAT.5: /[+-]?((\d+\.\d+|\.\d+)([eE][+-]?\d+)?|\d{2,}[eE][+-]?\d+)/
INT.4: /[+-]?\d+/
_X.3: /[x×]/
// Most units are just combinations of letters with no numbers, but there
// are a few special ones (\h is Planck constant) and three that end in 0.
UNIT.2: /%|°|\\h|(a|eps|mu)0|((?!\d)\w)+/
// These are separate from UNIT since they cannot have a prefactor.
DIMENSIONLESS.1: /---|-/
"""


@v_args(wrapper=token_values)
class _CDSTransformer(Transformer):
    """Turn a parsed CDS unit string into a unit.

    Terminal methods are applied as soon as a token is created and set its
    value, rule methods are applied when the rule is reduced.
    """

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

    t_DIMENSIONLESS = t_UNIT

    def main(self, unit):
        return Unit(unit)

    def scaled_unit(self, factor, unit):
        return CompositeUnit(factor * unit.scale, unit.bases, unit.powers)

    def dex_unit(self, unit):
        from astropy.units import dex

        return dex(unit)

    def product_of_units(self, unit, *other):
        return unit * other[0] if other else unit

    def inverse_unit(self, unit):
        return unit**-1

    def division_of_units(self, numerator, denominator):
        return numerator / denominator

    def scaled_factor(self, factor, base, exponent):
        return factor * self.power_of_ten(base, exponent)

    def power_of_ten(self, base, exponent):
        if base != 10:
            raise ValueError("Only base ten exponents are allowed in CDS")
        return 10.0**exponent

    def factor(self, factor):
        return factor

    def unit_with_power(self, unit, power=None):
        return unit if power is None else unit**power


class CDS(Base, _ParsingFormatMixin):
    """
    Support the `Centre de Données astronomiques de Strasbourg
    <https://cds.unistra.fr/>`_ `Standards for Astronomical
    Catalogues 2.0 <https://vizier.unistra.fr/vizier/doc/catstd-3.2.htx>`_
    format, and the `complete set of supported units
    <https://vizier.unistra.fr/viz-bin/Unit>`_.  This format is used
    by VOTable up to version 1.2.
    """

    _space: ClassVar[str] = "."
    _times: ClassVar[str] = "x"
    _scale_unit_separator: ClassVar[str] = ""

    @classproperty(lazy=True)
    def _units(cls) -> dict[str, UnitBase]:
        from astropy import units as u
        from astropy.units import cds

        return {k: v for k, v in cds.__dict__.items() if isinstance(v, u.UnitBase)}

    _grammar: ClassVar[str] = _GRAMMAR
    _transformer: ClassVar[type[Transformer]] = _CDSTransformer

    @classmethod
    def parse(cls, s: str, debug: bool = False) -> UnitBase:
        if " " in s:
            raise ValueError("CDS unit must not contain whitespace")
        if not isinstance(s, str):
            s = s.decode("ascii")

        return cls._do_parse(s, debug)

    @classmethod
    def _format_mantissa(cls, m: str) -> str:
        return "" if m == "1" else m

    @classmethod
    def _format_superscript(cls, number: str) -> str:
        return number if number.startswith("-") else "+" + number

    @classmethod
    def to_string(
        cls,
        unit: UnitBase,
        fraction: bool | Literal["inline", "multiline"] = False,
        deprecations: DeprecatedUnitAction = DeprecatedUnitAction.WARN,
    ) -> str:
        # Remove units that aren't known to the format
        unit = cls._decompose_to_known_units(unit)

        if not unit.bases:
            if unit.scale == 1:
                return "---"
            elif is_effectively_unity(unit.scale * 100.0):
                return "%"

        return super().to_string(unit, fraction=fraction)

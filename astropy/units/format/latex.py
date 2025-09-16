# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "LaTeX" unit format.
"""

import re
from typing import ClassVar, Literal

from astropy.units.core import NamedUnit, UnitBase
from astropy.units.enums import DeprecatedUnitAction
from astropy.units.typing import UnitPower

from . import console


class Latex(console.Console):
    """
    Output LaTeX to display the unit based on IAU style guidelines.

    Attempts to follow the `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_.
    """

    _space: ClassVar[str] = r"\,"
    _scale_unit_separator: ClassVar[str] = r"\,"
    _times: ClassVar[str] = r" \times "

    @classmethod
    def _format_mantissa(cls, m: str) -> str:
        return m.replace("nan", r"{\rm NaN}").replace("inf", r"\infty")

    @classmethod
    def _format_superscript(cls, number: str) -> str:
        return f"^{{{number}}}"

    @classmethod
    def _format_unit_power(cls, unit: NamedUnit, power: UnitPower = 1) -> str:
        name = unit._get_format_name("latex")
        if name == unit.name:
            # This doesn't escape arbitrary LaTeX strings, but it should
            # be good enough for unit names which are required to be alpha
            # + "_" anyway.
            name = name.replace("_", r"\_")
        if power != 1:
            # If the LaTeX representation of the base unit already ends with
            # a superscript, we need to spell out the unit to avoid double
            # superscripts. For example, the logic below ensures that
            # `u.deg**2` returns `deg^{2}` instead of `{}^{\circ}^{2}`.
            if re.match(r".*\^{[^}]*}$", name):  # ends w/ superscript?
                name = unit.short_names[0]
            name += cls._format_power(power)
        return name

    @classmethod
    def _format_multiline_fraction(
        cls, scale: str, numerator: str, denominator: str
    ) -> str:
        return rf"{scale}\frac{{{numerator}}}{{{denominator}}}"

    @classmethod
    def to_string(
        cls,
        unit: UnitBase,
        fraction: bool | Literal["inline", "multiline"] = "multiline",
        deprecations: DeprecatedUnitAction = DeprecatedUnitAction.WARN,
    ) -> str:
        s = super().to_string(unit, fraction=fraction)
        return rf"$\mathrm{{{s}}}$"


class LatexInline(Latex):
    """
    Output LaTeX to display the unit based on IAU style guidelines with negative
    powers.

    Attempts to follow the `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_ and the
    `ApJ and AJ style guide
    <https://journals.aas.org/manuscript-preparation/>`_.
    """

    name: ClassVar[str] = "latex_inline"

    @classmethod
    def to_string(
        cls,
        unit: UnitBase,
        fraction: bool | Literal["inline", "multiline"] = False,
        deprecations: DeprecatedUnitAction = DeprecatedUnitAction.WARN,
    ) -> str:
        return super().to_string(unit, fraction=fraction)

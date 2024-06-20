# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "Console" unit format.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from . import base, utils

if TYPE_CHECKING:
    from numbers import Real
    from typing import ClassVar, Literal

    from astropy.units import UnitBase


class Console(base.Base):
    """
    Output-only format for to display pretty formatting at the
    console.

    For example::

      >>> import astropy.units as u
      >>> print(u.Ry.decompose().to_string('console'))  # doctest: +FLOAT_CMP
      2.1798721*10^-18 m^2 kg s^-2
      >>> print(u.Ry.decompose().to_string('console', fraction='multiline'))  # doctest: +FLOAT_CMP
                       m^2 kg
      2.1798721*10^-18 ------
                        s^2
      >>> print(u.Ry.decompose().to_string('console', fraction='inline'))  # doctest: +FLOAT_CMP
      2.1798721*10^-18 m^2 kg / s^2
    """

    _times: ClassVar[str] = "*"
    _line: ClassVar[str] = "-"
    _space: ClassVar[str] = " "

    @classmethod
    def _format_mantissa(cls, m: str) -> str:
        return m

    @classmethod
    def _format_superscript(cls, number: str) -> str:
        return f"^{number}"

    @classmethod
    def format_exponential_notation(cls, val: Real, format_spec: str = ".8g") -> str:
        m, ex = utils.split_mantissa_exponent(val, format_spec)

        parts = []
        if m:
            parts.append(cls._format_mantissa(m))

        if ex:
            parts.append(f"10{cls._format_superscript(ex)}")

        return cls._times.join(parts)

    @classmethod
    def _format_fraction(
        cls,
        scale: str,
        numerator: str,
        denominator: str,
        *,
        fraction: Literal[True, "inline", "multiline"] = "multiline",
    ) -> str:
        if fraction != "multiline":
            return super()._format_fraction(
                scale, numerator, denominator, fraction=fraction
            )

        fraclength = max(len(numerator), len(denominator))
        f = f"{{0:<{len(scale)}s}}{{1:^{fraclength}s}}"

        return "\n".join(
            (
                f.format("", numerator),
                f.format(scale, cls._line * fraclength),
                f.format("", denominator),
            )
        )

    @classmethod
    def to_string(
        cls, unit: UnitBase, fraction: bool | Literal["inline", "multiline"] = False
    ) -> str:
        # Change default of fraction to False, i.e., we typeset
        # without a fraction by default.
        return super().to_string(unit, fraction=fraction)

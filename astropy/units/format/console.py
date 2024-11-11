# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "Console" unit format.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from . import base

if TYPE_CHECKING:
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

    _line: ClassVar[str] = "-"
    _space: ClassVar[str] = " "

    @classmethod
    def _format_superscript(cls, number: str) -> str:
        return f"^{number}"

    @classmethod
    def _format_multiline_fraction(
        cls, scale: str, numerator: str, denominator: str
    ) -> str:
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

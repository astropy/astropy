# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "Console" unit format.
"""


from . import base, utils


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

    _times = "*"
    _line = "-"
    _space = " "

    @classmethod
    def _format_mantissa(cls, m):
        return m

    @classmethod
    def _format_superscript(cls, number):
        return f"^{number}"

    @classmethod
    def format_exponential_notation(cls, val, format_spec=".8g"):
        m, ex = utils.split_mantissa_exponent(val, format_spec)

        parts = []
        if m:
            parts.append(cls._format_mantissa(m))

        if ex:
            parts.append(f"10{cls._format_superscript(ex)}")

        return cls._times.join(parts)

    @classmethod
    def _format_fraction(cls, scale, numerator, denominator, fraction="multiline"):
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
    def to_string(cls, unit, fraction=False):
        # Change default of fraction to False, i.e., we typeset
        # without a fraction by default.
        return super().to_string(unit, fraction=fraction)

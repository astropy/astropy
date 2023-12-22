# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "Unicode" unit format.
"""
from . import console, utils


class Unicode(console.Console):
    """
    Output-only format to display pretty formatting at the console
    using Unicode characters.

    For example::

      >>> import astropy.units as u
      >>> print(u.bar.decompose().to_string('unicode'))
      100000 kg m⁻¹ s⁻²
      >>> print(u.bar.decompose().to_string('unicode', fraction='multiline'))
              kg
      100000 ────
             m s²
      >>> print(u.bar.decompose().to_string('unicode', fraction='inline'))
      100000 kg / (m s²)
    """

    _times = "×"
    _line = "─"

    @classmethod
    def _format_mantissa(cls, m):
        return m.replace("-", "−")

    @classmethod
    def _format_unit_power(cls, unit, power=1):
        name = cls._get_unit_name(unit)
        # Check for superscript units
        if power != 1:
            if name in ("°", "e⁻", "″", "′", "ʰ"):
                name = unit.short_names[0]
            name += cls._format_superscript(utils.format_power(power))
        return name

    @classmethod
    def _format_superscript(cls, number):
        mapping = str.maketrans(
            {
                "0": "⁰",
                "1": "¹",
                "2": "²",
                "3": "³",
                "4": "⁴",
                "5": "⁵",
                "6": "⁶",
                "7": "⁷",
                "8": "⁸",
                "9": "⁹",
                "-": "⁻",
                "−": "⁻",
                # This is actually a "raised omission bracket", but it's
                # the closest thing I could find to a superscript solidus.
                "/": "⸍",
            }
        )
        return number.translate(mapping)

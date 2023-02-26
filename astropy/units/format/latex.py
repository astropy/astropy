# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "LaTeX" unit format.
"""

import re

from . import console, utils


class Latex(console.Console):
    """
    Output LaTeX to display the unit based on IAU style guidelines.

    Attempts to follow the `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_.
    """

    _space = r"\,"
    _scale_unit_separator = r"\,"
    _times = r" \times "

    @classmethod
    def _get_unit_name(cls, unit):
        # Do not use super() to help latex_inline subclass.
        name = unit.get_format_name("latex")
        if name == unit.name:
            # This doesn't escape arbitrary LaTeX strings, but it should
            # be good enough for unit names which are required to be alpha
            # + "_" anyway.
            return name.replace("_", r"\_")
        else:
            return name

    @classmethod
    def _format_mantissa(cls, m):
        return m.replace("nan", r"{\rm NaN}").replace("inf", r"\infty")

    @classmethod
    def _format_superscript(cls, number):
        return f"^{{{number}}}"

    @classmethod
    def _format_unit_power(cls, unit, power=1):
        name = cls._get_unit_name(unit)
        if power != 1:
            # If the LaTeX representation of the base unit already ends with
            # a superscript, we need to spell out the unit to avoid double
            # superscripts. For example, the logic below ensures that
            # `u.deg**2` returns `deg^{2}` instead of `{}^{\circ}^{2}`.
            if re.match(r".*\^{[^}]*}$", name):  # ends w/ superscript?
                name = unit.short_names[0]
            name += cls._format_superscript(utils.format_power(power))
        return name

    @classmethod
    def _format_fraction(cls, scale, numerator, denominator, *, fraction="multiline"):
        if fraction != "multiline":
            return super()._format_fraction(
                scale, numerator, denominator, fraction=fraction
            )

        return rf"{scale}\frac{{{numerator}}}{{{denominator}}}"

    @classmethod
    def to_string(cls, unit, fraction="multiline"):
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

    name = "latex_inline"

    @classmethod
    def to_string(cls, unit, fraction=False):
        return super().to_string(unit, fraction=fraction)

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
    def _format_unit_list(cls, units):
        out = []
        for base_, power in units:
            base_latex = cls._get_unit_name(base_)
            if power == 1:
                out.append(base_latex)
            else:
                # If the LaTeX representation of the base unit already ends with
                # a superscript, we need to spell out the unit to avoid double
                # superscripts. For example, the logic below ensures that
                # `u.deg**2` returns `deg^{2}` instead of `{}^{\circ}^{2}`.
                if re.match(r".*\^{[^}]*}$", base_latex):  # ends w/ superscript?
                    base_latex = base_.short_names[0]
                out.append(f"{base_latex}^{{{utils.format_power(power)}}}")
        return r"\,".join(out)

    @classmethod
    def _format_mantissa(cls, m):
        return m.replace("nan", r"{\rm NaN}").replace("inf", r"\infty")

    @classmethod
    def _format_superscript(cls, number):
        return f"^{{{number}}}"

    @classmethod
    def _format_fraction(cls, scale, nominator, denominator):
        return rf"{scale}\frac{{{nominator}}}{{{denominator}}}"

    @classmethod
    def to_string(cls, unit, inline=False):
        s = super().to_string(unit, inline=inline)
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
    def to_string(cls, unit, inline=True):
        return super().to_string(unit, inline)

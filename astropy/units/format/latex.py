# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "LaTeX" unit format.
"""


import numpy as np

from . import base, core, utils


class Latex(base.Base):
    """
    Output LaTeX to display the unit based on IAU style guidelines.

    Attempts to follow the `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_.
    """

    @classmethod
    def _latex_escape(cls, name):
        # This doesn't escape arbitrary LaTeX strings, but it should
        # be good enough for unit names which are required to be alpha
        # + "_" anyway.
        return name.replace('_', r'\_')

    @classmethod
    def _get_unit_name(cls, unit):
        name = unit.get_format_name('latex')
        if name == unit.name:
            return cls._latex_escape(name)
        return name

    @classmethod
    def _format_unit_list(cls, units):
        out = []
        for base, power in units:
            if power == 1:
                out.append(cls._get_unit_name(base))
            else:
                out.append('{0}^{{{1}}}'.format(
                    cls._get_unit_name(base),
                    utils.format_power(power)))
        return r'\,'.join(out)

    @classmethod
    def _format_bases(cls, unit):
        positives, negatives = utils.get_grouped_by_powers(
                unit.bases, unit.powers)

        if len(negatives):
            if len(positives):
                positives = cls._format_unit_list(positives)
            else:
                positives = '1'
            negatives = cls._format_unit_list(negatives)
            s = r'\frac{{{0}}}{{{1}}}'.format(positives, negatives)
        else:
            positives = cls._format_unit_list(positives)
            s = positives

        return s

    @classmethod
    def to_string(cls, unit):
        latex_name = None
        if hasattr(unit, '_format'):
            latex_name = unit._format.get('latex')

        if latex_name is not None:
            s = latex_name
        elif isinstance(unit, core.CompositeUnit):
            if unit.scale == 1:
                s = ''
            else:
                s = cls.format_exponential_notation(unit.scale) + r'\,'

            if len(unit.bases):
                s += cls._format_bases(unit)

        elif isinstance(unit, core.NamedUnit):
            s = cls._latex_escape(unit.name)

        return r'$\mathrm{{{0}}}$'.format(s)

    @classmethod
    def format_exponential_notation(cls, val, format_spec=".8g"):
        """
        Formats a value in exponential notation for LaTeX.

        Parameters
        ----------
        val : number
            The value to be formatted

        format_spec : str, optional
            Format used to split up mantissa and exponent

        Returns
        -------
        latex_string : str
            The value in exponential notation in a format suitable for LaTeX.
        """
        if np.isfinite(val):
            m, ex = utils.split_mantissa_exponent(val, format_spec)

            parts = []
            if m:
                parts.append(m)
            if ex:
                parts.append("10^{{{0}}}".format(ex))

            return r" \times ".join(parts)
        else:
            if np.isnan(val):
                return r'{\rm NaN}'
            elif val > 0:
                # positive infinity
                return r'\infty'
            else:
                # negative infinity
                return r'-\infty'


class LatexInline(Latex):
    """
    Output LaTeX to display the unit based on IAU style guidelines with negative
    powers.

    Attempts to follow the `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_ and the
    `ApJ and AJ style guide
    <http://journals.aas.org/authors/manuscript.html>`_.
    """
    name = 'latex_inline'

    @classmethod
    def _format_bases(cls, unit):
        return cls._format_unit_list(zip(unit.bases, unit.powers))

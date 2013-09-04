# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "LaTeX" unit format.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import base
from . import utils


class Latex(base.Base):
    """
    Output LaTeX to display the unit based on IAU style guidelines.

    Attempts to follow the `IAU Style Manual
    <http://www.iau.org/static/publications/stylemanual1989.pdf>`_.
    """
    def __init__(self):
        pass

    def _latex_escape(self, name):
        # This doesn't escape arbitrary LaTeX strings, but it should
        # be good enough for unit names which are required to be alpha
        # + "_" anyway.
        return name.replace('_', r'\_')

    def _get_unit_name(self, unit):
        name = unit.get_format_name('latex')
        if name == unit.name:
            return self._latex_escape(name)
        return name

    def _format_unit_list(self, units):
        out = []
        for base, power in units:
            if power == 1:
                out.append(self._get_unit_name(base))
            else:
                out.append('{0}^{{{1}}}'.format(
                    self._get_unit_name(base),
                    utils.format_power(power)))
        return r'\,'.join(out)

    def _format_exponential_notation(self, val):
        m, ex = utils.split_mantissa_exponent(val)

        parts = []
        if m:
            parts.append(m)
        if ex:
            parts.append("10^{{{0}}}".format(ex))

        return r" \times ".join(parts)

    def to_string(self, unit):
        from .. import core

        latex_name = None
        if hasattr(unit, '_format'):
            latex_name = unit._format.get('latex')

        if latex_name is not None:
            s = latex_name
        elif isinstance(unit, core.CompositeUnit):
            if unit.scale == 1:
                s = ''
            else:
                s = self._format_exponential_notation(unit.scale) + r'\,'

            if len(unit.bases):
                positives, negatives = utils.get_grouped_by_powers(
                    unit.bases, unit.powers)
                if len(negatives):
                    if len(positives):
                        positives = self._format_unit_list(positives)
                    else:
                        positives = '1'
                    negatives = self._format_unit_list(negatives)
                    s += r'\frac{{{0}}}{{{1}}}'.format(positives, negatives)
                else:
                    positives = self._format_unit_list(positives)
                    s += positives
        elif isinstance(unit, core.NamedUnit):
            s = self._latex_escape(unit.name)

        return r'$\mathrm{{{0}}}$'.format(s)

# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "Console" unit format.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import base
from . import utils


class Console(base.Base):
    """
    Output-only format for to display pretty formatting at the
    console.

    For example::

      >>> import astropy.units as u
      >>> print(u.Ry.decompose().to_string('console'))  # doctest: +FLOAT_CMP
                       m^2 kg
      2.1798721*10^-18 ------
                        s^2
    """
    def __init__(self):
        pass

    _times = "*"
    _line = "-"

    def _get_unit_name(self, unit):
        return unit.get_format_name('console')

    def _format_superscript(self, number):
        return '^{0}'.format(number)

    def _format_unit_list(self, units):
        out = []
        for base, power in units:
            if power == 1:
                out.append(self._get_unit_name(base))
            else:
                out.append('{0}{1}'.format(
                    self._get_unit_name(base),
                    self._format_superscript(
                            utils.format_power(power))))
        return ' '.join(out)

    def _format_exponential_notation(self, val):
        m, ex = utils.split_mantissa_exponent(val)

        parts = []
        if m:
            parts.append(m)

        if ex:
            parts.append("10{0}".format(
                self._format_superscript(ex)))

        return self._times.join(parts)

    def to_string(self, unit):
        from .. import core

        if isinstance(unit, core.CompositeUnit):
            if unit.scale == 1:
                s = ''
            else:
                s = self._format_exponential_notation(unit.scale)

            if len(unit.bases):
                positives, negatives = utils.get_grouped_by_powers(
                    unit.bases, unit.powers)
                if len(negatives):
                    if len(positives):
                        positives = self._format_unit_list(positives)
                    else:
                        positives = '1'
                    negatives = self._format_unit_list(negatives)
                    l = len(s)
                    r = max(len(positives), len(negatives))
                    f = "{{0:^{0}s}} {{1:^{1}s}}".format(l, r)

                    lines = [
                        f.format('', positives),
                        f.format(s, self._line * r),
                        f.format('', negatives)
                    ]

                    s = '\n'.join(lines)
                else:
                    positives = self._format_unit_list(positives)
                    s += positives
        elif isinstance(unit, core.NamedUnit):
            s = self._get_unit_name(unit)

        return s

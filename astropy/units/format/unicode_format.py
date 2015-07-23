# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "Unicode" unit format.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import console, utils


class Unicode(console.Console):
    """
    Output-only format to display pretty formatting at the console
    using Unicode characters.

    For example::

      >>> import astropy.units as u
      >>> print(u.Ry.decompose().to_string('unicode'))  # doctest: +FLOAT_CMP
                      m² kg
      2.1798721×10⁻¹⁸ ─────
                       s²
    """

    _times = "×"
    _line = "─"

    @classmethod
    def _get_unit_name(cls, unit):
        return unit.get_format_name('unicode')

    @classmethod
    def format_exponential_notation(cls, val):
        m, ex = utils.split_mantissa_exponent(val)

        parts = []
        if m:
            parts.append(m.replace('-', '−'))

        if ex:
            parts.append("10{0}".format(
                cls._format_superscript(ex)))

        return cls._times.join(parts)

    @classmethod
    def _format_superscript(cls, number):
        mapping = {
            '0': '⁰',
            '1': '¹',
            '2': '²',
            '3': '³',
            '4': '⁴',
            '5': '⁵',
            '6': '⁶',
            '7': '⁷',
            '8': '⁸',
            '9': '⁹',
            '-': '⁻',
            '−': '⁻',
            # This is actually a "raised omission bracket", but it's
            # the closest thing I could find to a superscript solidus.
            '/': '⸍',
            }
        output = []
        for c in number:
            output.append(mapping[c])
        return ''.join(output)

# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "Unicode" unit format.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import console


class Unicode(console.Console):
    """
    Output-only format to display pretty formatting at the console
    using Unicode characters.

    For example::

      >>> import astropy.units as u
      >>> print(u.Ry.decompose().to_string('unicode'))
                   m² kg
      2.1...×10⁻¹⁸ ─────
                    s²
    """

    def __init__(self):
        pass

    _times = "×"
    _line = "─"

    def _get_unit_name(self, unit):
        return unit.get_format_name('unicode')

    @staticmethod
    def _format_superscript(number):
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
            '-': '⁻'}
        output = []
        for c in number:
            output.append(mapping[c])
        return ''.join(output)

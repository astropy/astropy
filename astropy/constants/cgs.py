# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in cgs units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Only constants that cannot be converted directly from S.I. are defined here.

from .constant import EMConstant
from . import si

# PHYSICAL CONSTANTS

# Electron charge

for cver, cc in [(si.ASTROPYCONST20, si.CODATA2014),
                    (si.ASTROPYCONST13, si.CODATA2010)]:
    cver.e_esu = EMConstant(cc.e.abbrev, cc.e.name, cc.e.value * cc.c.value * 10.0,
                   'statC', cc.e.uncertainty * cc.c.value * 10.0,
                   cc.e.reference, system='esu')
    cver.e_emu = EMConstant(cc.e.abbrev, cc.e.name, cc.e.value / 10, 'abC',
                   cc.e.uncertainty / 10, cc.e.reference, system='emu')
    cver.e_gauss = EMConstant(cc.e.abbrev, cc.e.name, cc.e.value * cc.c.value * 10.0,
                     'Fr', cc.e.uncertainty * cc.c.value * 10.0,
                     cc.e.reference, system='gauss')

e_esu = si.ASTROPYCONST20.e_esu

e_emu = si.ASTROPYCONST20.e_emu

e_gauss = si.ASTROPYCONST20.e_gauss

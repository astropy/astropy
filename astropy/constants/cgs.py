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
e_esu = EMConstant(si.e.abbrev, si.e.name, si.e.value * si.c.value * 10.0,
                   'statC', si.e.uncertainty * si.c.value * 10.0,
                   si.e.reference, system='esu')


e_emu = EMConstant(si.e.abbrev, si.e.name, si.e.value / 10, 'abC',
                   si.e.uncertainty / 10, si.e.reference, system='emu')


e_gauss = EMConstant(si.e.abbrev, si.e.name, si.e.value * si.c.value * 10.0,
                     'Fr', si.e.uncertainty * si.c.value * 10.0,
                     si.e.reference, system='gauss')

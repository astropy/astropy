# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in cgs units.  The constants
available (with approximate values) are:
"""
# This docstring is extended by __init__.py

# The values of constants are defined here instead of `astropy/constants` in
# order to avoid circular dependencies, since Constant depends on Quantity,
# Quantity depends on Unit, and Unit depends on the value of some of the
# constants, so they need to be kept separately.

# Only constants that cannot be converted directly from S.I. are defined here.

from .constant import EMConstant
from . import si

# PHYSICAL CONSTANTS

# Electron charge

e_esu = EMConstant(si.e.abbrev, si.e.name, si.e.value * si.c.value * 10.0,
                   'Fr', si.e.uncertainty * si.c.value * 10.0, si.e.reference,
                   system='esu')

e_gauss = EMConstant(si.e.abbrev, si.e.name, si.e.value * si.c.value * 10.0,
                     'statC', si.e.uncertainty * si.c.value * 10.0,
                     si.e.reference, system='gauss')

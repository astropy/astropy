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

from .definition import ConstantDefinition
from . import si

# PHYSICAL CONSTANTS

# Electron charge

e_esu = ConstantDefinition(si.e * si.c * 10.,
                           si.e.uncertainty * si.c * 10.,
                           si.e.name, si.e.reference, 'statC', system='esu')

e_gauss = ConstantDefinition(si.e * si.c * 10.,
                             si.e.uncertainty * si.c * 10.,
                             si.e.name, si.e.reference, 'Fr', system='gauss')

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from .constant import Constant

# ASTRONOMICAL CONSTANTS


class IAU2012(Constant):
    default_reference = 'IAU 2012'
    _registry = {}
    _has_incompatible_units = set()


# DISTANCE

# Astronomical Unit
au = IAU2012('au', "Astronomical Unit", 1.49597870700e11, 'm', 0.0,
              "IAU 2012 Resolution B2", system='si')

# Parsec

pc = IAU2012('pc', "Parsec", au.value / np.tan(np.radians(1. / 3600.)), 'm',
              au.uncertainty / np.tan(np.radians(1. / 3600.)),
              "Derived from au", system='si')

# Kiloparsec
kpc = IAU2012('kpc', "Kiloparsec",
               1000. * au.value / np.tan(np.radians(1. / 3600.)), 'm',
               1000. * au.uncertainty / np.tan(np.radians(1. / 3600.)),
               "Derived from au", system='si')

# Luminosity
L_bol0 = IAU2012('L_bol0', "Luminosity for absolute bolometric magnitude 0",
                  3.0128e28, "W", 0.0, "IAU 2015 Resolution B 2", system='si')


# SOLAR QUANTITIES

# Solar luminosity
L_sun = IAU2012('L_sun', "Solar luminosity", 3.846e26, 'W', 0.0005e26,
                 "Allen's Astrophysical Quantities 4th Ed.", system='si')

# Solar mass
M_sun = IAU2012('M_sun', "Solar mass", 1.9891e30, 'kg', 0.00005e30,
                 "Allen's Astrophysical Quantities 4th Ed.", system='si')

# Solar radius
R_sun = IAU2012('R_sun', "Solar radius", 6.95508e8, 'm', 0.00026e8,
                 "Allen's Astrophysical Quantities 4th Ed.", system='si')


# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = IAU2012('M_jup', "Jupiter mass", 1.8987e27, 'kg', 0.00005e27,
                 "Allen's Astrophysical Quantities 4th Ed.", system='si')

# Jupiter equatorial radius
R_jup = IAU2012('R_jup', "Jupiter equatorial radius", 7.1492e7, 'm',
                 0.00005e7, "Allen's Astrophysical Quantities 4th Ed.",
                 system='si')

# Earth mass
M_earth = IAU2012('M_earth', "Earth mass", 5.9742e24, 'kg', 0.00005e24,
                   "Allen's Astrophysical Quantities 4th Ed.", system='si')

# Earth equatorial radius
R_earth = IAU2012('R_earth', "Earth equatorial radius", 6.378136e6, 'm',
                   0.0000005e6, "Allen's Astrophysical Quantities 4th Ed.",
                   system='si')

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

import numpy as np

from .constant import Constant


class IERS2010(Constant):
    default_reference = 'IERS 2010'
    _registry = {}
    _has_incompatible_units = set()


# NATURAL DEFINING CONSTANTS
c = IERS2010('c', "Speed of light in vacuum", 299792458.,
             'm / (s)', 0.0, system='si')

# AUXILIARY DEFINING CONSTANTS
k = IERS2010('k', "Gaussian gravitational constant", 1.720209895e-2,
             '', 0.0, system='si')

L_G = IERS2010('L_G', "1−d(TT)/d(TCG)", 6.969290134e-10,
               '', 0.0, system='si')

L_B = IERS2010('L_B', "1−d(TDB)/d(TCB)", 1.550519768e-8,
               '', 0.0, system='si')

TDB0 = IERS2010('TDB0', "TDB-TCB at JD 2443144.5 TAI", -6.55e-5,
                's', 0.0, system='si')

theta0 = IERS2010('theta0', "Earth Rotation Angle (ERA) at J2000.0",
                  0.7790572732640 * 2 * np.pi, 'rad', 0.0, system='si')

dtheta = IERS2010('dtheta', "Rate of advance of ERA",
                  1.00273781191135448 * 2 * np.pi, 'rad', 0.0, system='si')

# NATURAL MEASURABLE CONSTANT
G = IERS2010('G', "Gravitational constant", 6.67428e-11,
             'm3 / (kg s2)', 0.00067e-11, system='si')

# BODY CONSTANTS
# Solar mass parameter
GM_sun = IERS2010('GM_sun', 'Nominal solar mass parameter', 1.32712442099e20,
                  'm3 / (s2)', 1e10, system='si')

# Solar dynamical form factor
J2_sun = IERS2010('J2_sun', 'Dynamical form factor of the Sun', 2e-7,
                  '', 0.0, system='si')

# Moon-Earth mass ratio
mu = IERS2010('mu', 'Moon-Earth mass ratio', 0.0123000371,
              '', 4e-10, system='si')

# EARTH CONSTANTS

# Earth mass parameter
GM_earth = IERS2010('GM_earth', 'Nominal Earth mass parameter', 3.986004418e14,
                    'm3 / (s2)', 8e5, system='si')

# Earth equatorial radius
R_earth = IERS2010('R_earth', "Nominal Earth equatorial radius", 6378136.6,
                   'm', 0.1, system='si')

# Earth dynamical form factor
J2_earth = IERS2010('J2_earth', 'Dynamical form factor of the Earth', 1.0826359e-3,
                    '', 1e-10, system='si')

# Earth flattening factor
rf_earth = IERS2010('rf_earth', 'Flattening factor of the Earth', 298.25642,
                    '', 0.00001, system='si')

# Mean equatorial gravity
ge = IERS2010('ge', 'Mean equatorial gravity', 9.7803278,
              'm / s2', 1e-6, system='si')

# Potential of the geoid
W0 = IERS2010('W0', 'Potential of the geoid', 62636856.0,
              'm2 / s2', 0.5, system='si')

# Geoppotential scale factor
R0_earth = IERS2010('R0_earth', "Geopotential scale factor",
                    GM_earth.value / W0.value, 'm', 0.1, system='si')

# Dynamical flattening
H_earth = IERS2010('H_earth', "Dynamical flattening", 3273795e-9,
                   'm', 1e-9, system='si')

# INITIAL VALUE AT J2000.0
eps0 = IERS2010('eps0', "Obliquity of the ecliptic at J2000.0",
                np.radians(84381.406 / 3600), 'rad', np.radians(0.001 / 3600),
                system='si')

# OTHER CONSTANTS
# Astronomical Unit
au = IERS2010('au', "Astronomical Unit", 1.49597870700e11, 'm', 3.0,
        system='si')

L_C = IERS2010('L_C', "Average value of 1−d(TCG)/d(TCB)", 1.48082686741e-8,
               '', 2e-7, system='si')

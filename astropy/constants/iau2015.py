# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

import numpy as np

from .config import codata
from .constant import Constant

# ASTRONOMICAL CONSTANTS


class IAU2015(Constant):
    default_reference = "IAU 2015"
    _registry = {}
    _has_incompatible_units = set()


# DISTANCE

# Astronomical Unit (did not change from 2012)
au = IAU2015(
    "au",
    "Astronomical Unit",
    1.49597870700e11,
    "m",
    0.0,
    "IAU 2012 Resolution B2",
    system="si",
)

# Parsec

pc = IAU2015(
    "pc",
    "Parsec",
    au.value / np.radians(1.0 / 3600.0),
    "m",
    au.uncertainty / np.radians(1.0 / 3600.0),
    "Derived from au + IAU 2015 Resolution B 2 note [4]",
    system="si",
)

# Kiloparsec
kpc = IAU2015(
    "kpc",
    "Kiloparsec",
    1000.0 * au.value / np.radians(1.0 / 3600.0),
    "m",
    1000.0 * au.uncertainty / np.radians(1.0 / 3600.0),
    "Derived from au + IAU 2015 Resolution B 2 note [4]",
    system="si",
)

# Luminosity
L_bol0 = IAU2015(
    "L_bol0",
    "Luminosity for absolute bolometric magnitude 0",
    3.0128e28,
    "W",
    0.0,
    "IAU 2015 Resolution B 2",
    system="si",
)


# SOLAR QUANTITIES

# Solar luminosity
L_sun = IAU2015(
    "L_sun",
    "Nominal solar luminosity",
    3.828e26,
    "W",
    0.0,
    "IAU 2015 Resolution B 3",
    system="si",
)

# Solar mass parameter
GM_sun = IAU2015(
    "GM_sun",
    "Nominal solar mass parameter",
    1.3271244e20,
    "m3 / (s2)",
    0.0,
    "IAU 2015 Resolution B 3",
    system="si",
)

# Solar mass (derived from mass parameter and gravitational constant)
M_sun = IAU2015(
    "M_sun",
    "Solar mass",
    GM_sun.value / codata.G.value,
    "kg",
    ((codata.G.uncertainty / codata.G.value) * (GM_sun.value / codata.G.value)),
    f"IAU 2015 Resolution B 3 + {codata.G.reference}",
    system="si",
)

# Solar radius
R_sun = IAU2015(
    "R_sun",
    "Nominal solar radius",
    6.957e8,
    "m",
    0.0,
    "IAU 2015 Resolution B 3",
    system="si",
)


# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass parameter
GM_jup = IAU2015(
    "GM_jup",
    "Nominal Jupiter mass parameter",
    1.2668653e17,
    "m3 / (s2)",
    0.0,
    "IAU 2015 Resolution B 3",
    system="si",
)

# Jupiter mass (derived from mass parameter and gravitational constant)
M_jup = IAU2015(
    "M_jup",
    "Jupiter mass",
    GM_jup.value / codata.G.value,
    "kg",
    ((codata.G.uncertainty / codata.G.value) * (GM_jup.value / codata.G.value)),
    f"IAU 2015 Resolution B 3 + {codata.G.reference}",
    system="si",
)

# Jupiter equatorial radius
R_jup = IAU2015(
    "R_jup",
    "Nominal Jupiter equatorial radius",
    7.1492e7,
    "m",
    0.0,
    "IAU 2015 Resolution B 3",
    system="si",
)

# Earth mass parameter
GM_earth = IAU2015(
    "GM_earth",
    "Nominal Earth mass parameter",
    3.986004e14,
    "m3 / (s2)",
    0.0,
    "IAU 2015 Resolution B 3",
    system="si",
)

# Earth mass (derived from mass parameter and gravitational constant)
M_earth = IAU2015(
    "M_earth",
    "Earth mass",
    GM_earth.value / codata.G.value,
    "kg",
    ((codata.G.uncertainty / codata.G.value) * (GM_earth.value / codata.G.value)),
    f"IAU 2015 Resolution B 3 + {codata.G.reference}",
    system="si",
)

# Earth equatorial radius
R_earth = IAU2015(
    "R_earth",
    "Nominal Earth equatorial radius",
    6.3781e6,
    "m",
    0.0,
    "IAU 2015 Resolution B 3",
    system="si",
)

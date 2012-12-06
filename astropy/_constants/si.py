# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units. The constants
available (with approximate values) are:
"""
# This docstring is extended by __init__.py

# The values of constants are defined here instead of `astropy/constants` in
# order to avoid circular dependencies, since Constant depends on Quantity,
# Quantity depends on Unit, and Unit depends on the value of some of the
# units, so they need to be kept separately.

import numpy as np

from .definition import ConstantDefinition

# PHYSICAL CONSTANTS

# Planck constant
h = ConstantDefinition(6.62606957e-34, 0.00000029e-34,
                       "Planck constant",
                       'CODATA 2010', 'J.s')

# Reduced Planck constant
hbar = ConstantDefinition(h * 0.5 / np.pi, h.uncertainty * 0.5 / np.pi,
                          "Reduced Planck constant",
                          'CODATA 2010', 'J.s')

# Boltzmann constant
k_B = ConstantDefinition(1.3806488e-23, 0.0000013e-23,
                         "Boltzmann constant",
                         'CODATA 2010', 'J/K')

# Speed of light
c = ConstantDefinition(2.99792458e8, 0.,
                       "Speed of light in vacuum",
                       'CODATA 2010', 'm/s')

# Gravitional constant
G = ConstantDefinition(6.67384e-11, 0.00080e-11,
                       "Gravitational constant",
                       'CODATA 2010', 'm3.kg-1.s-2')

# Proton mass
m_p = ConstantDefinition(1.672621777e-27, 0.000000074e-27,
                         "Proton mass",
                         'CODATA 2010', 'kg')

# Neutron mass
m_n = ConstantDefinition(1.674927351e-27, 0.000000074e-27,
                         "Neutron mass",
                         'CODATA 2010', 'kg')

# Electron mass
m_e = ConstantDefinition(9.10938291e-31, 0.00000040e-31,
                         "Electron mass",
                         'CODATA 2010', 'kg')

# Stefan-Boltzmann constant
sigma_sb = ConstantDefinition(5.670373e-8, 0.000021e-8,
                              "Stefan-Boltzmann constant",
                              'CODATA 2010', 'W.m-2.K-4')

# Electron charge
e = ConstantDefinition(1.602176565e-19, 0.000000035e-19,
                       "Electron charge",
                       'CODATA 2010', 'C')

# Avogadro's number
N_A = ConstantDefinition(6.02214129e23, 0.00000027e23,
                         "Avogadro's number",
                         'CODATA 2010', 'mol-1')

# Gas constant
R = ConstantDefinition(8.3144621, 0.0000075,
                       "Gas constant",
                       'CODATA 2010', 'J.mol-1.K-1')

# Rydberg constant
Ryd = ConstantDefinition(10973731.568539, 0.000055,
                         'Rydberg constant',
                         'CODATA 2010', 'm-1')

# DISTANCE

# Astronomical Unit
au = ConstantDefinition(1.49597870700e11, 0.0,
                        "Astronomical Unit",
                        "IAU 2012 Resolution B2", 'm')

# Parsec

pc = ConstantDefinition(au / np.tan(np.radians(1. / 3600.)),
                        au.uncertainty / np.tan(np.radians(1. / 3600.)),
                        "Parsec",
                        "Derived from au", 'm')

# Kiloparsec
kpc = ConstantDefinition(1000. * au / np.tan(np.radians(1. / 3600.)),
                         1000. * au.uncertainty / np.tan(np.radians(1. / 3600.)),
                         "Kiloparsec",
                         "Derived from au", 'm')

# SOLAR QUANTITIES

# Solar luminosity
L_sun = ConstantDefinition(3.846e26, 0.0005e26,
                           "Solar luminosity",
                           "Allen's Astrophysical Quantities 4th Ed.", 'W')

# Solar mass
M_sun = ConstantDefinition(1.9891e30, 0.00005e30,
                           "Solar mass",
                           "Allen's Astrophysical Quantities 4th Ed.", 'kg')

# Solar radius
R_sun = ConstantDefinition(6.95508e8, 0.00026e8,
                           "Solar radius",
                           "Allen's Astrophysical Quantities 4th Ed.", 'm')


# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = ConstantDefinition(1.8987e27, 0.00005e27,
                           "Jupiter mass",
                           "Allen's Astrophysical Quantities 4th Ed.", 'kg')

# Jupiter equatorial radius
R_jup = ConstantDefinition(7.1492e7, 0.00005e7,
                           "Jupiter equatorial radius",
                           "Allen's Astrophysical Quantities 4th Ed.", 'm')

# Earth mass
M_earth = ConstantDefinition(5.9742e24, 0.00005e24,
                             "Earth mass",
                             "Allen's Astrophysical Quantities 4th Ed.", 'kg')

# Earth equatorial radius
R_earth = ConstantDefinition(6.378136e6, 0.0000005e6,
                             "Earth equatorial radius",
                             "Allen's Astrophysical Quantities 4th Ed.", 'm')

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

from .definition import ConstantDefinition
from . import si

# PHYSICAL CONSTANTS

# Planck constant
h = ConstantDefinition(si.h * 1.e7, si.h.uncertainty * 1.e7,
                       si.h.name, si.h.reference, 'erg.s')

# Reduced Planck constant
hbar = ConstantDefinition(si.hbar * 1.e7, si.hbar.uncertainty * 1.e7,
                          si.hbar.name, si.hbar.reference, 'erg.s')

# Boltzmann constant [CODATA]
k_B = ConstantDefinition(si.k_B * 1.e7, si.k_B.uncertainty * 1.e7,
                         si.k_B.name, si.k_B.reference, 'erg/K')

# Speed of light [CODATA]
c = ConstantDefinition(si.c * 1.e2, si.c.uncertainty * 1.e2,
                       si.c.name, si.c.reference, 'cm/s')

# Gravitional constant [CODATA]
G = ConstantDefinition(si.G * 1.e3, si.G.uncertainty * 1.e3,
                       si.G.name, si.G.reference, 'cm3.g-1.s-2')

# Proton mass [CODATA]
m_p = ConstantDefinition(si.m_p * 1.e3, si.m_p.uncertainty * 1.e3,
                         si.m_p.name, si.m_p.reference, 'g')

# Neutron mass [CODATA]
m_n = ConstantDefinition(si.m_n * 1.e3, si.m_n.uncertainty * 1.e3,
                         si.m_n.name, si.m_n.reference, 'g')

# Electron mass [CODATA]
m_e = ConstantDefinition(si.m_e * 1.e3, si.m_e.uncertainty * 1.e3,
                         si.m_e.name, si.m_e.reference, 'g')

# Stefan-Boltzmann constant
sigma_sb = ConstantDefinition(si.sigma_sb * 1.e3, si.sigma_sb.uncertainty * 1.e3,
                              si.sigma_sb.name, si.sigma_sb.reference,
                              'erg.cm-2.K-4.s-1')

# Electron charge
e = ConstantDefinition(si.e * si.c * 10., si.e.uncertainty * si.c * 10.,
                       si.e.name, si.e.reference, 'statC')

# Avogadro's number
N_A = ConstantDefinition(si.N_A, si.N_A.uncertainty,
                         si.N_A.name, si.N_A.reference, 'mol-1')

# Gas constant
R = ConstantDefinition(si.R * 1.e7, si.R.uncertainty * 1.e7,
                       si.R.name, si.R.reference, 'erg.K-1.mol-1')

# Rydberg constant
Ryd = ConstantDefinition(si.Ryd * 1.e-2, si.Ryd.uncertainty * 1.e-2,
                         si.Ryd.name, si.Ryd.reference, 'cm-1')

# DISTANCE

# Astronomical Unit
au = ConstantDefinition(si.au * 1.e2, si.au.uncertainty * 1.e2,
                        si.au.name, si.au.reference, 'cm')

# Parsec
pc = ConstantDefinition(si.pc * 1.e2, si.pc.uncertainty * 1.e2,
                        si.pc.name, si.pc.reference, 'cm')

# Kiloparsec
kpc = ConstantDefinition(si.kpc * 1.e2, si.kpc.uncertainty * 1.e2,
                         si.kpc.name, si.kpc.reference, 'cm')

# SOLAR QUANTITIES

# Solar luminosity
L_sun = ConstantDefinition(si.L_sun * 1.e7, si.L_sun.uncertainty * 1.e7,
                           si.L_sun.name, si.L_sun.reference, 'erg/s')

# Solar mass
M_sun = ConstantDefinition(si.M_sun * 1.e3, si.M_sun.uncertainty * 1.e3,
                           si.M_sun.name, si.M_sun.reference, 'g')

# Solar radius
R_sun = ConstantDefinition(si.R_sun * 1.e2, si.R_sun.uncertainty * 1.e2,
                           si.R_sun.name, si.R_sun.reference, 'cm')

# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = ConstantDefinition(si.M_jup * 1.e3, si.M_jup.uncertainty * 1.e3,
                           si.M_jup.name, si.M_jup.reference, 'g')

# Jupiter equatorial radius
R_jup = ConstantDefinition(si.R_jup * 1.e2, si.R_jup.uncertainty * 1.e2,
                           si.R_jup.name, si.R_jup.reference, 'cm')

# Earth mass
M_earth = ConstantDefinition(si.M_earth * 1.e3, si.M_earth.uncertainty * 1.e3,
                             si.M_earth.name, si.M_earth.reference, 'g')

# Earth equatorial radius
R_earth = ConstantDefinition(si.R_earth * 1.e2, si.R_earth.uncertainty * 1.e2,
                             si.R_earth.name, si.R_earth.reference, 'cm')

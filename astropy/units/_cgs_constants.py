# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in cgs units.  The constants
available (with approximate values) are:
"""
# This docstring is extended by __init__.py

from ._constant_definition import ConstantDefinition
from . import _si_constants as si

# PHYSICAL CONSTANTS

# Planck constant
h = ConstantDefinition(si.h * 1.e7, si.h.error * 1.e7,
             si.h.name, si.h.origin, 'erg.s')

# Reduced Planck constant
hbar = ConstantDefinition(si.hbar * 1.e7, si.hbar.error * 1.e7,
                si.hbar.name, si.hbar.origin, 'erg.s')

# Boltzmann constant [CODATA]
k_B = ConstantDefinition(si.k_B * 1.e7, si.k_B.error * 1.e7,
               si.k_B.name, si.k_B.origin, 'erg.s')

# Speed of light [CODATA]
c = ConstantDefinition(si.c * 1.e2, si.c.error * 1.e2,
             si.c.name, si.c.origin, 'cm/s')

# Gravitional constant [CODATA]
G = ConstantDefinition(si.G * 1.e3, si.G.error * 1.e3,
             si.G.name, si.G.origin, 'cm3.g-1.s-2')

# Proton mass [CODATA]
m_p = ConstantDefinition(si.m_p * 1.e3, si.m_p.error * 1.e3,
               si.m_p.name, si.m_p.origin, 'g')

# Neutron mass [CODATA]
m_n = ConstantDefinition(si.m_n * 1.e3, si.m_n.error * 1.e3,
               si.m_n.name, si.m_n.origin, 'g')

# Electron mass [CODATA]
m_e = ConstantDefinition(si.m_e * 1.e3, si.m_e.error * 1.e3,
               si.m_e.name, si.m_e.origin, 'g')

# Stefan-Boltzmann constant
sigma_sb = ConstantDefinition(si.sigma_sb * 1.e3, si.sigma_sb.error * 1.e3,
                    si.sigma_sb.name, si.sigma_sb.origin,
                    'erg.cm-2.K-4.s-1')

# Electron charge
e = ConstantDefinition(si.e * si.c * 10., si.e.error * si.c * 10.,
             si.e.name, si.e.origin, 'statC')

# Avogadro's number
N_A = ConstantDefinition(si.N_A, si.N_A.error,
               si.N_A.name, si.N_A.origin, 'mol-1')

# Gas constant
R = ConstantDefinition(si.R * 1.e7, si.R.error * 1.e7,
             si.R.name, si.R.origin, 'erg.K-1.mol-1')

# Rydberg constant
Ryd = ConstantDefinition(si.Ryd * 1.e-2, si.Ryd.error * 1.e-2,
               si.Ryd.name, si.Ryd.origin, 'cm-1')

# DISTANCE

# Astronomical Unit
au = ConstantDefinition(si.au * 1.e2, si.au.error * 1.e2,
              si.au.name, si.au.origin, 'cm')

# Parsec
pc = ConstantDefinition(si.pc * 1.e2, si.pc.error * 1.e2,
              si.pc.name, si.pc.origin, 'cm')

# Kiloparsec
kpc = ConstantDefinition(si.kpc * 1.e2, si.kpc.error * 1.e2,
               si.kpc.name, si.kpc.origin, 'cm')

# SOLAR QUANTITIES

# Solar luminosity
L_sun = ConstantDefinition(si.L_sun * 1.e7, si.L_sun.error * 1.e7,
                 si.L_sun.name, si.L_sun.origin, 'erg/s')

# Solar mass
M_sun = ConstantDefinition(si.M_sun * 1.e3, si.M_sun.error * 1.e3,
                 si.M_sun.name, si.M_sun.origin, 'g')

# Solar radius
R_sun = ConstantDefinition(si.R_sun * 1.e2, si.R_sun.error * 1.e2,
                 si.R_sun.name, si.R_sun.origin, 'cm')

# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = ConstantDefinition(si.M_jup * 1.e3, si.M_jup.error * 1.e3,
                 si.M_jup.name, si.M_jup.origin, 'g')

# Jupiter equatorial radius
R_jup = ConstantDefinition(si.R_jup * 1.e2, si.R_jup.error * 1.e2,
                 si.R_jup.name, si.R_jup.origin, 'cm')

# Earth mass
M_earth = ConstantDefinition(si.M_earth * 1.e3, si.M_earth.error * 1.e3,
                 si.M_earth.name, si.M_earth.origin, 'g')

# Earth equatorial radius
R_earth = ConstantDefinition(si.R_earth * 1.e2, si.R_earth.error * 1.e2,
                 si.R_earth.name, si.R_earth.origin, 'cm')

# Numerical constants, in S.I. units

from numpy import pi

import cgs

# PHYSICAL CONSTANTS

# Planck constant [CODATA]
h = cgs.h * 1.e-7  # J.s

# Reduced Planck constant
hbar = cgs.hbar * 1.e-7  # J.s

# Boltzmann constant [CODATA]
k_B = cgs.k_B * 1.e-7  # J/K

# Speed of light [CODATA]
c = cgs.c * 1.e-2  # m/s

# Gravitional constant [CODATA]
G = cgs.G * 1.e-3  # m^3/kg/s

# Proton mass [CODATA]
m_p = cgs.m_p * 1.e-3  # kg

# Neutron mass [CODATA]
m_n = cgs.m_n * 1.e-3  # kg

# Electron mass [CODATA]
m_e = cgs.m_e * 1.e-3  # kg

# Stefan-Boltzmann constant
stef_boltz = 2. * pi ** 5 * k_B ** 4 / 15. / h ** 3. / c ** 2.

# Electron charge
e = cgs.e / 2997924580.  # C

# Avogadro's number
N_A = cgs.N_A

# Gas constant
R = N_A * k_B

# DISTANCE

# Astronomical Unit
au = cgs.au * 1.e-2  # m

# Parsec
pc = 3600. * 180. / pi * au  # m

# Kiloparsec
kpc = 1000. * pc  # m

# TIME

# Year
year = cgs.year

# SOLAR QUANTITIES

# Solar luminosity
L_sun = cgs.L_sun * 1.e-7  # W

# Solar mass
M_sun = cgs.M_sun * 1.e-3  # kg

# Solar radius
R_sun = cgs.R_sun * 1.e-2  # m

# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = cgs.M_jup * 1.e-3  # kg

# Jupiter radius
R_jup = cgs.R_jup * 1.e-2  # m

# Earth mass
M_earth = cgs.M_earth * 1.e-3  # kg

# Earth radius
R_earth = cgs.R_earth * 1.e-2  # m

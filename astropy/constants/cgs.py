# Numerical constants, in cgs units

from numpy import pi

# PHYSICAL CONSTANTS

# Planck constant [CODATA]
h = 6.62606957e-27  # erg.s

# Reduced Planck constant
hbar = h * 0.5 / pi  # erg.s

# Boltzmann constant [CODATA]
k_B = 1.3806488e-16  # erg/K

# Speed of light [CODATA]
c = 2.99792458e10  # cm/s

# Gravitional constant [CODATA]
G = 6.67384e-8  # cm^3/g/s

# Proton mass [CODATA]
m_p = 1.672621777e-24  # g

# Neutron mass [CODATA]
m_n = 1.674927351e-24  # g

# Electron mass [CODATA]
m_e = 9.10938291e-28  # g

# Stefan-Boltzmann constant
stef_boltz = 2. * pi ** 5 * k_B ** 4 / 15. / h ** 3. / c ** 2.

# Electron charge
e = 4.8032068e-10  # statcoulombs

# Avogadro's number
N_A = 6.0221367e23

# Gas constant
R = N_A * k_B

# DISTANCE

# Astronomical Unit
au = 1.49598e13  # cm

# Parsec
pc = 3600. * 180. / pi * au  # cm

# Kiloparsec
kpc = 1000. * pc  # cm

# TIME

# Year
year = 31556925.9936

# SOLAR QUANTITIES

# Solar luminosity
L_sun = 3.846e33  # erg/s

# Solar mass
M_sun = 1.989e33  # g

# Solar radius
R_sun = 6.95508e10  # cm

# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = 1.8986e30  # g

# Jupiter radius
R_jup = 6.9911e9  # cm

# Earth mass
M_earth = 5.9722e27  # g

# Earth radius
R_earth = 6.371e8  # cm

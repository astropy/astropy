# Numerical constants, in cgs units

from numpy import pi

from .constant import Constant
from . import si

# PHYSICAL CONSTANTS

# Planck constant
h = Constant(si.h * 1.e7, si.h.error * 1.e7,
             si.h.name, si.h.origin, 'cgs')  # erg.s

# Reduced Planck constant
hbar = Constant(si.hbar * 1.e7, si.hbar.error * 1.e7,
                si.hbar.name, si.hbar.origin, 'cgs')  # erg.s

# Boltzmann constant [CODATA]
k_B = Constant(si.k_B * 1.e7, si.k_B.error * 1.e7,
               si.k_B.name, si.k_B.origin, 'cgs')  # erg.s

# Speed of light [CODATA]
c = Constant(si.c * 1.e2, si.c.error * 1.e2,
             si.c.name, si.c.origin, 'cgs')  # cm/s

# Gravitional constant [CODATA]
G = Constant(si.G * 1.e3, si.G.error * 1.e3,
             si.G.name, si.G.origin, 'cgs')  # cm^3/g/s

# Proton mass [CODATA]
m_p = Constant(si.m_p * 1.e3, si.m_p.error * 1.e3,
               si.m_p.name, si.m_p.origin, 'cgs')  # g

# Neutron mass [CODATA]
m_n = Constant(si.m_n * 1.e3, si.m_n.error * 1.e3,
               si.m_n.name, si.m_n.origin, 'cgs')  # g

# Electron mass [CODATA]
m_e = Constant(si.m_e * 1.e3, si.m_e.error * 1.e3,
               si.m_e.name, si.m_e.origin, 'cgs')  # g

# Stefan-Boltzmann constant
stef_boltz = Constant(si.stef_boltz * 1.e5, si.stef_boltz.error * 1.e5,
                      si.stef_boltz.name, si.stef_boltz.origin,
                      'cgs')  # erg/cm^2/K^4/s

# Electron charge
e = Constant(si.e * si.c * 10., si.e.error * si.c * 10.,
             si.e.name, si.e.origin, 'cgs')  # statC

# Avogadro's number
N_A = Constant(si.N_A, si.N_A.error,
               si.N_A.name, si.N_A.origin, 'cgs')

# Gas constant
R = Constant(si.R * 1.e7, si.R.error * 1.e7,
             si.R.name, si.R.origin, 'cgs')  # erg/K/mol

# DISTANCE

# Astronomical Unit
au = Constant(si.au * 1.e2, si.au.error * 1.e2,
              si.au.name, si.au.origin, 'cgs')  # cm

# Parsec
pc = Constant(si.pc * 1.e2, si.pc.error * 1.e2,
              si.pc.name, si.pc.origin, 'cgs')  # cm

# Kiloparsec
kpc = Constant(si.kpc * 1.e2, si.kpc.error * 1.e2,
               si.kpc.name, si.kpc.origin, 'cgs')  # cm

# SOLAR QUANTITIES

# Solar luminosity
L_sun = Constant(si.L_sun * 1.e7, si.L_sun.error * 1.e7,
                 si.L_sun.name, si.L_sun.origin, 'cgs')  # erg/s

# Solar mass
M_sun = Constant(si.M_sun * 1.e3, si.M_sun.error * 1.e3,
                 si.M_sun.name, si.M_sun.origin, 'cgs')  # g

# Solar radius
R_sun = Constant(si.R_sun * 1.e2, si.R_sun.error * 1.e2,
                 si.R_sun.name, si.R_sun.origin, 'cgs')  # cm

# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = Constant(si.M_jup * 1.e3, si.M_jup.error * 1.e3,
                 si.M_jup.name, si.M_jup.origin, 'cgs')  # g

# Jupiter equatorial radius
R_jup = Constant(si.R_jup * 1.e2, si.R_jup.error * 1.e2,
                 si.R_jup.name, si.R_jup.origin, 'cgs')  # cm

# Earth mass
M_earth = Constant(si.M_earth * 1.e3, si.M_earth.error * 1.e3,
                 si.M_earth.name, si.M_earth.origin, 'cgs')  # g

# Earth equatorial radius
R_earth = Constant(si.R_earth * 1.e2, si.R_earth.error * 1.e2,
                 si.R_earth.name, si.R_earth.origin, 'cgs')  # cm

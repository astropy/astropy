# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in cgs units.

G        6.67384e-08     cm^3/g/s^2      Gravitational constant
Jy       1e-23           erg/s/cm^2/Hz   Jansky
L_sun    3.846e+33       erg/s           Solar luminosity
M_earth  5.9742e+27      g               Earth mass
M_jup    1.8987e+30      g               Jupiter mass
M_sun    1.9891e+33      g               Solar mass
Mpc      3.0856776e+24   cm              Megaparsec
N_A      6.02214129e+23  /mol            Avogadro's number
R        83144621        erg/K/mol       Gas constant
R_earth  637813600       cm              Earth equatorial radius
R_jup    7149200000      cm              Jupiter equatorial radius
R_sun    6.95508e+10     cm              Solar radius
Ryd      109737.3157     cm^-1           Rydberg constant
au       1.495978707e+13 cm              Astronomical Unit
c        2.99792458e+10  cm/s            Speed of light in vacuum
e        4.803204506e-10 statC           Electron charge
eV       1.602176565e-12 erg             Electron volt
h        6.62606957e-27  erg.s           Planck constant
hbar     1.054571725e-27 erg.s           Reduced Planck constant
k_B      1.3806488e-16   erg.s           Boltzmann constant
kpc      3.0856776e+21   cm              Kiloparsec
m_e      9.10938291e-28  g               Electron mass
m_n      1.674927351e-24 g               Neutron mass
m_p      1.672621777e-24 g               Proton mass
pc       3.0856776e+18   cm              Parsec
sigma_sb 5.670373e-05    erg/cm^2/K^4/s  Stefan-Boltzmann constant

This list is generated using the cgs.info() function and the values
printed are approximate. Look at the individual constants for precise
values, errors, and their source.
"""

from .constant import Constant
from . import si

def info():
    """ Return a string giving an overview of constants in this module.

    The values printed are approximate. Look at the individual
    constants for precise values, errors, and their source.
    """
    s = []
    for key in sorted(globals()):
        val = globals()[key]
        try:
            s.append('%-8s %-15.10g %-14s  %-s' % (
                key, val.real, val.units, val.name))
        except AttributeError:
            pass

    return '\n'.join(s)

# PHYSICAL CONSTANTS

# Planck constant
h = Constant(si.h * 1.e7, si.h.error * 1.e7,
             si.h.name, si.h.origin, 'erg.s')

# Reduced Planck constant
hbar = Constant(si.hbar * 1.e7, si.hbar.error * 1.e7,
                si.hbar.name, si.hbar.origin, 'erg.s')

# Boltzmann constant [CODATA]
k_B = Constant(si.k_B * 1.e7, si.k_B.error * 1.e7,
               si.k_B.name, si.k_B.origin, 'erg.s')

# Speed of light [CODATA]
c = Constant(si.c * 1.e2, si.c.error * 1.e2,
             si.c.name, si.c.origin, 'cm/s')

# Gravitional constant [CODATA]
G = Constant(si.G * 1.e3, si.G.error * 1.e3,
             si.G.name, si.G.origin, 'cm^3/g/s^2')

# Proton mass [CODATA]
m_p = Constant(si.m_p * 1.e3, si.m_p.error * 1.e3,
               si.m_p.name, si.m_p.origin, 'g')

# Neutron mass [CODATA]
m_n = Constant(si.m_n * 1.e3, si.m_n.error * 1.e3,
               si.m_n.name, si.m_n.origin, 'g')

# Electron mass [CODATA]
m_e = Constant(si.m_e * 1.e3, si.m_e.error * 1.e3,
               si.m_e.name, si.m_e.origin, 'g')

# Stefan-Boltzmann constant
sigma_sb = Constant(si.sigma_sb * 1.e3, si.sigma_sb.error * 1.e3,
                    si.sigma_sb.name, si.sigma_sb.origin,
                    'erg/cm^2/K^4/s')

# Electron charge
e = Constant(si.e * si.c * 10., si.e.error * si.c * 10.,
             si.e.name, si.e.origin, 'statC')

# Avogadro's number
N_A = Constant(si.N_A, si.N_A.error,
               si.N_A.name, si.N_A.origin, '/mol')

# Gas constant
R = Constant(si.R * 1.e7, si.R.error * 1.e7,
             si.R.name, si.R.origin, 'erg/K/mol')

# Rydberg constant
Ryd = Constant(si.Ryd * 1.e-2, si.Ryd.error * 1.e-2,
               si.Ryd.name, si.Ryd.origin, 'cm^-1')

# Electron volt
eV = Constant(si.eV * 1e7, si.eV.error * 1e7, si.eV.name,
              si.eV.origin, 'erg')

# Jansky
Jy = Constant(si.Jy * 1e3, si.Jy.error * 1e3, si.Jy.name,
              si.Jy.origin, 'erg/s/cm^2/Hz')

# DISTANCE

# Astronomical Unit
au = Constant(si.au * 1.e2, si.au.error * 1.e2,
              si.au.name, si.au.origin, 'cm')

# Parsec
pc = Constant(si.pc * 1.e2, si.pc.error * 1.e2,
              si.pc.name, si.pc.origin, 'cm')

# Kiloparsec
kpc = Constant(si.kpc * 1.e2, si.kpc.error * 1.e2,
               si.kpc.name, si.kpc.origin, 'cm')

# Megaparsec
Mpc = Constant(si.Mpc * 1.e2, si.Mpc.error * 1.e2,
               si.Mpc.name, si.Mpc.origin, 'cm')


# SOLAR QUANTITIES

# Solar luminosity
L_sun = Constant(si.L_sun * 1.e7, si.L_sun.error * 1.e7,
                 si.L_sun.name, si.L_sun.origin, 'erg/s')

# Solar mass
M_sun = Constant(si.M_sun * 1.e3, si.M_sun.error * 1.e3,
                 si.M_sun.name, si.M_sun.origin, 'g')

# Solar radius
R_sun = Constant(si.R_sun * 1.e2, si.R_sun.error * 1.e2,
                 si.R_sun.name, si.R_sun.origin, 'cm')

# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = Constant(si.M_jup * 1.e3, si.M_jup.error * 1.e3,
                 si.M_jup.name, si.M_jup.origin, 'g')

# Jupiter equatorial radius
R_jup = Constant(si.R_jup * 1.e2, si.R_jup.error * 1.e2,
                 si.R_jup.name, si.R_jup.origin, 'cm')

# Earth mass
M_earth = Constant(si.M_earth * 1.e3, si.M_earth.error * 1.e3,
                 si.M_earth.name, si.M_earth.origin, 'g')

# Earth equatorial radius
R_earth = Constant(si.R_earth * 1.e2, si.R_earth.error * 1.e2,
                 si.R_earth.name, si.R_earth.origin, 'cm')

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.

G        6.67384e-11     m^3/kg/s^2      Gravitational constant
Jy       1e-26           W/m^2/Hz        Jansky
L_sun    3.846e+26       W               Solar luminosity
M_earth  5.9742e+24      kg              Earth mass
M_jup    1.8987e+27      kg              Jupiter mass
M_sun    1.9891e+30      kg              Solar mass
Mpc      3.0856776e+22   m               Megaparsec
N_A      6.02214129e+23  /mol            Avogadro's number
R        8.3144621       J/mol/K         Gas constant
R_earth  6378136         m               Earth equatorial radius
R_jup    71492000        m               Jupiter equatorial radius
R_sun    695508000       m               Solar radius
Ryd      10973731.57     m^-1            Rydberg constant
au       1.495978707e+11 m               Astronomical Unit
c        299792458       m/s             Speed of light in vacuum
e        1.602176565e-19 C               Electron charge
eV       1.602176565e-19 J               Electron volt
h        6.62606957e-34  J.s             Planck constant
hbar     1.054571725e-34 J.s             Reduced Planck constant
k_B      1.3806488e-23   J/K             Boltzmann constant
kpc      3.0856776e+19   m               Kiloparsec
m_e      9.10938291e-31  kg              Electron mass
m_n      1.674927351e-27 kg              Neutron mass
m_p      1.672621777e-27 kg              Proton mass
pc       3.0856776e+16   m               Parsec
sigma_sb 5.670373e-08    W/m^2/K^4       Stefan-Boltzmann constant

This list is generated using the si.info() function and the values
printed are approximate. Look at the individual constants for precise
values, errors, and their source.
"""

import numpy as np
from .constant import Constant

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
h = Constant(6.62606957e-34, 0.00000029e-34,
             "Planck constant",
             'CODATA 2010', 'J.s')

# Reduced Planck constant
hbar = Constant(h * 0.5 / np.pi, h.error * 0.5 / np.pi,
                "Reduced Planck constant",
                'CODATA 2010', 'J.s')

# Boltzmann constant
k_B = Constant(1.3806488e-23, 0.0000013e-23,
               "Boltzmann constant",
               'CODATA 2010', 'J/K')

# Speed of light
c = Constant(2.99792458e8, 0.,
             "Speed of light in vacuum",
             'CODATA 2010', 'm/s')

# Gravitional constant
G = Constant(6.67384e-11, 0.00080e-11,
             "Gravitational constant",
             'CODATA 2010', 'm^3/kg/s^2')

# Proton mass
m_p = Constant(1.672621777e-27, 0.000000074e-27,
               "Proton mass",
               'CODATA 2010', 'kg')

# Neutron mass
m_n = Constant(1.674927351e-27, 0.000000074e-27,
               "Neutron mass",
               'CODATA 2010', 'kg')

# Electron mass
m_e = Constant(9.10938291e-31, 0.00000040e-31,
               "Electron mass",
               'CODATA 2010', 'kg')

# Stefan-Boltzmann constant
sigma_sb = Constant(5.670373e-8, 0.000021e-8,
                    "Stefan-Boltzmann constant",
                    'CODATA 2010', 'W/m^2/K^4')

# Electron charge
e = Constant(1.602176565e-19, 0.000000035e-19,
             "Electron charge",
             'CODATA 2010', 'C')

# Avogadro's number
N_A = Constant(6.02214129e23, 0.00000027e23,
               "Avogadro's number",
               'CODATA 2010', '/mol')

# Gas constant
R = Constant(8.3144621, 0.0000075,
             "Gas constant",
             'CODATA 2010', 'J/mol/K')

# Rydberg constant
Ryd = Constant(10973731.568539, 0.000055,
               'Rydberg constant',
               'CODATA 2010', 'm^-1')

# Electron volt
eV = Constant(1.602176565e-19, 0.000000035e-19,
              'Electron volt',
              'CODATA 2010', 'J')

# Jansky
Jy = Constant(1e-26, 0, 'Jansky',
              "Allen's Astrophysical Quantities 4th", 'W/m^2/Hz')


# DISTANCE

# Astronomical Unit
au = Constant(1.4959787066e11, 0.00000000005e11,
              "Astronomical Unit",
              "Allen's Astrophysical Quantities 4th Ed.", 'm')

# Parsec
pc = Constant(3.0856776e16, 0.00000005e16,
              "Parsec",
              "Allen's Astrophysical Quantities 4th Ed.", 'm')

# Kiloparsec
kpc = Constant(3.0856776e19, 0.00000005e19,
              "Kiloparsec",
              "Allen's Astrophysical Quantities 4th Ed.", 'm')

# Megaparsec
Mpc = Constant(3.0856776e22, 0.00000005e22,
              "Megaparsec",
              "Allen's Astrophysical Quantities 4th Ed.", 'm')

# SOLAR QUANTITIES

# Solar luminosity
L_sun = Constant(3.846e26, 0.0005e26,
                 "Solar luminosity",
                 "Allen's Astrophysical Quantities 4th Ed.", 'W')

# Solar mass
M_sun = Constant(1.9891e30, 0.00005e30,
                 "Solar mass",
                 "Allen's Astrophysical Quantities 4th Ed.", 'kg')

# Solar radius
R_sun = Constant(6.95508e8, 0.00026e8,
                 "Solar radius",
                 "Allen's Astrophysical Quantities 4th Ed.", 'm')


# OTHER SOLAR SYSTEM QUANTITIES

# Jupiter mass
M_jup = Constant(1.8987e27, 0.00005e27,
                 "Jupiter mass",
                 "Allen's Astrophysical Quantities 4th Ed.", 'kg')

# Jupiter equatorial radius
R_jup = Constant(7.1492e7, 0.00005e7,
                 "Jupiter equatorial radius",
                 "Allen's Astrophysical Quantities 4th Ed.", 'm')

# Earth mass
M_earth = Constant(5.9742e24, 0.00005e24,
                  "Earth mass",
                  "Allen's Astrophysical Quantities 4th Ed.", 'kg')

# Earth equatorial radius
R_earth = Constant(6.378136e6, 0.0000005e6,
                   "Earth equatorial radius",
                   "Allen's Astrophysical Quantities 4th Ed.", 'm')

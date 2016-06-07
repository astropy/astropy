# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from .constant import Constant, EMConstant

# PHYSICAL CONSTANTS

# Planck constant
h = Constant('h', "Planck constant", 6.626070040e-34, 'J s', 0.000000081e-34,
             'CODATA 2014', system='si')

# Reduced Planck constant
hbar = Constant('hbar', "Reduced Planck constant", h.value * 0.5 / np.pi,
                'J s', h.uncertainty * 0.5 / np.pi, h.reference, system='si')

# Boltzmann constant
k_B = Constant('k_B', "Boltzmann constant", 1.38064852e-23, 'J / (K)',
               0.00000079e-23, 'CODATA 2014', system='si')

# Speed of light
c = Constant('c', "Speed of light in vacuum", 299792458., 'm / (s)', 0.,
             'CODATA 2014', system='si')

# Gravitional constant
G = Constant('G', "Gravitational constant", 6.67408e-11, 'm3 / (kg s2)',
             0.00031e-11, 'CODATA 2014', system='si')

# Standard acceleration of gravity
g0 = Constant('g0', "Standard acceleration of gravity", 9.80665, 'm / s2', 0.0,
              'CODATA 2010', system='si')

# Proton mass
m_p = Constant('m_p', "Proton mass", 1.672621898e-27, 'kg', 0.000000021e-27,
               'CODATA 2014', system='si')

# Neutron mass
m_n = Constant('m_n', "Neutron mass", 1.674927471e-27, 'kg', 0.000000021e-27,
               'CODATA 2014', system='si')

# Electron mass
m_e = Constant('m_e', "Electron mass", 9.10938356e-31, 'kg', 0.00000011e-31,
               'CODATA 2014', system='si')

# Atomic mass
u = Constant('u', "Atomic mass", 1.660539040e-27, 'kg', 0.000000020e-27,
             'CODATA 2014', system='si')

# Stefan-Boltzmann constant
sigma_sb = Constant('sigma_sb', "Stefan-Boltzmann constant", 5.670367e-8,
                    'W / (K4 m2)', 0.000013e-8, 'CODATA 2014', system='si')

# Electron charge; EM constants require a system to be specified
e = EMConstant('e', 'Electron charge', 1.6021766208e-19, 'C', 0.0000000098e-19,
               'CODATA 2014', system='si')

# Electric constant
eps0 = EMConstant('eps0', 'Electric constant', 8.854187817e-12, 'F/m', 0.0,
                  'CODATA 2014', system='si')

# Avogadro's number
N_A = Constant('N_A', "Avogadro's number", 6.022140857, '1 / (mol)',
               0.000000074, 'CODATA 2014', system='si')

# Gas constant
R = Constant('R', "Gas constant", 8.3144598, 'J / (K mol)', 0.0000048,
             'CODATA 2014', system='si')

# Rydberg constant
Ryd = Constant('Ryd', 'Rydberg constant', 10973731.568508, '1 / (m)', 0.000065,
               'CODATA 2014', system='si')

# Bohr radius
a0 = Constant('a0', "Bohr radius", 0.52917721067e-10, 'm', 0.00000000012e-10,
              'CODATA 2014', system='si')

# Bohr magneton
muB = Constant('muB', "Bohr magneton", 927.4009994-26, 'J/T', 0.0000057e-26,
               'CODATA 2014', system='si')

# Fine structure constant
alpha = Constant('alpha', "Fine-structure constant", 7.2973525664e-3, '',
                 0.0000000017e-3, 'CODATA 2014', system='si')

# Atmosphere
atm = Constant('atmosphere', "Atmosphere", 101325, 'Pa', 0.0,
               'CODATA 2010', system='si')

# Magnetic constant
mu0 = Constant('mu0', "Magnetic constant", 4.0e-7 * np.pi, 'N/A2', 0.0,
               'CODATA 2014', system='si')

# Thomson scattering cross-section
sigma_T = Constant('sigma_T', "Thomson scattering cross-section",
                   0.66524587158e-28, 'm2', 0.00000000091e-28, 
                   'CODATA 2014', system='si')

# DISTANCE

# Astronomical Unit
au = Constant('au', "Astronomical Unit", 1.49597870700e11, 'm', 0.0,
              "IAU 2012 Resolution B2", system='si')

# Parsec

pc = Constant('pc', "Parsec", au.value / np.tan(np.radians(1. / 3600.)), 'm',
              au.uncertainty / np.tan(np.radians(1. / 3600.)),
              "Derived from au", system='si')

# Kiloparsec
kpc = Constant('kpc', "Kiloparsec",
               1000. * au.value / np.tan(np.radians(1. / 3600.)), 'm',
               1000. * au.uncertainty / np.tan(np.radians(1. / 3600.)),
               "Derived from au", system='si')

# Luminosity
L_bol0 = Constant('L_bol0', "Luminosity for absolute bolometric magnitude 0",
                  3.0128e28, "W", 0.0, "IAU 2015 Resolution B2", system='si')

# Wien wavelength displacement law constant
b_wien = Constant('b_wien', 'Wien wavelength displacement law constant',
                  2.8977729e-3, 'm K', 0.0000017e-3, 'CODATA 2014', system='si')

# SOLAR QUANTITIES

# Solar luminosity
L_Sun = Constant('L_Sun', "Solar luminosity", 3.828e26, 'W', 0.0,
                 "IAU 2015 Resolution B3", system='si')

# Solar radius
R_Sun = Constant('R_Sun', "Solar radius", 6.957e8, 'm', 0.0,
                 "IAU 2015 Resolution B3", system='si')

# Solar irradiance
S_Sun = Constant('S_Sun', "Solar irradiance", 1361., 'W', 0.0,
                 "IAU 2015 Resolution B3", system='si')

# Solar effective temperature
Teff_Sun = Constant('Teff_Sun', "Solar effective temperature", 5772., 'W', 0.0,
                 "IAU 2015 Resolution B3", system='si')

# Solar mass parameter
GM_Sun = Constant('GM_Sun', "Solar mass parameter", 1.3271244e20, 'm3/(s2)', 0.0,
                 "IAU 2015 Resolution B3", system='si')

# Solar mass 
M_Sun = Constant('M_Sun', "Solar mass", 
                 GM_Sun.value/G.value, 'kg', 0.0,
                 "Derived from GM_Sun and G", system='si')

# Solar density
rho_Sun = Constant('rho_Sun', "Solar mean density",
    M_Sun.value/(4.*np.pi*(R_Sun.value**3)/3.), 'kg/(m3)', 0.0,
                 "Derived from M_Sun and R_Sun", system='si')

# PLANETARY CONVERSION CONSTANTS

# Earth equatorial radius
Re_Earth = Constant('Re_Earth', "Earth equatorial radius", 6.3781e6, 'm',
                   0., "IAU 2015 Resolution B3", system='si')

# Earth polar radius
Rp_Earth = Constant('Rp_Earth', "Earth polar radius", 6.3568e6, 'm',
                   0., "IAU 2015 Resolution B3", system='si')

# Earth volume
V_Earth = Constant('V_Earth', "Earth volume", 
                   4*np.pi*(Re_Earth.value)**2*Rp_Earth.value/3., 'm3',
                   0., "IAU 2015 Resolution B3", system='si')

# Earth mass parameter
GM_Earth = Constant('GM_Earth', "Earth mass parameter", 3.986004e14, 'm3/(s2)',
                   0., "IAU 2015 Resolution B3", system='si')

# Earth mass 
M_Earth = Constant('M_Earth', "Earth mass", GM_Earth.value/G.value, 'kg',
                   0., "Derived from GM_Earth and G", system='si')

# Earth density
rho_Earth = Constant('rho_Earth', "Earth mean density", 
                   M_Earth.value/V_Earth.value, 'kg/(m3)',
                   0., "Derived from M_Earth and V_Earth", system='si')

# Jupiter equatorial radius
Re_Jup = Constant('Re_Jup', "Jupiter equatorial radius", 7.1492e7, 'm',
                   0., "IAU 2015 Resolution B3", system='si')

# Jupiter polar radius
Rp_Jup = Constant('Rp_Jup', "Jupiter polar radius", 6.6854e7, 'm',
                   0., "IAU 2015 Resolution B3", system='si')

# Jupiter volume
V_Jup = Constant('V_Jup', "Jupiter volume", 
                   4*np.pi*(Re_Jup.value)**2*Rp_Jup.value/3., 'm3',
                   0., "IAU 2015 Resolution B3", system='si')

# Jupiter mass parameter
GM_Jup = Constant('GM_Jup', "Jupiter mass parameter", 1.2668653e17, 'm3/(s2)',
                   0., "IAU 2015 Resolution B3", system='si')

# Jupiter mass 
M_Jup = Constant('M_Jup', "Jupiter mass", GM_Jup.value/G.value, 'kg',
                   0., "Derived from GM_Jup and G", system='si')

# Jupiter density
rho_Jup = Constant('rho_Jup', "Jupiter mean density", M_Jup.value/V_Jup.value, 
                    'kg/(m3)', 0., "Derived from M_Jup and V_Jup", system='si')


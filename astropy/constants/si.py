# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import inspect
import numpy as np

from .constant import Constant, EMConstant


def get_attributes(a_class):
    attributes = inspect.getmembers(a_class,
                               lambda a:not(inspect.isroutine(a)))
    return [a for a in attributes
                if not(a[0].startswith('__') and a[0].endswith('__'))]

# PHYSICAL CONSTANTS

class CODATA2014:
    h = Constant('h', "Planck constant", 6.626070040e-34,
                            'J s', 0.000000081e-34, 'CODATA 2014',
                            system='si')

    hbar = Constant('hbar', "Reduced Planck constant", 1.054571800e-34,
                        'J s', 0.000000013e-34, 'CODATA 2014',
                        system='si')

    k_B = Constant('k_B', "Boltzmann constant", 1.38064852e-23,
                        'J / (K)', 0.00000079e-23, 'CODATA 2014',
                        system='si')

    c = Constant('c', "Speed of light in vacuum", 299792458.,
                        'm / (s)', 0.0, 'CODATA 2014', system='si')


    G= Constant('G', "Gravitational constant", 6.67408e-11,
                        'm3 / (kg s2)', 0.00031e-11, 'CODATA 2014',
                        system='si')

    g0 = Constant('g0', "Standard acceleration of gravity", 9.80665,
                             'm / s2', 0.0, 'CODATA 2014', system='si')

    m_p = Constant('m_p', "Proton mass", 1.672621898e-27,
                   'kg', 0.000000021e-27, 'CODATA 2014', system='si')

    m_n = Constant('m_n', "Neutron mass", 1.674927471e-27,
                   'kg', 0.000000021e-27, 'CODATA 2014', system='si')

    m_e = Constant('m_e', "Electron mass", 9.10938356e-31,
                   'kg', 0.00000011e-31, 'CODATA 2014', system='si')

    u = Constant('u', "Atomic mass", 1.660539040e-27,
                 'kg', 0.000000020e-27, 'CODATA 2014', system='si')

    sigma_sb = Constant('sigma_sb', "Stefan-Boltzmann constant", 5.670367e-8,
                    'W / (K4 m2)', 0.000013e-8, 'CODATA 2014', system='si')

    e = EMConstant('e', 'Electron charge', 1.6021766208e-19,
                   'C', 0.0000000098e-19, 'CODATA 2014', system='si')

    eps0 = EMConstant('eps0', 'Electric constant', 8.854187817e-12,
                      'F/m', 0.0, 'CODATA 2014', system='si')

    N_A = Constant('N_A', "Avogadro's number", 6.022140857e23,
                   '1 / (mol)', 0.000000074e23, 'CODATA 2014', system='si')

    R = Constant('R', "Gas constant", 8.3144598,
                 'J / (K mol)', 0.0000048, 'CODATA 2014', system='si')

    Ryd = Constant('Ryd', 'Rydberg constant', 10973731.568508,
                   '1 / (m)', 0.000065, 'CODATA 2014', system='si')

    a0 = Constant('a0', "Bohr radius", 0.52917721067e-10,
                  'm', 0.00000000012e-10, 'CODATA 2014', system='si')

    muB = Constant('muB', "Bohr magneton", 927.4009994e-26,
                   'J/T', 0.00002e-26, 'CODATA 2014', system='si')

    alpha = Constant('alpha', "Fine-structure constant", 7.2973525664e-3,
                     '', 0.0000000017e-3, 'CODATA 2014', system='si')

    atm = Constant('atmosphere', "Atmosphere", 101325,
                   'Pa', 0.0, 'CODATA 2014', system='si')

    mu0 = Constant('mu0', "Magnetic constant", 4.0e-7 * np.pi, 'N/A2', 0.0,
                   'CODATA 2014', system='si')

    sigma_T = Constant('sigma_T', "Thomson scattering cross-section",
                       0.66524587158e-28, 'm2', 0.00000000091e-28,
                       'CODATA 2014', system='si')


class CODATA2010:

    h = Constant('h', "Planck constant", 6.62606957e-34, 'J s',
                        0.00000029e-34, 'CODATA 2010', system='si')

    hbar = Constant('hbar', "Reduced Planck constant",
                        h.value * 0.5 / np.pi, 'J s',
                        h.uncertainty * 0.5 / np.pi,
                        h.reference, system='si')

    k_B = Constant('k_B', "Boltzmann constant", 1.3806488e-23, 'J / (K)',
                   0.0000013e-23, 'CODATA 2010', system='si')

    c = Constant('c', "Speed of light in vacuum", 2.99792458e8, 'm / (s)', 0.,
                 'CODATA 2010', system='si')

    G = Constant('G', "Gravitational constant", 6.67384e-11, 'm3 / (kg s2)',
                 0.00080e-11, 'CODATA 2010', system='si')

    g0 = Constant('g0', "Standard acceleration of gravity", 9.80665, 'm / s2', 0.0,
                  'CODATA 2010', system='si')

    m_p = Constant('m_p', "Proton mass", 1.672621777e-27, 'kg', 0.000000074e-27,
                   'CODATA 2010', system='si')

    m_n = Constant('m_n', "Neutron mass", 1.674927351e-27, 'kg', 0.000000074e-27,
                   'CODATA 2010', system='si')

    m_e = Constant('m_e', "Electron mass", 9.10938291e-31, 'kg', 0.00000040e-31,
                   'CODATA 2010', system='si')

    u = Constant('u', "Atomic mass", 1.660538921e-27, 'kg', 0.000000073e-27,
                 'CODATA 2010', system='si')

    sigma_sb = Constant('sigma_sb', "Stefan-Boltzmann constant", 5.670373e-8,
                        'W / (K4 m2)', 0.000021e-8, 'CODATA 2010', system='si')

    e = EMConstant('e', 'Electron charge', 1.602176565e-19, 'C', 0.000000035e-19,
                   'CODATA 2010', system='si')

    eps0 = EMConstant('eps0', 'Electric constant', 8.854187817e-12, 'F/m', 0.0,
                      'CODATA 2010', system='si')

    N_A = Constant('N_A', "Avogadro's number", 6.02214129e23, '1 / (mol)',
                   0.00000027e23, 'CODATA 2010', system='si')

    R = Constant('R', "Gas constant", 8.3144621, 'J / (K mol)', 0.0000075,
                 'CODATA 2010', system='si')

    Ryd = Constant('Ryd', 'Rydberg constant', 10973731.568539, '1 / (m)', 0.000055,
                   'CODATA 2010', system='si')

    a0 = Constant('a0', "Bohr radius", 0.52917721092e-10, 'm', 0.00000000017e-10,
                  'CODATA 2010', system='si')

    muB = Constant('muB', "Bohr magneton", 927.400968e-26, 'J/T', 0.00002e-26,
                   'CODATA 2010', system='si')

    alpha = Constant('alpha', "Fine-structure constant", 7.2973525698e-3,
                        '', 0.0000000024e-3, 'CODATA 2010', system='si')

    atm = Constant('atmosphere', "Atmosphere", 101325, 'Pa', 0.0,
                   'CODATA 2010', system='si')

    mu0 = Constant('mu0', "Magnetic constant", 4.0e-7 * np.pi, 'N/A2', 0.0,
                   'CODATA 2010', system='si')

    sigma_T = Constant('sigma_T', "Thomson scattering cross-section",
                       0.6652458734e-28, 'm2', 0.0000000013e-28, 'CODATA 2010',
                       system='si')


for name, value in get_attributes(CODATA2014):
    locals()[name] = value


class ASTRON_V2_0:

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
                      3.0128e28, "W", 0.0, "IAU 2015 Resolution B 2", system='si')

    # Wien wavelength displacement law constant
    b_wien = Constant('b_wien', 'Wien wavelength displacement law constant',
                      2.8977729e-3, 'm K', 00.0000017e-3, 'CODATA 2014',
                      system='si')

    # SOLAR QUANTITIES

    # Solar luminosity
    L_sun = Constant('L_sun', "Nominal solar luminosity", 3.828e26,
                     'W', 0.0, "IAU 2015 Resolution B 3", system='si')

    # Solar mass parameter
    GM_sun = Constant('GM_sun', 'Nominal solar mass parameter', 1.3271244e20,
                      'm3 / (s2)', 0.0, "IAU 2015 Resolution B 3", system='si')

    # Solar mass (derived from mass parameter and gravitational constant)
    #M_sun = Constant('M_sun', "Solar mass", 1.98848e+30,
    #                 'kg', 0.000092e+30,
    #
    M_sun = Constant('M_sun', "Solar mass", GM_sun.value / CODATA2014.G.value,
                     'kg', ((CODATA2014.G.uncertainty / CODATA2014.G.value) *
                            (GM_sun.value / CODATA2014.G.value)),
                     "IAU 2015 Resolution B 3 + CODATA 2014", system='si')

    # Solar radius
    R_sun = Constant('R_sun', "Nominal solar radius", 6.957e8, 'm', 0.0,
                     "IAU 2015 Resolution B 3", system='si')


    # OTHER SOLAR SYSTEM QUANTITIES

    # Jupiter mass parameter
    GM_jup = Constant('GM_jup', 'Nominal Jupiter mass parameter', 1.2668653e17,
                      'm3 / (s2)', 0.0, "IAU 2015 Resolution B 3", system='si')

    # Jupiter mass (derived from mass parameter and gravitational constant)
    #M_jup = Constant('M_earth', "Jupiter mass", 1.898187e+27,
    #                   'kg', 0.000088e+24,
    #                   "IAU 2015 Resolution B 3 + CODATA 2014", system='si')
    M_jup = Constant('M_jup', "Jupiter mass", GM_jup.value / CODATA2014.G.value,
                     'kg', ((CODATA2014.G.uncertainty / CODATA2014.G.value) *
                            (GM_jup.value / CODATA2014.G.value)),
                     "IAU 2015 Resolution B 3 + CODATA 2014", system='si')

    # Jupiter equatorial radius
    R_jup = Constant('R_jup', "Nominal Jupiter equatorial radius", 7.1492e7,
                     'm', 0.0, "IAU 2015 Resolution B 3", system='si')

    # Earth mass parameter
    GM_earth = Constant('GM_earth', 'Nominal Earth mass parameter', 3.986004e14,
                      'm3 / (s2)', 0.0, "IAU 2015 Resolution B 3", system='si')

    # Earth mass (derived from mass parameter and gravitational constant)
    #M_earth = Constant('M_earth', "Earth mass", 5.97236e+24,
    #                   'kg', 0.00028e+24,
    #                   "IAU 2015 Resolution B 3 + CODATA 2014", system='si')
    M_earth = Constant('M_earth', "Earth mass",
                       GM_earth.value / CODATA2014.G.value,
                     'kg', ((CODATA2014.G.uncertainty / CODATA2014.G.value) *
                            (GM_earth.value / CODATA2014.G.value)),
                     "IAU 2015 Resolution B 3 + CODATA 2014", system='si')

    # Earth equatorial radius
    R_earth = Constant('R_earth', "Nominal Earth equatorial radius", 6.3568e6,
                       'm', 0.0, "IAU 2015 Resolution B 3", system='si')

class ASTRON_V1_3:

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
                      3.0128e28, "W", 0.0, "IAU 2015 Resolution B 2", system='si')

    # Wien wavelength displacement law constant
    b_wien = Constant('b_wien', 'Wien wavelength displacement law constant',
                      2.8977721e-3, 'm K', 0.0000026e-3, 'CODATA 2010', system='si')

    # SOLAR QUANTITIES

    # Solar luminosity
    L_sun = Constant('L_sun', "Solar luminosity", 3.846e26, 'W', 0.0005e26,
                     "Allen's Astrophysical Quantities 4th Ed.", system='si')

    # Solar mass
    M_sun = Constant('M_sun', "Solar mass", 1.9891e30, 'kg', 0.00005e30,
                     "Allen's Astrophysical Quantities 4th Ed.", system='si')

    # Solar radius
    R_sun = Constant('R_sun', "Solar radius", 6.95508e8, 'm', 0.00026e8,
                     "Allen's Astrophysical Quantities 4th Ed.", system='si')


    # OTHER SOLAR SYSTEM QUANTITIES

    # Jupiter mass
    M_jup = Constant('M_jup', "Jupiter mass", 1.8987e27, 'kg', 0.00005e27,
                     "Allen's Astrophysical Quantities 4th Ed.", system='si')

    # Jupiter equatorial radius
    R_jup = Constant('R_jup', "Jupiter equatorial radius", 7.1492e7, 'm',
                     0.00005e7, "Allen's Astrophysical Quantities 4th Ed.",
                     system='si')

    # Earth mass
    M_earth = Constant('M_earth', "Earth mass", 5.9742e24, 'kg', 0.00005e24,
                       "Allen's Astrophysical Quantities 4th Ed.", system='si')

    # Earth equatorial radius
    R_earth = Constant('R_earth', "Earth equatorial radius", 6.378136e6, 'm',
                       0.0000005e6, "Allen's Astrophysical Quantities 4th Ed.",
                       system='si')


for name, value in get_attributes(ASTRON_V2_0):
    locals()[name] = value

class ASTROPYCONST20(CODATA2014, ASTRON_V2_0):
    pass

class ASTROPYCONST13(CODATA2010, ASTRON_V1_3):
    pass

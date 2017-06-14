# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from .constant import Constant


# PHYSICAL CONSTANTS

class Codata2014(Constant):
    default_reference = 'CODATA 2014'

    def __new__(cls, abbrev, name, value, unit, uncertainty,
                reference=default_reference, system=None):
        super().__new__(cls, abbrev, name, value, unit, uncertainty,
                        reference, system)


class EMCodata2014(Codata2014):
    """An electromagnetic constant."""

    @property
    def cgs(self):
        """Overridden to raise a `TypeError`
        emphasizing that there are multiple EM extensions to CGS.
        """

        raise TypeError("Cannot convert EM constants to cgs because there "
                        "are different systems for E.M constants within the "
                        "c.g.s system (ESU, Gaussian, etc.). Instead, "
                        "directly use the constant with the appropriate "
                        "suffix (e.g. e.esu, e.gauss, etc.).")


h = Codata2014('h', "Planck constant", 6.626070040e-34,
                        'J s', 0.000000081e-34, system='si')

hbar = Codata2014('hbar', "Reduced Planck constant", 1.054571800e-34,
                    'J s', 0.000000013e-34, system='si')

k_B = Codata2014('k_B', "Boltzmann constant", 1.38064852e-23,
                    'J / (K)', 0.00000079e-23, system='si')

c = Codata2014('c', "Speed of light in vacuum", 299792458.,
                    'm / (s)', 0.0, system='si')


G= Codata2014('G', "Gravitational constant", 6.67408e-11,
                    'm3 / (kg s2)', 0.00031e-11, system='si')

g0 = Codata2014('g0', "Standard acceleration of gravity", 9.80665,
                         'm / s2', 0.0, system='si')

m_p = Codata2014('m_p', "Proton mass", 1.672621898e-27,
               'kg', 0.000000021e-27, system='si')

m_n = Codata2014('m_n', "Neutron mass", 1.674927471e-27,
               'kg', 0.000000021e-27, system='si')

m_e = Codata2014('m_e', "Electron mass", 9.10938356e-31,
               'kg', 0.00000011e-31, system='si')

u = Codata2014('u', "Atomic mass", 1.660539040e-27,
             'kg', 0.000000020e-27, system='si')

sigma_sb = Codata2014('sigma_sb', "Stefan-Boltzmann constant", 5.670367e-8,
                'W / (K4 m2)', 0.000013e-8, system='si')

e = EMCodata2014('e', 'Electron charge', 1.6021766208e-19,
               'C', 0.0000000098e-19, system='si')

eps0 = EMCodata2014('eps0', 'Electric constant', 8.854187817e-12,
                  'F/m', 0.0, system='si')

N_A = Codata2014('N_A', "Avogadro's number", 6.022140857e23,
               '1 / (mol)', 0.000000074e23, system='si')

R = Codata2014('R', "Gas constant", 8.3144598,
             'J / (K mol)', 0.0000048, system='si')

Ryd = Codata2014('Ryd', 'Rydberg constant', 10973731.568508,
               '1 / (m)', 0.000065, system='si')

a0 = Codata2014('a0', "Bohr radius", 0.52917721067e-10,
              'm', 0.00000000012e-10, system='si')

muB = Codata2014('muB', "Bohr magneton", 927.4009994e-26,
               'J/T', 0.00002e-26, system='si')

alpha = Codata2014('alpha', "Fine-structure constant", 7.2973525664e-3,
                 '', 0.0000000017e-3, system='si')

atm = Codata2014('atmosphere', "Atmosphere", 101325,
               'Pa', 0.0, system='si')

mu0 = Codata2014('mu0', "Magnetic constant", 4.0e-7 * np.pi, 'N/A2', 0.0,
               system='si')

sigma_T = Codata2014('sigma_T', "Thomson scattering cross-section",
                   0.66524587158e-28, 'm2', 0.00000000091e-28,
                   system='si')

b_wien = Constant('b_wien', 'Wien wavelength displacement law constant',
                  2.8977729e-3, 'm K', 00.0000017e-3, 'CODATA 2014',
                  system='si')

# cgs constants
# Only constants that cannot be converted directly from S.I. are defined here.

e_esu = EMCodata2014(e.abbrev, e.name, e.value * c.value * 10.0,
                     'statC', e.uncertainty * c.value * 10.0, system='esu')

e_emu = EMCodata2014(e.abbrev, e.name, e.value / 10, 'abC',
                     e.uncertainty / 10, system='emu')

e_gauss = EMCodata2014(e.abbrev, e.name, e.value * c.value * 10.0,
                     'Fr', e.uncertainty * c.value * 10.0, system='gauss')

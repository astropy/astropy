# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
import math

from .constant import Constant, EMConstant


# PHYSICAL CONSTANTS

class CODATA2017(Constant):
    default_reference = 'CODATA 2017'
    _registry = {}
    _has_incompatible_units = set()


class EMCODATA2017(CODATA2017, EMConstant):
    _registry = CODATA2017._registry


# Table 3 from Newell, D. B., et al. 2018, Metrologia, 55, L13
# http://iopscience.iop.org/article/10.1088/1681-7575/aa950a/meta

h = CODATA2017('h', "Planck constant", 6.62607015e-34,
               'J s', 0.0, system='si')

e = EMCODATA2017('e', 'Electron charge', 1.602176634e-19,
                 'C', 0.0, system='si')

k_B = CODATA2017('k_B', "Boltzmann constant", 1.380649e-23,
                 'J / (K)', 0.0, system='si')

N_A = CODATA2017('N_A', "Avogadro's number", 6.02214076e23,
                 '1 / (mol)', 0.0, system='si')

# https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units

c = CODATA2017('c', "Speed of light in vacuum", 299792458.,
               'm / (s)', 0.0, system='si')

# Derived from above.

hbar = CODATA2017('hbar', "Reduced Planck constant", h.value / (2 * math.pi),
                  'J s', 0.0, system='si')

sigma_sb = CODATA2017(
    'sigma_sb', "Stefan-Boltzmann constant",
    2 * math.pi ** 5 * k_B.value ** 4 / (15 * h.value ** 3 * c.value ** 2),
    'W / (K4 m2)', 0.0, system='si')

R = CODATA2017('R', "Gas constant", k_B.value * N_A.value,
               'J / (K mol)', 0.0, system='si')

# Same as CODATA2014.

G = CODATA2017('G', "Gravitational constant", 6.67408e-11,
               'm3 / (kg s2)', 0.00031e-11, system='si')

g0 = CODATA2017('g0', "Standard acceleration of gravity", 9.80665,
                'm / s2', 0.0, system='si')

m_p = CODATA2017('m_p', "Proton mass", 1.672621898e-27,
                 'kg', 0.000000021e-27, system='si')

m_n = CODATA2017('m_n', "Neutron mass", 1.674927471e-27,
                 'kg', 0.000000021e-27, system='si')

m_e = CODATA2017('m_e', "Electron mass", 9.10938356e-31,
                 'kg', 0.00000011e-31, system='si')

u = CODATA2017('u', "Atomic mass", 1.660539040e-27,
               'kg', 0.000000020e-27, system='si')

eps0 = EMCODATA2017('eps0', 'Electric constant', 8.854187817e-12,
                    'F/m', 0.0, system='si')

Ryd = CODATA2017('Ryd', 'Rydberg constant', 10973731.568508,
                 '1 / (m)', 0.000065, system='si')

a0 = CODATA2017('a0', "Bohr radius", 0.52917721067e-10,
                'm', 0.00000000012e-10, system='si')

muB = CODATA2017('muB', "Bohr magneton", 927.4009994e-26,
                 'J/T', 0.00002e-26, system='si')

alpha = CODATA2017('alpha', "Fine-structure constant", 7.2973525664e-3,
                   '', 0.0000000017e-3, system='si')

atm = CODATA2017('atm', "Standard atmosphere", 101325,
                 'Pa', 0.0, system='si')

# Except this one -- from Wikipedia link above.
mu0 = CODATA2017('mu0', "Magnetic constant",
                 2 * h.value * alpha.value / (c.value * e.value ** 2),
                 'N/A2', alpha.uncertainty, system='si')

sigma_T = CODATA2017('sigma_T', "Thomson scattering cross-section",
                     0.66524587158e-28, 'm2', 0.00000000091e-28,
                     system='si')

b_wien = CODATA2017('b_wien', 'Wien wavelength displacement law constant',
                    2.8977729e-3, 'm K', 00.0000017e-3, system='si')

# CGS constants.
# Only constants that cannot be converted directly from S.I. are defined here.
# Because both e and c are exact, these are also exact by definition.

e_esu = EMCODATA2017(e.abbrev, e.name, e.value * c.value * 10.0,
                     'statC', 0.0, system='esu')

e_emu = EMCODATA2017(e.abbrev, e.name, e.value / 10, 'abC',
                     0.0, system='emu')

e_gauss = EMCODATA2017(e.abbrev, e.name, e.value * c.value * 10.0,
                       'Fr', 0.0, system='gauss')

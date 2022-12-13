# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

import numpy as np

from .constant import Constant, EMConstant

# PHYSICAL CONSTANTS


class CODATA2014(Constant):
    default_reference = "CODATA 2014"
    _registry = {}
    _has_incompatible_units = set()


class EMCODATA2014(CODATA2014, EMConstant):
    _registry = CODATA2014._registry


h = CODATA2014(
    "h", "Planck constant", 6.626070040e-34, "J s", 0.000000081e-34, system="si"
)

hbar = CODATA2014(
    "hbar",
    "Reduced Planck constant",
    1.054571800e-34,
    "J s",
    0.000000013e-34,
    system="si",
)

k_B = CODATA2014(
    "k_B", "Boltzmann constant", 1.38064852e-23, "J / (K)", 0.00000079e-23, system="si"
)

c = CODATA2014(
    "c", "Speed of light in vacuum", 299792458.0, "m / (s)", 0.0, system="si"
)


G = CODATA2014(
    "G", "Gravitational constant", 6.67408e-11, "m3 / (kg s2)", 0.00031e-11, system="si"
)

g0 = CODATA2014(
    "g0", "Standard acceleration of gravity", 9.80665, "m / s2", 0.0, system="si"
)

m_p = CODATA2014(
    "m_p", "Proton mass", 1.672621898e-27, "kg", 0.000000021e-27, system="si"
)

m_n = CODATA2014(
    "m_n", "Neutron mass", 1.674927471e-27, "kg", 0.000000021e-27, system="si"
)

m_e = CODATA2014(
    "m_e", "Electron mass", 9.10938356e-31, "kg", 0.00000011e-31, system="si"
)

u = CODATA2014("u", "Atomic mass", 1.660539040e-27, "kg", 0.000000020e-27, system="si")

sigma_sb = CODATA2014(
    "sigma_sb",
    "Stefan-Boltzmann constant",
    5.670367e-8,
    "W / (K4 m2)",
    0.000013e-8,
    system="si",
)

e = EMCODATA2014(
    "e", "Electron charge", 1.6021766208e-19, "C", 0.0000000098e-19, system="si"
)

eps0 = EMCODATA2014(
    "eps0", "Electric constant", 8.854187817e-12, "F/m", 0.0, system="si"
)

N_A = CODATA2014(
    "N_A", "Avogadro's number", 6.022140857e23, "1 / (mol)", 0.000000074e23, system="si"
)

R = CODATA2014("R", "Gas constant", 8.3144598, "J / (K mol)", 0.0000048, system="si")

Ryd = CODATA2014(
    "Ryd", "Rydberg constant", 10973731.568508, "1 / (m)", 0.000065, system="si"
)

a0 = CODATA2014(
    "a0", "Bohr radius", 0.52917721067e-10, "m", 0.00000000012e-10, system="si"
)

muB = CODATA2014(
    "muB", "Bohr magneton", 927.4009994e-26, "J/T", 0.00002e-26, system="si"
)

alpha = CODATA2014(
    "alpha",
    "Fine-structure constant",
    7.2973525664e-3,
    "",
    0.0000000017e-3,
    system="si",
)

atm = CODATA2014("atm", "Standard atmosphere", 101325, "Pa", 0.0, system="si")

mu0 = CODATA2014("mu0", "Magnetic constant", 4.0e-7 * np.pi, "N/A2", 0.0, system="si")

sigma_T = CODATA2014(
    "sigma_T",
    "Thomson scattering cross-section",
    0.66524587158e-28,
    "m2",
    0.00000000091e-28,
    system="si",
)

b_wien = CODATA2014(
    "b_wien",
    "Wien wavelength displacement law constant",
    2.8977729e-3,
    "m K",
    0.0000017e-3,
    system="si",
)

# cgs constants
# Only constants that cannot be converted directly from S.I. are defined here.

e_esu = EMCODATA2014(
    e.abbrev,
    e.name,
    e.value * c.value * 10.0,
    "statC",
    e.uncertainty * c.value * 10.0,
    system="esu",
)

e_emu = EMCODATA2014(
    e.abbrev, e.name, e.value / 10, "abC", e.uncertainty / 10, system="emu"
)

e_gauss = EMCODATA2014(
    e.abbrev,
    e.name,
    e.value * c.value * 10.0,
    "Fr",
    e.uncertainty * c.value * 10.0,
    system="gauss",
)

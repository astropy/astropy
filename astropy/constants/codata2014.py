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


h = CODATA2014("h", "Planck constant", 6.626070040e-34, "J s", 8.1e-42, system="si")

hbar = CODATA2014(
    "hbar",
    "Reduced Planck constant",
    1.054571800e-34,
    "J s",
    1.3e-42,
    system="si",
)

k_B = CODATA2014(
    "k_B", "Boltzmann constant", 1.38064852e-23, "J / (K)", 7.9e-30, system="si"
)

c = CODATA2014(
    "c", "Speed of light in vacuum", 2.99792458e8, "m / (s)", 0.0, system="si"
)


G = CODATA2014(
    "G", "Gravitational constant", 6.67408e-11, "m3 / (kg s2)", 3.1e-15, system="si"
)

g0 = CODATA2014(
    "g0", "Standard acceleration of gravity", 9.80665, "m / s2", 0.0, system="si"
)

m_p = CODATA2014("m_p", "Proton mass", 1.672621898e-27, "kg", 2.1e-35, system="si")

m_n = CODATA2014("m_n", "Neutron mass", 1.674927471e-27, "kg", 2.1e-35, system="si")

m_e = CODATA2014("m_e", "Electron mass", 9.10938356e-31, "kg", 1.1e-38, system="si")

u = CODATA2014("u", "Atomic mass", 1.660539040e-27, "kg", 2.0e-35, system="si")

sigma_sb = CODATA2014(
    "sigma_sb",
    "Stefan-Boltzmann constant",
    5.670367e-8,
    "W / (K4 m2)",
    1.3e-13,
    system="si",
)

e = EMCODATA2014("e", "Electron charge", 1.6021766208e-19, "C", 9.8e-28, system="si")

eps0 = EMCODATA2014(
    "eps0", "Electric constant", 8.854187817e-12, "F/m", 0.0, system="si"
)

N_A = CODATA2014(
    "N_A", "Avogadro's number", 6.022140857e23, "1 / (mol)", 7.4e15, system="si"
)

R = CODATA2014("R", "Gas constant", 8.3144598, "J / (K mol)", 0.0000048, system="si")

Ryd = CODATA2014(
    "Ryd", "Rydberg constant", 1.0973731568508e7, "1 / (m)", 0.000065, system="si"
)

a0 = CODATA2014("a0", "Bohr radius", 5.2917721067e-11, "m", 1.2e-20, system="si")

muB = CODATA2014("muB", "Bohr magneton", 9.274009994e-24, "J/T", 2.0e-31, system="si")

alpha = CODATA2014(
    "alpha",
    "Fine-structure constant",
    7.2973525664e-3,
    "",
    1.7e-12,
    system="si",
)

atm = CODATA2014("atm", "Standard atmosphere", 101325, "Pa", 0.0, system="si")

mu0 = CODATA2014("mu0", "Magnetic constant", 4.0e-7 * np.pi, "N/A2", 0.0, system="si")

sigma_T = CODATA2014(
    "sigma_T",
    "Thomson scattering cross-section",
    6.6524587158e-29,
    "m2",
    9.1e-38,
    system="si",
)

b_wien = CODATA2014(
    "b_wien",
    "Wien wavelength displacement law constant",
    2.8977729e-3,
    "m K",
    1.7e-9,
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

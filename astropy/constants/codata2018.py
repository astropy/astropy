# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

import math

from .constant import Constant, EMConstant

# PHYSICAL CONSTANTS
# https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units


class CODATA2018(Constant):
    default_reference = "CODATA 2018"
    _registry = {}
    _has_incompatible_units = set()


class EMCODATA2018(CODATA2018, EMConstant):
    _registry = CODATA2018._registry


h = CODATA2018("h", "Planck constant", 6.62607015e-34, "J s", 0.0, system="si")

hbar = CODATA2018(
    "hbar", "Reduced Planck constant", h.value / (2 * math.pi), "J s", 0.0, system="si"
)

k_B = CODATA2018("k_B", "Boltzmann constant", 1.380649e-23, "J / (K)", 0.0, system="si")

c = CODATA2018(
    "c", "Speed of light in vacuum", 2.99792458e8, "m / (s)", 0.0, system="si"
)


G = CODATA2018(
    "G", "Gravitational constant", 6.67430e-11, "m3 / (kg s2)", 1.5e-15, system="si"
)

g0 = CODATA2018(
    "g0", "Standard acceleration of gravity", 9.80665, "m / s2", 0.0, system="si"
)

m_p = CODATA2018("m_p", "Proton mass", 1.67262192369e-27, "kg", 5.1e-37, system="si")

m_n = CODATA2018("m_n", "Neutron mass", 1.67492749804e-27, "kg", 9.5e-37, system="si")

m_e = CODATA2018("m_e", "Electron mass", 9.1093837015e-31, "kg", 2.8e-40, system="si")

u = CODATA2018("u", "Atomic mass", 1.66053906660e-27, "kg", 5.0e-37, system="si")

sigma_sb = CODATA2018(
    "sigma_sb",
    "Stefan-Boltzmann constant",
    2 * math.pi**5 * k_B.value**4 / (15 * h.value**3 * c.value**2),
    "W / (K4 m2)",
    0.0,
    system="si",
)

e = EMCODATA2018("e", "Electron charge", 1.602176634e-19, "C", 0.0, system="si")

eps0 = EMCODATA2018(
    "eps0",
    "Vacuum electric permittivity",
    8.8541878128e-12,
    "F/m",
    1.3e-21,
    system="si",
)

N_A = CODATA2018(
    "N_A", "Avogadro's number", 6.02214076e23, "1 / (mol)", 0.0, system="si"
)

R = CODATA2018(
    "R", "Gas constant", k_B.value * N_A.value, "J / (K mol)", 0.0, system="si"
)

Ryd = CODATA2018(
    "Ryd", "Rydberg constant", 1.097373156816e7, "1 / (m)", 0.000021, system="si"
)

a0 = CODATA2018("a0", "Bohr radius", 5.29177210903e-11, "m", 8.0e-21, system="si")

muB = CODATA2018("muB", "Bohr magneton", 9.2740100783e-24, "J/T", 2.8e-33, system="si")

alpha = CODATA2018(
    "alpha",
    "Fine-structure constant",
    7.2973525693e-3,
    "",
    1.1e-12,
    system="si",
)

atm = CODATA2018("atm", "Standard atmosphere", 101325, "Pa", 0.0, system="si")

mu0 = CODATA2018(
    "mu0",
    "Vacuum magnetic permeability",
    1.25663706212e-6,
    "N/A2",
    1.9e-16,
    system="si",
)

sigma_T = CODATA2018(
    "sigma_T",
    "Thomson scattering cross-section",
    6.6524587321e-29,
    "m2",
    6.0e-38,
    system="si",
)

# Formula taken from NIST wall chart.
# The numerical factor is from a numerical solution to the equation for the
# maximum. See https://en.wikipedia.org/wiki/Wien%27s_displacement_law
b_wien = CODATA2018(
    "b_wien",
    "Wien wavelength displacement law constant",
    h.value * c.value / (k_B.value * 4.965114231744276),
    "m K",
    0.0,
    system="si",
)

# CGS constants.
# Only constants that cannot be converted directly from S.I. are defined here.
# Because both e and c are exact, these are also exact by definition.

e_esu = EMCODATA2018(
    e.abbrev, e.name, e.value * c.value * 10.0, "statC", 0.0, system="esu"
)

e_emu = EMCODATA2018(e.abbrev, e.name, e.value / 10, "abC", 0.0, system="emu")

e_gauss = EMCODATA2018(
    e.abbrev, e.name, e.value * c.value * 10.0, "Fr", 0.0, system="gauss"
)

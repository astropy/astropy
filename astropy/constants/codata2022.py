# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

import math

from .constant import Constant, EMConstant

# PHYSICAL CONSTANTS
# https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units


class CODATA2022(Constant):
    default_reference = "CODATA 2022"
    _registry = {}
    _has_incompatible_units = set()


class EMCODATA2022(CODATA2022, EMConstant):
    _registry = CODATA2022._registry


h = CODATA2022("h", "Planck constant", 6.62607015e-34, "J s", 0.0, system="si")

hbar = CODATA2022(
    "hbar", "Reduced Planck constant", h.value / (2 * math.pi), "J s", 0.0, system="si"
)

k_B = CODATA2022("k_B", "Boltzmann constant", 1.380649e-23, "J / (K)", 0.0, system="si")

c = CODATA2022(
    "c", "Speed of light in vacuum", 299792458.0, "m / (s)", 0.0, system="si"
)

G = CODATA2022(
    "G", "Gravitational constant", 6.67430e-11, "m3 / (kg s2)", 0.00015e-11, system="si"
)

g0 = CODATA2022(
    "g0", "Standard acceleration of gravity", 9.80665, "m / s2", 0.0, system="si"
)

m_p = CODATA2022(
    "m_p", "Proton mass", 1.67262192595e-27, "kg", 0.00000000052e-27, system="si"
)

m_n = CODATA2022(
    "m_n", "Neutron mass", 1.67492750056e-27, "kg", 0.00000000085e-27, system="si"
)

m_e = CODATA2022(
    "m_e", "Electron mass", 9.1093837139e-31, "kg", 0.0000000028e-31, system="si"
)

u = CODATA2022(
    "u", "Atomic mass", 1.66053906892e-27, "kg", 0.00000000052e-27, system="si"
)

sigma_sb = CODATA2022(
    "sigma_sb",
    "Stefan-Boltzmann constant",
    2 * math.pi**5 * k_B.value**4 / (15 * h.value**3 * c.value**2),
    "W / (K4 m2)",
    0.0,
    system="si",
)

e = EMCODATA2022("e", "Electron charge", 1.602176634e-19, "C", 0.0, system="si")

eps0 = EMCODATA2022(
    "eps0",
    "Vacuum electric permittivity",
    8.8541878188e-12,
    "F/m",
    0.0000000014e-12,
    system="si",
)

N_A = CODATA2022(
    "N_A", "Avogadro's number", 6.02214076e23, "1 / (mol)", 0.0, system="si"
)

R = CODATA2022(
    "R", "Gas constant", k_B.value * N_A.value, "J / (K mol)", 0.0, system="si"
)

Ryd = CODATA2022(
    "Ryd", "Rydberg constant", 10973731.568157, "1 / (m)", 0.000012, system="si"
)

a0 = CODATA2022(
    "a0", "Bohr radius", 5.29177210544e-11, "m", 0.00000000082e-11, system="si"
)

muB = CODATA2022(
    "muB", "Bohr magneton", 9.2740100657e-24, "J/T", 0.0000000029e-24, system="si"
)

alpha = CODATA2022(
    "alpha",
    "Fine-structure constant",
    7.2973525643e-3,
    "",
    0.0000000011e-3,
    system="si",
)

atm = CODATA2022("atm", "Standard atmosphere", 101325, "Pa", 0.0, system="si")

mu0 = CODATA2022(
    "mu0",
    "Vacuum magnetic permeability",
    1.25663706127e-6,
    "N/A2",
    0.00000000020e-6,
    system="si",
)

sigma_T = CODATA2022(
    "sigma_T",
    "Thomson scattering cross-section",
    6.6524587051e-29,
    "m2",
    0.0000000062e-29,
    system="si",
)

# Formula taken from NIST wall chart.
# The numerical factor is from a numerical solution to the equation for the
# maximum. See https://en.wikipedia.org/wiki/Wien%27s_displacement_law
b_wien = CODATA2022(
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

e_esu = EMCODATA2022(
    e.abbrev, e.name, e.value * c.value * 10.0, "statC", 0.0, system="esu"
)

e_emu = EMCODATA2022(e.abbrev, e.name, e.value / 10, "abC", 0.0, system="emu")

e_gauss = EMCODATA2022(
    e.abbrev, e.name, e.value * c.value * 10.0, "Fr", 0.0, system="gauss"
)

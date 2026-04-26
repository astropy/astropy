# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

import numpy as np

from .constant import Constant, EMConstant

# PHYSICAL CONSTANTS


class CODATA2010(Constant):
    default_reference = "CODATA 2010"
    _registry = {}
    _has_incompatible_units = set()

    def __new__(
        cls,
        abbrev,
        name,
        value,
        unit,
        uncertainty,
        reference=default_reference,
        system=None,
    ):
        return super().__new__(
            cls, abbrev, name, value, unit, uncertainty, reference, system
        )


class EMCODATA2010(CODATA2010, EMConstant):
    _registry = CODATA2010._registry


h = CODATA2010("h", "Planck constant", 6.62606957e-34, "J s", 2.9e-41, system="si")

hbar = CODATA2010(
    "hbar",
    "Reduced Planck constant",
    h.value * 0.5 / np.pi,
    "J s",
    h.uncertainty * 0.5 / np.pi,
    h.reference,
    system="si",
)

k_B = CODATA2010(
    "k_B", "Boltzmann constant", 1.3806488e-23, "J / (K)", 1.3e-29, system="si"
)

c = CODATA2010(
    "c", "Speed of light in vacuum", 2.99792458e8, "m / (s)", 0.0, system="si"
)

G = CODATA2010(
    "G", "Gravitational constant", 6.67384e-11, "m3 / (kg s2)", 8.0e-15, system="si"
)

g0 = CODATA2010(
    "g0", "Standard acceleration of gravity", 9.80665, "m / s2", 0.0, system="si"
)

m_p = CODATA2010("m_p", "Proton mass", 1.672621777e-27, "kg", 7.4e-35, system="si")

m_n = CODATA2010("m_n", "Neutron mass", 1.674927351e-27, "kg", 7.4e-35, system="si")

m_e = CODATA2010("m_e", "Electron mass", 9.10938291e-31, "kg", 4.0e-38, system="si")

u = CODATA2010("u", "Atomic mass", 1.660538921e-27, "kg", 7.3e-35, system="si")

sigma_sb = CODATA2010(
    "sigma_sb",
    "Stefan-Boltzmann constant",
    5.670373e-8,
    "W / (K4 m2)",
    2.1e-13,
    system="si",
)

e = EMCODATA2010("e", "Electron charge", 1.602176565e-19, "C", 3.5e-27, system="si")

eps0 = EMCODATA2010(
    "eps0", "Electric constant", 8.854187817e-12, "F/m", 0.0, system="si"
)

N_A = CODATA2010(
    "N_A", "Avogadro's number", 6.02214129e23, "1 / (mol)", 2.7e16, system="si"
)

R = CODATA2010("R", "Gas constant", 8.3144621, "J / (K mol)", 0.0000075, system="si")

Ryd = CODATA2010(
    "Ryd", "Rydberg constant", 1.0973731568539e7, "1 / (m)", 0.000055, system="si"
)

a0 = CODATA2010("a0", "Bohr radius", 5.2917721092e-11, "m", 1.7e-20, system="si")

muB = CODATA2010("muB", "Bohr magneton", 9.27400968e-24, "J/T", 2.0e-31, system="si")

alpha = CODATA2010(
    "alpha",
    "Fine-structure constant",
    7.2973525698e-3,
    "",
    2.4e-12,
    system="si",
)

atm = CODATA2010("atm", "Standard atmosphere", 101325, "Pa", 0.0, system="si")

mu0 = CODATA2010("mu0", "Magnetic constant", 4.0e-7 * np.pi, "N/A2", 0.0, system="si")

sigma_T = CODATA2010(
    "sigma_T",
    "Thomson scattering cross-section",
    6.652458734e-29,
    "m2",
    1.3e-37,
    system="si",
)

b_wien = Constant(
    "b_wien",
    "Wien wavelength displacement law constant",
    2.8977721e-3,
    "m K",
    2.6e-9,
    "CODATA 2010",
    system="si",
)

# cgs constants
# Only constants that cannot be converted directly from S.I. are defined here.

e_esu = EMCODATA2010(
    e.abbrev,
    e.name,
    e.value * c.value * 10.0,
    "statC",
    e.uncertainty * c.value * 10.0,
    system="esu",
)

e_emu = EMCODATA2010(
    e.abbrev, e.name, e.value / 10, "abC", e.uncertainty / 10, system="emu"
)

e_gauss = EMCODATA2010(
    e.abbrev,
    e.name,
    e.value * c.value * 10.0,
    "Fr",
    e.uncertainty * c.value * 10.0,
    system="gauss",
)

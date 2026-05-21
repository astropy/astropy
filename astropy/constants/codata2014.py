# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

import math
from pathlib import Path

from .codata_parser import parse_codata_file
from .constant import Constant, EMConstant

# PHYSICAL CONSTANTS


class CODATA2014(Constant):
    default_reference = "CODATA 2014"
    _registry = {}
    _has_incompatible_units = set()


class EMCODATA2014(CODATA2014, EMConstant):
    _registry = CODATA2014._registry


CODATA2014txt = parse_codata_file(Path(__file__).with_name("data") / "codata2014.txt")


def get_codata2014_constant(name):
    return CODATA2014txt[name]


h_entry = get_codata2014_constant("Planck constant")
h = CODATA2014(
    "h", h_entry.name, h_entry.value, h_entry.unit, h_entry.uncertainty, system="si"
)

hbar = CODATA2014(
    "hbar", "Reduced Planck constant", h.value / (2 * math.pi), "J s", 0.0, system="si"
)

k_B_entry = get_codata2014_constant("Boltzmann constant")
k_B = CODATA2014(
    "k_B",
    k_B_entry.name,
    k_B_entry.value,
    k_B_entry.unit,
    k_B_entry.uncertainty,
    system="si",
)

c_entry = get_codata2014_constant("speed of light in vacuum")
c = CODATA2014(
    "c", c_entry.name, c_entry.value, c_entry.unit, c_entry.uncertainty, system="si"
)


G_entry = get_codata2014_constant("Newtonian constant of gravitation")
G = CODATA2014(
    "G", G_entry.name, G_entry.value, G_entry.unit, G_entry.uncertainty, system="si"
)

g0_entry = get_codata2014_constant("standard acceleration of gravity")
g0 = CODATA2014(
    "g0",
    g0_entry.name,
    g0_entry.value,
    g0_entry.unit,
    g0_entry.uncertainty,
    system="si",
)

m_p_entry = get_codata2014_constant("proton mass")
m_p = CODATA2014(
    "m_p",
    m_p_entry.name,
    m_p_entry.value,
    m_p_entry.unit,
    m_p_entry.uncertainty,
    system="si",
)

m_n_entry = get_codata2014_constant("neutron mass")
m_n = CODATA2014(
    "m_n",
    m_n_entry.name,
    m_n_entry.value,
    m_n_entry.unit,
    m_n_entry.uncertainty,
    system="si",
)

m_e_entry = get_codata2014_constant("electron mass")
m_e = CODATA2014(
    "m_e",
    m_e_entry.name,
    m_e_entry.value,
    m_e_entry.unit,
    m_e_entry.uncertainty,
    system="si",
)

u_entry = get_codata2014_constant("atomic mass constant")
u = CODATA2014(
    "u", u_entry.name, u_entry.value, u_entry.unit, u_entry.uncertainty, system="si"
)

sigma_sb_entry = get_codata2014_constant("Stefan-Boltzmann constant")
sigma_sb = CODATA2014(
    "sigma_sb",
    sigma_sb_entry.name,
    sigma_sb_entry.value,
    sigma_sb_entry.unit,
    sigma_sb_entry.uncertainty,
    system="si",
)

e_entry = get_codata2014_constant("elementary charge")
e = EMCODATA2014(
    "e", e_entry.name, e_entry.value, e_entry.unit, e_entry.uncertainty, system="si"
)

eps0_entry = get_codata2014_constant("electric constant")
eps0 = EMCODATA2014(
    "eps0",
    eps0_entry.name,
    eps0_entry.value,
    eps0_entry.unit,
    eps0_entry.uncertainty,
    system="si",
)

N_A_entry = get_codata2014_constant("Avogadro constant")
N_A = CODATA2014(
    "N_A",
    N_A_entry.name,
    N_A_entry.value,
    N_A_entry.unit,
    N_A_entry.uncertainty,
    system="si",
)

R_entry = get_codata2014_constant("molar gas constant")
R = CODATA2014(
    "R", R_entry.name, R_entry.value, R_entry.unit, R_entry.uncertainty, system="si"
)

Ryd_entry = get_codata2014_constant("Rydberg constant")
Ryd = CODATA2014(
    "Ryd",
    Ryd_entry.name,
    Ryd_entry.value,
    Ryd_entry.unit,
    Ryd_entry.uncertainty,
    system="si",
)

a0_entry = get_codata2014_constant("Bohr radius")
a0 = CODATA2014(
    "a0",
    a0_entry.name,
    a0_entry.value,
    a0_entry.unit,
    a0_entry.uncertainty,
    system="si",
)

muB_entry = get_codata2014_constant("Bohr magneton")
muB = CODATA2014(
    "muB",
    muB_entry.name,
    muB_entry.value,
    muB_entry.unit,
    muB_entry.uncertainty,
    system="si",
)

alpha_entry = get_codata2014_constant("fine-structure constant")
alpha = CODATA2014(
    "alpha",
    alpha_entry.name,
    alpha_entry.value,
    alpha_entry.unit,
    alpha_entry.uncertainty,
    system="si",
)

atm_entry = get_codata2014_constant("standard atmosphere")
atm = CODATA2014(
    "atm",
    atm_entry.name,
    atm_entry.value,
    atm_entry.unit,
    atm_entry.uncertainty,
    system="si",
)

mu0_entry = get_codata2014_constant("mag. constant")
mu0 = CODATA2014(
    "mu0",
    mu0_entry.name,
    mu0_entry.value,
    mu0_entry.unit,
    mu0_entry.uncertainty,
    system="si",
)

sigma_T_entry = get_codata2014_constant("Thomson cross section")
sigma_T = CODATA2014(
    "sigma_T",
    sigma_T_entry.name,
    sigma_T_entry.value,
    sigma_T_entry.unit,
    sigma_T_entry.uncertainty,
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

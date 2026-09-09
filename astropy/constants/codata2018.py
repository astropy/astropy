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
# https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units


class CODATA2018(Constant):
    default_reference = "CODATA 2018"
    _registry = {}
    _has_incompatible_units = set()


class EMCODATA2018(CODATA2018, EMConstant):
    _registry = CODATA2018._registry


CODATA2018txt = parse_codata_file(Path(__file__).with_name("data") / "codata2018.txt")


def get_codata2018_constant(name):
    return CODATA2018txt[name]


h_entry = get_codata2018_constant("Planck constant")
h = CODATA2018(
    "h", h_entry.name, h_entry.value, h_entry.unit, h_entry.uncertainty, system="si"
)

hbar = CODATA2018(
    "hbar", "Reduced Planck constant", h.value / (2 * math.pi), "J s", 0.0, system="si"
)

k_B_entry = get_codata2018_constant("Boltzmann constant")
k_B = CODATA2018(
    "k_B",
    k_B_entry.name,
    k_B_entry.value,
    k_B_entry.unit,
    k_B_entry.uncertainty,
    system="si",
)

c_entry = get_codata2018_constant("speed of light in vacuum")
c = CODATA2018(
    "c", c_entry.name, c_entry.value, c_entry.unit, c_entry.uncertainty, system="si"
)


G_entry = get_codata2018_constant("Newtonian constant of gravitation")
G = CODATA2018(
    "G", G_entry.name, G_entry.value, G_entry.unit, G_entry.uncertainty, system="si"
)

g0_entry = get_codata2018_constant("standard acceleration of gravity")
g0 = CODATA2018(
    "g0",
    g0_entry.name,
    g0_entry.value,
    g0_entry.unit,
    g0_entry.uncertainty,
    system="si",
)

m_p_entry = get_codata2018_constant("proton mass")
m_p = CODATA2018(
    "m_p",
    m_p_entry.name,
    m_p_entry.value,
    m_p_entry.unit,
    m_p_entry.uncertainty,
    system="si",
)

m_n_entry = get_codata2018_constant("neutron mass")
m_n = CODATA2018(
    "m_n",
    m_n_entry.name,
    m_n_entry.value,
    m_n_entry.unit,
    m_n_entry.uncertainty,
    system="si",
)

m_e_entry = get_codata2018_constant("electron mass")
m_e = CODATA2018(
    "m_e",
    m_e_entry.name,
    m_e_entry.value,
    m_e_entry.unit,
    m_e_entry.uncertainty,
    system="si",
)

u_entry = get_codata2018_constant("atomic mass constant")
u = CODATA2018(
    "u", u_entry.name, u_entry.value, u_entry.unit, u_entry.uncertainty, system="si"
)

sigma_sb = CODATA2018(
    "sigma_sb",
    "Stefan-Boltzmann constant",
    2 * math.pi**5 * k_B.value**4 / (15 * h.value**3 * c.value**2),
    "W / (K4 m2)",
    0.0,
    system="si",
)

e_entry = get_codata2018_constant("elementary charge")
e = EMCODATA2018(
    "e", e_entry.name, e_entry.value, e_entry.unit, e_entry.uncertainty, system="si"
)

eps0_entry = get_codata2018_constant("vacuum electric permittivity")
eps0 = EMCODATA2018(
    "eps0",
    eps0_entry.name,
    eps0_entry.value,
    eps0_entry.unit,
    eps0_entry.uncertainty,
    system="si",
)

N_A_entry = get_codata2018_constant("Avogadro constant")
N_A = CODATA2018(
    "N_A",
    N_A_entry.name,
    N_A_entry.value,
    N_A_entry.unit,
    N_A_entry.uncertainty,
    system="si",
)

R = CODATA2018(
    "R", "Gas constant", k_B.value * N_A.value, "J / (K mol)", 0.0, system="si"
)

Ryd_entry = get_codata2018_constant("Rydberg constant")
Ryd = CODATA2018(
    "Ryd",
    Ryd_entry.name,
    Ryd_entry.value,
    Ryd_entry.unit,
    Ryd_entry.uncertainty,
    system="si",
)

a0_entry = get_codata2018_constant("Bohr radius")
a0 = CODATA2018(
    "a0",
    a0_entry.name,
    a0_entry.value,
    a0_entry.unit,
    a0_entry.uncertainty,
    system="si",
)

muB_entry = get_codata2018_constant("Bohr magneton")
muB = CODATA2018(
    "muB",
    muB_entry.name,
    muB_entry.value,
    muB_entry.unit,
    muB_entry.uncertainty,
    system="si",
)

alpha_entry = get_codata2018_constant("fine-structure constant")
alpha = CODATA2018(
    "alpha",
    alpha_entry.name,
    alpha_entry.value,
    alpha_entry.unit,
    alpha_entry.uncertainty,
    system="si",
)

atm_entry = get_codata2018_constant("standard atmosphere")
atm = CODATA2018(
    "atm",
    atm_entry.name,
    atm_entry.value,
    atm_entry.unit,
    atm_entry.uncertainty,
    system="si",
)

mu0_entry = get_codata2018_constant("vacuum mag. permeability")
mu0 = CODATA2018(
    "mu0",
    mu0_entry.name,
    mu0_entry.value,
    mu0_entry.unit,
    mu0_entry.uncertainty,
    system="si",
)

sigma_T_entry = get_codata2018_constant("Thomson cross section")
sigma_T = CODATA2018(
    "sigma_T",
    sigma_T_entry.name,
    sigma_T_entry.value,
    sigma_T_entry.unit,
    sigma_T_entry.uncertainty,
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

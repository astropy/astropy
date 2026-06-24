# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

from pathlib import Path

import numpy as np

from .codata_parser import parse_codata_file
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


CODATA2010txt = parse_codata_file(Path(__file__).with_name("data") / "codata2010.txt")


def get_codata2010_constant(name):
    return CODATA2010txt[name]


h_entry = get_codata2010_constant("Planck constant")
h = CODATA2010(
    "h", h_entry.name, h_entry.value, h_entry.unit, h_entry.uncertainty, system="si"
)

hbar = CODATA2010(
    "hbar",
    "Reduced Planck constant",
    h.value * 0.5 / np.pi,
    "J s",
    h.uncertainty * 0.5 / np.pi,
    h.reference,
    system="si",
)

k_B_entry = get_codata2010_constant("Boltzmann constant")
k_B = CODATA2010(
    "k_B",
    k_B_entry.name,
    k_B_entry.value,
    k_B_entry.unit,
    k_B_entry.uncertainty,
    system="si",
)

c_entry = get_codata2010_constant("speed of light in vacuum")
c = CODATA2010(
    "c", c_entry.name, c_entry.value, c_entry.unit, c_entry.uncertainty, system="si"
)

G_entry = get_codata2010_constant("Newtonian constant of gravitation")
G = CODATA2010(
    "G", G_entry.name, G_entry.value, G_entry.unit, G_entry.uncertainty, system="si"
)

g0_entry = get_codata2010_constant("standard acceleration of gravity")
g0 = CODATA2010(
    "g0",
    g0_entry.name,
    g0_entry.value,
    g0_entry.unit,
    g0_entry.uncertainty,
    system="si",
)

m_p_entry = get_codata2010_constant("proton mass")
m_p = CODATA2010(
    "m_p",
    m_p_entry.name,
    m_p_entry.value,
    m_p_entry.unit,
    m_p_entry.uncertainty,
    system="si",
)

m_n_entry = get_codata2010_constant("neutron mass")
m_n = CODATA2010(
    "m_n",
    m_n_entry.name,
    m_n_entry.value,
    m_n_entry.unit,
    m_n_entry.uncertainty,
    system="si",
)

m_e_entry = get_codata2010_constant("electron mass")
m_e = CODATA2010(
    "m_e",
    m_e_entry.name,
    m_e_entry.value,
    m_e_entry.unit,
    m_e_entry.uncertainty,
    system="si",
)

u_entry = get_codata2010_constant("atomic mass constant")
u = CODATA2010(
    "u", u_entry.name, u_entry.value, u_entry.unit, u_entry.uncertainty, system="si"
)

sigma_sb_entry = get_codata2010_constant("Stefan-Boltzmann constant")
sigma_sb = CODATA2010(
    "sigma_sb",
    sigma_sb_entry.name,
    sigma_sb_entry.value,
    sigma_sb_entry.unit,
    sigma_sb_entry.uncertainty,
    system="si",
)

e_entry = get_codata2010_constant("elementary charge")
e = EMCODATA2010(
    "e", e_entry.name, e_entry.value, e_entry.unit, e_entry.uncertainty, system="si"
)

eps0_entry = get_codata2010_constant("electric constant")
eps0 = EMCODATA2010(
    "eps0",
    eps0_entry.name,
    eps0_entry.value,
    eps0_entry.unit,
    eps0_entry.uncertainty,
    system="si",
)

N_A_entry = get_codata2010_constant("Avogadro constant")
N_A = CODATA2010(
    "N_A",
    N_A_entry.name,
    N_A_entry.value,
    N_A_entry.unit,
    N_A_entry.uncertainty,
    system="si",
)

R_entry = get_codata2010_constant("molar gas constant")
R = CODATA2010(
    "R", R_entry.name, R_entry.value, R_entry.unit, R_entry.uncertainty, system="si"
)

Ryd_entry = get_codata2010_constant("Rydberg constant")
Ryd = CODATA2010(
    "Ryd",
    Ryd_entry.name,
    Ryd_entry.value,
    Ryd_entry.unit,
    Ryd_entry.uncertainty,
    system="si",
)

a0_entry = get_codata2010_constant("Bohr radius")
a0 = CODATA2010(
    "a0",
    a0_entry.name,
    a0_entry.value,
    a0_entry.unit,
    a0_entry.uncertainty,
    system="si",
)

muB_entry = get_codata2010_constant("Bohr magneton")
muB = CODATA2010(
    "muB",
    muB_entry.name,
    muB_entry.value,
    muB_entry.unit,
    muB_entry.uncertainty,
    system="si",
)

alpha_entry = get_codata2010_constant("fine-structure constant")
alpha = CODATA2010(
    "alpha",
    alpha_entry.name,
    alpha_entry.value,
    alpha_entry.unit,
    alpha_entry.uncertainty,
    system="si",
)

atm_entry = get_codata2010_constant("standard atmosphere")
atm = CODATA2010(
    "atm",
    atm_entry.name,
    atm_entry.value,
    atm_entry.unit,
    atm_entry.uncertainty,
    system="si",
)

mu0_entry = get_codata2010_constant("mag. constant")
mu0 = CODATA2010(
    "mu0",
    mu0_entry.name,
    mu0_entry.value,
    mu0_entry.unit,
    mu0_entry.uncertainty,
    system="si",
)

sigma_T_entry = get_codata2010_constant("Thomson cross section")
sigma_T = CODATA2010(
    "sigma_T",
    sigma_T_entry.name,
    sigma_T_entry.value,
    sigma_T_entry.unit,
    sigma_T_entry.uncertainty,
    system="si",
)

b_wien = Constant(
    "b_wien",
    "Wien wavelength displacement law constant",
    2.8977721e-3,
    "m K",
    0.0000026e-3,
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

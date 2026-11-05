# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pathlib import Path

import numpy as np
import pytest

from astropy.constants import codata2010, codata2014, codata2018
from astropy.constants.codata_parser import parse_codata_file
from astropy.constants.constant import Constant, EMConstant


@pytest.mark.parametrize(
    ("module", "required_symbols"),
    [
        (
            codata2010,
            ("h", "c", "G", "m_e", "eps0", "mu0", "sigma_T", "e", "e_esu"),
        ),
        (
            codata2014,
            ("h", "c", "G", "m_e", "eps0", "mu0", "sigma_T", "e", "e_esu"),
        ),
        (
            codata2018,
            ("h", "c", "G", "m_e", "eps0", "mu0", "sigma_T", "e", "e_esu"),
        ),
    ],
)
# this test ensures that the constants modules expose the expected symbols and that they are instances of Constant
def test_codata_modules_expose_expected_symbols(module, required_symbols):
    for symbol in required_symbols:
        value = getattr(module, symbol)
        assert isinstance(value, Constant)


@pytest.mark.parametrize(
    ("module", "txt_name", "mapping"),
    [
        (
            codata2010,
            "codata2010.txt",
            {
                "h": "Planck constant",
                "c": "speed of light in vacuum",
                "G": "Newtonian constant of gravitation",
                "m_e": "electron mass",
                "eps0": "electric constant",
                "mu0": "mag. constant",
                "sigma_T": "Thomson cross section",
            },
        ),
        (
            codata2014,
            "codata2014.txt",
            {
                "h": "Planck constant",
                "c": "speed of light in vacuum",
                "G": "Newtonian constant of gravitation",
                "m_e": "electron mass",
                "eps0": "electric constant",
                "mu0": "mag. constant",
                "sigma_T": "Thomson cross section",
            },
        ),
        (
            codata2018,
            "codata2018.txt",
            {
                "h": "Planck constant",
                "c": "speed of light in vacuum",
                "G": "Newtonian constant of gravitation",
                "m_e": "electron mass",
                "eps0": "vacuum electric permittivity",
                "mu0": "vacuum mag. permeability",
                "sigma_T": "Thomson cross section",
            },
        ),
    ],
)
# this test ensures that the values in the constants modules match the entries in the corresponding CODATA text files
def test_codata_values_match_txt_entries(module, txt_name, mapping):
    txt_path = Path(module.__file__).with_name("data") / txt_name
    entries = parse_codata_file(txt_path)

    for symbol, entry_name in mapping.items():
        constant = getattr(module, symbol)
        entry = entries[entry_name]
        assert constant.value == entry.value
        assert constant.uncertainty == entry.uncertainty


@pytest.mark.parametrize("module", [codata2010, codata2014, codata2018])
def test_em_cgs_definitions_preserved(module):
    assert isinstance(module.e, EMConstant)
    assert module.e_esu.value == module.e.value * module.c.value * 10.0
    assert module.e_emu.value == module.e.value / 10.0
    assert module.e_gauss.value == module.e.value * module.c.value * 10.0


def test_codata2010_hbar_relationship_preserved():
    assert codata2010.hbar.value == codata2010.h.value * 0.5 / np.pi
    assert codata2010.hbar.uncertainty == codata2010.h.uncertainty * 0.5 / np.pi


def test_codata2018_derived_relationships_preserved():
    expected_sigma_sb = (
        2
        * np.pi**5
        * codata2018.k_B.value**4
        / (15 * codata2018.h.value**3 * codata2018.c.value**2)
    )
    assert codata2018.sigma_sb.value == expected_sigma_sb
    assert codata2018.R.value == codata2018.k_B.value * codata2018.N_A.value

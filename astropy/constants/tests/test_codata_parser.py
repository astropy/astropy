# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pathlib import Path

from astropy.constants.codata_parser import parse_codata_file

path = Path(__file__).resolve().parents[1] / "data" / "codata2022.txt"
file_content = parse_codata_file(path)


def test_parse_codata_file():
    # make sure the file is read and parsed without error and that some expected constants are present
    content = parse_codata_file(path)
    assert content
    assert "Planck constant" in content
    assert "speed of light in vacuum" in content


def test_parse_exact_constant():
    # Planck constant should parse with expected numeric value, uncertainty, and unit
    planck = file_content["Planck constant"]

    assert planck.value == 6.62607015e-34
    assert planck.uncertainty == 0.0
    assert planck.unit == "J Hz^-1"


def test_parse_truncated_exact_constant():
    # Stefan-Boltzmann constant should parse with expected numeric value, uncertainty, and unit
    sigma = file_content["Stefan-Boltzmann constant"]

    assert sigma.value == 5.670374419e-8
    assert sigma.uncertainty == 0.0
    assert sigma.unit == "W m^-2 K^-4"


def test_parse_nonexact_constant():
    # Newtonian constant of gravitation should parse with expected value, uncertainty, and unit
    g = file_content["Newtonian constant of gravitation"]

    assert g.value == 6.67430e-11
    assert g.uncertainty == 0.00015e-11
    assert g.unit == "m^3 kg^-1 s^-2"


def test_parse_constant_without_unit():
    # W to Z mass ratio is dimensionless and has no unit.
    ratio = file_content["W to Z mass ratio"]

    assert ratio.value == 0.88145
    assert ratio.uncertainty == 0.00013
    assert ratio.unit == ""

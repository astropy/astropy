# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import astropy.cosmology.constants as cosmology_constants


def test_c():
    assert cosmology_constants.c.unit == "km / s"
    assert cosmology_constants.c.value == pytest.approx(299792.458)


def test_G():
    assert cosmology_constants.G.unit == "pc km2 / (s2 solMass)"
    assert cosmology_constants.G.value == pytest.approx(4.30091727e-3, rel=1e-9)

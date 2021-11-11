# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for the photometric module.  Note that this is shorter than
might be expected because a lot of the relevant tests that deal
with magnidues are in `test_logarithmic.py`
"""

from astropy.tests.helper import assert_quantity_allclose
from astropy.units import AA, ABflux, Jy, Magnitude, STflux, cm, erg, mgy, nmgy, s, zero_point_flux


def test_maggies():
    assert_quantity_allclose(1e-9*mgy, 1*nmgy)
    assert_quantity_allclose(Magnitude((1*nmgy).to(mgy)).value, 22.5)


def test_maggies_zpts():
    assert_quantity_allclose((1*nmgy).to(ABflux, zero_point_flux(1*ABflux)), 3631e-9*Jy, rtol=1e-3)

    ST_base_unit = erg * cm**-2 / s / AA

    stmgy = (10*mgy).to(STflux, zero_point_flux(1*ST_base_unit))
    assert_quantity_allclose(stmgy, 10*ST_base_unit)

    mgyst = (2*ST_base_unit).to(mgy, zero_point_flux(0.5*ST_base_unit))
    assert_quantity_allclose(mgyst, 4*mgy)

    nmgyst = (5.e-10*ST_base_unit).to(mgy, zero_point_flux(0.5*ST_base_unit))
    assert_quantity_allclose(nmgyst, 1*nmgy)

# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ...tests.helper import assert_quantity_allclose

from .. import Magnitude, mgy, nmgy, ABflux, STflux, phot_zero_point, Jy
from .. import erg, cm, s, AA


def test_maggies():
    assert_quantity_allclose(1e-9*mgy, 1*nmgy)
    assert_quantity_allclose(Magnitude((1*nmgy).to(mgy)).value, 22.5)


def test_maggies_zpts():
    assert_quantity_allclose((1*nmgy).to(ABflux, phot_zero_point()), 3631e-9*Jy)

    ST_base_unit = erg * cm**-2 / s / AA

    stmgy = (10*mgy).to(STflux, phot_zero_point(1*ST_base_unit))
    assert_quantity_allclose(stmgy, 10*ST_base_unit)

    mgyst = (2*ST_base_unit).to(mgy, phot_zero_point(0.5*ST_base_unit))
    assert_quantity_allclose(mgyst, 4*mgy)

    nmgyst = (5.e-10*ST_base_unit).to(mgy, phot_zero_point(0.5*ST_base_unit))
    assert_quantity_allclose(nmgyst, 1*nmgy)

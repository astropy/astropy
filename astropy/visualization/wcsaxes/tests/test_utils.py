# Licensed under a 3-clause BSD style license - see LICENSE.rst


from numpy.testing import assert_almost_equal

from astropy import units as u
from astropy.coordinates import Angle, Galactic, HADec
from astropy.tests.helper import (
    assert_quantity_allclose as assert_almost_equal_quantity,
)
from astropy.visualization.wcsaxes.utils import (
    get_coord_meta,
    select_step_degree,
    select_step_hour,
    select_step_scalar,
)


def test_select_step_degree():
    assert_almost_equal_quantity(select_step_degree(127 * u.deg), 180.0 * u.deg)
    assert_almost_equal_quantity(select_step_degree(44 * u.deg), 45.0 * u.deg)
    assert_almost_equal_quantity(select_step_degree(18 * u.arcmin), 15 * u.arcmin)
    assert_almost_equal_quantity(select_step_degree(3.4 * u.arcmin), 3 * u.arcmin)
    assert_almost_equal_quantity(select_step_degree(2 * u.arcmin), 2 * u.arcmin)
    assert_almost_equal_quantity(select_step_degree(59 * u.arcsec), 1 * u.arcmin)
    assert_almost_equal_quantity(select_step_degree(33 * u.arcsec), 30 * u.arcsec)
    assert_almost_equal_quantity(select_step_degree(2.2 * u.arcsec), 2 * u.arcsec)
    assert_almost_equal_quantity(select_step_degree(0.8 * u.arcsec), 1 * u.arcsec)
    assert_almost_equal_quantity(select_step_degree(0.2 * u.arcsec), 0.2 * u.arcsec)
    assert_almost_equal_quantity(select_step_degree(0.11 * u.arcsec), 0.1 * u.arcsec)
    assert_almost_equal_quantity(select_step_degree(0.022 * u.arcsec), 0.02 * u.arcsec)
    assert_almost_equal_quantity(
        select_step_degree(0.0043 * u.arcsec), 0.005 * u.arcsec
    )
    assert_almost_equal_quantity(
        select_step_degree(0.00083 * u.arcsec), 0.001 * u.arcsec
    )
    assert_almost_equal_quantity(
        select_step_degree(0.000027 * u.arcsec), 0.00002 * u.arcsec
    )


def test_select_step_hour():
    assert_almost_equal_quantity(select_step_hour(127 * u.deg), 8.0 * u.hourangle)
    assert_almost_equal_quantity(select_step_hour(44 * u.deg), 3.0 * u.hourangle)
    assert_almost_equal_quantity(select_step_hour(18 * u.arcmin), 15 * u.arcmin)
    assert_almost_equal_quantity(select_step_hour(3.4 * u.arcmin), 3 * u.arcmin)
    assert_almost_equal_quantity(select_step_hour(2 * u.arcmin), 1.5 * u.arcmin)
    assert_almost_equal_quantity(select_step_hour(59 * u.arcsec), 1 * u.arcmin)
    assert_almost_equal_quantity(select_step_hour(33 * u.arcsec), 30 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(2.2 * u.arcsec), 3.0 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.8 * u.arcsec), 0.75 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.2 * u.arcsec), 0.15 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.11 * u.arcsec), 0.15 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.022 * u.arcsec), 0.03 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.0043 * u.arcsec), 0.003 * u.arcsec)
    assert_almost_equal_quantity(
        select_step_hour(0.00083 * u.arcsec), 0.00075 * u.arcsec
    )
    assert_almost_equal_quantity(
        select_step_hour(0.000027 * u.arcsec), 0.00003 * u.arcsec
    )


def test_select_step_scalar():
    assert_almost_equal(select_step_scalar(33122.0), 50000.0)
    assert_almost_equal(select_step_scalar(433.0), 500.0)
    assert_almost_equal(select_step_scalar(12.3), 10)
    assert_almost_equal(select_step_scalar(3.3), 5.0)
    assert_almost_equal(select_step_scalar(0.66), 0.5)
    assert_almost_equal(select_step_scalar(0.0877), 0.1)
    assert_almost_equal(select_step_scalar(0.00577), 0.005)
    assert_almost_equal(select_step_scalar(0.00022), 0.0002)
    assert_almost_equal(select_step_scalar(0.000012), 0.00001)
    assert_almost_equal(select_step_scalar(0.000000443), 0.0000005)


def test_get_coord_meta():
    galactic_meta = get_coord_meta(Galactic)
    assert galactic_meta["name"] == ["l", "b"]
    assert galactic_meta["wrap"] == (Angle(360 * u.deg), None)
    assert galactic_meta["unit"] == (u.deg, u.deg)

    hadec_meta = get_coord_meta(HADec)
    assert hadec_meta["name"] == ["ha", "dec"]
    assert hadec_meta["wrap"] == (Angle(180 * u.deg), None)
    assert hadec_meta["unit"] == (u.hourangle, u.deg)

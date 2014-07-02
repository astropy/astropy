from numpy.testing import assert_almost_equal

from astropy import units as u
from astropy.wcs import WCS

from ..utils import (select_step_degree, select_step_hour, select_step_scalar,
                     coord_type_from_ctype)


def assert_almost_equal_quantity(q1, q2):
    assert_almost_equal(q1, q2.to(q1.unit))


def test_select_step_degree():
    assert_almost_equal_quantity(select_step_degree(127 * u.deg), 180. * u.deg)
    assert_almost_equal_quantity(select_step_degree(44 * u.deg), 45. * u.deg)
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
    assert_almost_equal_quantity(select_step_degree(0.0043 * u.arcsec), 0.005 * u.arcsec)
    assert_almost_equal_quantity(select_step_degree(0.00083 * u.arcsec), 0.001 * u.arcsec)
    assert_almost_equal_quantity(select_step_degree(0.000027 * u.arcsec), 0.00002 * u.arcsec)


def test_select_step_hour():
    assert_almost_equal_quantity(select_step_hour(127 * u.deg), 8. * u.hourangle)
    assert_almost_equal_quantity(select_step_hour(44 * u.deg), 3. * u.hourangle)
    assert_almost_equal_quantity(select_step_hour(18 * u.arcmin), 15 * u.arcmin)
    assert_almost_equal_quantity(select_step_hour(3.4 * u.arcmin), 3 * u.arcmin)
    assert_almost_equal_quantity(select_step_hour(2 * u.arcmin), 1.5 * u.arcmin)
    assert_almost_equal_quantity(select_step_hour(59 * u.arcsec), 1 * u.arcmin)
    assert_almost_equal_quantity(select_step_hour(33 * u.arcsec), 30 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(2.2 * u.arcsec), 3. * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.8 * u.arcsec), 0.75 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.2 * u.arcsec), 0.15 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.11 * u.arcsec), 0.15 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.022 * u.arcsec), 0.03 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.0043 * u.arcsec), 0.003 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.00083 * u.arcsec), 0.00075 * u.arcsec)
    assert_almost_equal_quantity(select_step_hour(0.000027 * u.arcsec), 0.00003 * u.arcsec)


def test_select_step_scalar():
    assert_almost_equal(select_step_scalar(33122.), 50000.)
    assert_almost_equal(select_step_scalar(433.), 500.)
    assert_almost_equal(select_step_scalar(12.3), 10)
    assert_almost_equal(select_step_scalar(3.3), 5.)
    assert_almost_equal(select_step_scalar(0.66), 0.5)
    assert_almost_equal(select_step_scalar(0.0877), 0.1)
    assert_almost_equal(select_step_scalar(0.00577), 0.005)
    assert_almost_equal(select_step_scalar(0.00022), 0.0002)
    assert_almost_equal(select_step_scalar(0.000012), 0.00001)
    assert_almost_equal(select_step_scalar(0.000000443), 0.0000005)


def test_coord_type_from_ctype():
    assert coord_type_from_ctype(' LON') == ('longitude', None)
    assert coord_type_from_ctype(' LAT') == ('latitude', None)
    assert coord_type_from_ctype('HPLN') == ('longitude', 180.)
    assert coord_type_from_ctype('HPLT') == ('latitude', None)
    assert coord_type_from_ctype('RA--') == ('longitude', None)
    assert coord_type_from_ctype('DEC-') == ('latitude', None)
    assert coord_type_from_ctype('spam') == ('scalar', None)


def test_get_coordinate_frame():

    from ..utils import get_coordinate_frame, register_frame_identifier, reset_frame_identifiers
    from astropy.coordinates import FK5, Galactic
    from astropy.tests.helper import pytest

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    assert get_coordinate_frame(wcs) is FK5

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']

    assert get_coordinate_frame(wcs) is Galactic

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['SOLARX', 'SOLARY']

    with pytest.raises(ValueError) as exc:
        get_coordinate_frame(wcs)
    assert exc.value.args[0] == "Frame not supported: SOLARX/SOLARY"

    class SolarXY(object):
        pass

    def identify_solar(wcs):
        if wcs.wcs.ctype[0] == "SOLARX" and wcs.wcs.ctype[1] == "SOLARY":
            return SolarXY

    register_frame_identifier(identify_solar)

    assert get_coordinate_frame(wcs) is SolarXY

    reset_frame_identifiers()

    with pytest.raises(ValueError) as exc:
        get_coordinate_frame(wcs)
    assert exc.value.args[0] == "Frame not supported: SOLARX/SOLARY"

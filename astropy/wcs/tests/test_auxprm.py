# -*- coding: utf-8 -*-

# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Tests for the auxiliary parameters contained in wcsaux

from numpy.testing import assert_allclose

from astropy.io import fits
from astropy.wcs import WCS


STR_EXPECTED_EMPTY = """
rsun_ref:
dsun_obs:
crln_obs:
hgln_obs:
hglt_obs:""".lstrip()


def test_empty():
    w = WCS(naxis=1)
    assert w.wcs.aux.rsun_ref is None
    assert w.wcs.aux.dsun_obs is None
    assert w.wcs.aux.crln_obs is None
    assert w.wcs.aux.hgln_obs is None
    assert w.wcs.aux.hglt_obs is None
    assert str(w.wcs.aux) == STR_EXPECTED_EMPTY


HEADER_SOLAR = fits.Header.fromstring("""
WCSAXES =                    2 / Number of coordinate axes
CRPIX1  =                 64.5 / Pixel coordinate of reference point
CRPIX2  =                 64.5 / Pixel coordinate of reference point
PC1_1   =     0.99999994260024 / Coordinate transformation matrix element
PC1_2   = -0.00033882076120692 / Coordinate transformation matrix element
PC2_1   =  0.00033882076120692 / Coordinate transformation matrix element
PC2_2   =     0.99999994260024 / Coordinate transformation matrix element
CDELT1  =   0.0053287911111111 / [deg] Coordinate increment at reference point
CDELT2  =   0.0053287911111111 / [deg] Coordinate increment at reference point
CUNIT1  = 'deg'                / Units of coordinate increment and value
CUNIT2  = 'deg'                / Units of coordinate increment and value
CTYPE1  = 'HPLN-TAN'           / Coordinate type codegnomonic projection
CTYPE2  = 'HPLT-TAN'           / Coordinate type codegnomonic projection
CRVAL1  =  -0.0012589367249586 / [deg] Coordinate value at reference point
CRVAL2  =  0.00079599300143911 / [deg] Coordinate value at reference point
LONPOLE =                180.0 / [deg] Native longitude of celestial pole
LATPOLE =  0.00079599300143911 / [deg] Native latitude of celestial pole
DATE-OBS= '2011-02-15T00:00:00.34' / ISO-8601 time of observation
MJD-OBS =      55607.000003935 / [d] MJD at start of observation
RSUN_REF=          696000000.0 / [m] Solar radius
DSUN_OBS=       147724815128.0 / [m] Distance from centre of Sun to observer
CRLN_OBS=            22.814522 / [deg] Carrington heliographic lng of observer
CRLT_OBS=            -6.820544 / [deg] Heliographic latitude of observer
HGLN_OBS=             8.431123 / [deg] Stonyhurst heliographic lng of observer
HGLT_OBS=            -6.820544 / [deg] Heliographic latitude of observer
""".lstrip(), sep='\n')


STR_EXPECTED_GET = """
rsun_ref: 696000000.000000
dsun_obs: 147724815128.000000
crln_obs: 22.814522
hgln_obs: 8.431123
hglt_obs: -6.820544""".lstrip()


def test_solar_aux_get():
    w = WCS(HEADER_SOLAR)
    assert_allclose(w.wcs.aux.rsun_ref, 696000000)
    assert_allclose(w.wcs.aux.dsun_obs, 147724815128)
    assert_allclose(w.wcs.aux.crln_obs, 22.814522)
    assert_allclose(w.wcs.aux.hgln_obs, 8.431123)
    assert_allclose(w.wcs.aux.hglt_obs, -6.820544)
    assert str(w.wcs.aux) == STR_EXPECTED_GET


STR_EXPECTED_SET = """
rsun_ref: 698000000.000000
dsun_obs: 140000000000.000000
crln_obs: 10.000000
hgln_obs: 30.000000
hglt_obs: 40.000000""".lstrip()


def test_solar_aux_set():

    w = WCS(HEADER_SOLAR)

    w.wcs.aux.rsun_ref = 698000000
    assert_allclose(w.wcs.aux.rsun_ref, 698000000)

    w.wcs.aux.dsun_obs = 140000000000
    assert_allclose(w.wcs.aux.dsun_obs, 140000000000)

    w.wcs.aux.crln_obs = 10.
    assert_allclose(w.wcs.aux.crln_obs, 10.)

    w.wcs.aux.hgln_obs = 30.
    assert_allclose(w.wcs.aux.hgln_obs, 30.)

    w.wcs.aux.hglt_obs = 40.
    assert_allclose(w.wcs.aux.hglt_obs, 40.)

    assert str(w.wcs.aux) == STR_EXPECTED_SET

    header = w.to_header()
    assert_allclose(header['RSUN_REF'], 698000000)
    assert_allclose(header['DSUN_OBS'], 140000000000)
    assert_allclose(header['CRLN_OBS'], 10.)
    assert_allclose(header['HGLN_OBS'], 30.)
    assert_allclose(header['HGLT_OBS'], 40.)


def test_set_aux_on_empty():

    w = WCS(naxis=2)

    w.wcs.aux.rsun_ref = 698000000
    assert_allclose(w.wcs.aux.rsun_ref, 698000000)

    w.wcs.aux.dsun_obs = 140000000000
    assert_allclose(w.wcs.aux.dsun_obs, 140000000000)

    w.wcs.aux.crln_obs = 10.
    assert_allclose(w.wcs.aux.crln_obs, 10.)

    w.wcs.aux.hgln_obs = 30.
    assert_allclose(w.wcs.aux.hgln_obs, 30.)

    w.wcs.aux.hglt_obs = 40.
    assert_allclose(w.wcs.aux.hglt_obs, 40.)

    assert str(w.wcs.aux) == STR_EXPECTED_SET

    header = w.to_header()
    assert_allclose(header['RSUN_REF'], 698000000)
    assert_allclose(header['DSUN_OBS'], 140000000000)
    assert_allclose(header['CRLN_OBS'], 10.)
    assert_allclose(header['HGLN_OBS'], 30.)
    assert_allclose(header['HGLT_OBS'], 40.)

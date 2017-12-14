# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for putting velocity differentials into SkyCoord objects.
"""

import pytest

import numpy as np

from ... import units as u
from ...tests.helper import assert_quantity_allclose
from .. import (SkyCoord, ICRS, SphericalRepresentation, SphericalDifferential,
                SphericalCosLatDifferential, Galactic, PrecessedGeocentric)

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def test_creation_frameobjs():
    i = ICRS(1*u.deg, 2*u.deg, pm_ra_cosdec=.2*u.mas/u.yr, pm_dec=.1*u.mas/u.yr)
    sc = SkyCoord(i)

    for attrnm in ['ra', 'dec', 'pm_ra_cosdec', 'pm_dec']:
        assert_quantity_allclose(getattr(i, attrnm), getattr(sc, attrnm))

    sc_nod = SkyCoord(ICRS(1*u.deg, 2*u.deg))

    for attrnm in ['ra', 'dec']:
        assert_quantity_allclose(getattr(sc, attrnm), getattr(sc_nod, attrnm))


def test_creation_attrs():
    sc1 = SkyCoord(1*u.deg, 2*u.deg,
                   pm_ra_cosdec=.2*u.mas/u.yr, pm_dec=.1*u.mas/u.yr,
                   frame='fk5')
    assert_quantity_allclose(sc1.ra, 1*u.deg)
    assert_quantity_allclose(sc1.dec, 2*u.deg)
    assert_quantity_allclose(sc1.pm_ra_cosdec, .2*u.arcsec/u.kyr)
    assert_quantity_allclose(sc1.pm_dec, .1*u.arcsec/u.kyr)

    sc2 = SkyCoord(1*u.deg, 2*u.deg,
                   pm_ra=.2*u.mas/u.yr, pm_dec=.1*u.mas/u.yr,
                   differential_cls=SphericalDifferential)
    assert_quantity_allclose(sc2.ra, 1*u.deg)
    assert_quantity_allclose(sc2.dec, 2*u.deg)
    assert_quantity_allclose(sc2.pm_ra, .2*u.arcsec/u.kyr)
    assert_quantity_allclose(sc2.pm_dec, .1*u.arcsec/u.kyr)

    sc3 = SkyCoord('1:2:3 4:5:6',
                   pm_ra_cosdec=.2*u.mas/u.yr, pm_dec=.1*u.mas/u.yr,
                   unit=(u.hour, u.deg))

    assert_quantity_allclose(sc3.ra, 1*u.hourangle + 2*u.arcmin*15 + 3*u.arcsec*15)
    assert_quantity_allclose(sc3.dec, 4*u.deg + 5*u.arcmin + 6*u.arcsec)
    # might as well check with sillier units?
    assert_quantity_allclose(sc3.pm_ra_cosdec, 1.2776637006616473e-07 * u.arcmin / u.fortnight)
    assert_quantity_allclose(sc3.pm_dec, 6.388318503308237e-08 * u.arcmin / u.fortnight)

@pytest.mark.xfail  #TODO: FIX
def test_creation_copy():
    i = ICRS(1*u.deg, 2*u.deg, pm_ra_cosdec=.2*u.mas/u.yr, pm_dec=.1*u.mas/u.yr)
    sc = SkyCoord(i)
    sc_cpy = SkyCoord(sc)

    for attrnm in ['ra', 'dec', 'pm_ra_cosdec', 'pm_dec']:
        assert_quantity_allclose(getattr(sc, attrnm), getattr(sc_cpy, attrnm))

    sc2 = SkyCoord(1*u.deg, 2*u.deg,
                   pm_ra=.2*u.mas/u.yr, pm_dec=.1*u.mas/u.yr,
                   differential_cls=SphericalDifferential)

    sc2_cpy = SkyCoord(sc2)
    for attrnm in ['ra', 'dec', 'pm_ra', 'pm_dec']:
        assert_quantity_allclose(getattr(sc2, attrnm), getattr(sc2_cpy, attrnm))

    sc2_newdiff = SkyCoord(sc2, differential_cls=SphericalCosLatDifferential)
    reprepr = sc2.represent_as(SphericalRepresentation, SphericalCosLatDifferential)
    assert_quantity_allclose(sc2_newdiff.pm_ra_cosdec, reprepr.d_lon_coslat)


@pytest.mark.xfail  #TODO: FIX
def test_useful_error_missing():
    sc_nod = SkyCoord(ICRS(1*u.deg, 2*u.deg))
    try:
        sc_nod.l
    except AttributeError as e:
        # this is double-checking the *normal* behavior
        msg_l = e.args[0]

    try:
        sc_nod.pm_dec
    except AttributeError as e:
        msg_pm_dec = e.args[0]

    assert "has no attribute" in msg_l
    assert "has no attribute" in msg_pm_dec


# ----------------------Operations on SkyCoords w/ velocities-------------------

# define some fixtires to get baseline coordinates to try operations with
@pytest.fixture(scope="module")
def sc():
    return SkyCoord(1*u.deg, 2*u.deg,
                  pm_dec=1*u.mas/u.yr, pm_ra_cosdec=2*u.mas/u.yr)
@pytest.fixture(scope="module")
def sc2():
    return SkyCoord(1*u.deg, 2*u.deg,
                pm_dec=1*u.mas/u.yr, pm_ra_cosdec=2*u.mas/u.yr)
@pytest.fixture(scope="module")
def scmany():
    return SkyCoord(ICRS(ra=[1]*100*u.deg, dec=[2]*100*u.deg,
                     pm_ra_cosdec=np.random.randn(100)*u.mas/u.yr,
                     pm_dec=np.random.randn(100)*u.mas/u.yr,))


def test_separation(sc, sc2):
    sc.separation(sc2)


def test_accessors(sc, scmany):
    sc.data.differentials['s']
    sc.spherical
    sc.galactic

    scmany[0]
    scmany.spherical
    scmany.galactic


def test_transforms(sc):
    trans = sc.transform_to('galactic')
    assert isinstance(trans.frame, Galactic)


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_matching(sc, scmany):
    # just check that it works and yields something
    idx, d2d, d3d = sc.match_to_catalog_sky(scmany)


def test_position_angle(sc, sc2):
    sc.position_angle(sc2)


@pytest.mark.xfail  #TODO: FIX
def test_constellations(sc):
    const = sc.get_constellation()
    assert const == 'Pisces'


@pytest.mark.xfail  #TODO: FIX
def test_skyoffset_frame(sc):
    sc.skyoffset_frame()

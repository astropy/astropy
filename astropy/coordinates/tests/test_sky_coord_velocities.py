# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for putting velocity differentials into SkyCoord objects.

Note: the skyoffset velocity tests are in a different file, in
test_skyoffset_transformations.py
"""

import pytest

import numpy as np

from ... import units as u
from ...tests.helper import assert_quantity_allclose
from .. import (SkyCoord, ICRS, SphericalRepresentation, SphericalDifferential,
                SphericalCosLatDifferential, CartesianRepresentation,
                CartesianDifferential, Galactic, PrecessedGeocentric)

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
                   differential_type=SphericalDifferential)
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


def test_creation_copy_basic():
    i = ICRS(1*u.deg, 2*u.deg, pm_ra_cosdec=.2*u.mas/u.yr, pm_dec=.1*u.mas/u.yr)
    sc = SkyCoord(i)
    sc_cpy = SkyCoord(sc)

    for attrnm in ['ra', 'dec', 'pm_ra_cosdec', 'pm_dec']:
        assert_quantity_allclose(getattr(sc, attrnm), getattr(sc_cpy, attrnm))


def test_creation_copy_rediff():
    sc = SkyCoord(1*u.deg, 2*u.deg,
                  pm_ra=.2*u.mas/u.yr, pm_dec=.1*u.mas/u.yr,
                  differential_type=SphericalDifferential)

    sc_cpy = SkyCoord(sc)
    for attrnm in ['ra', 'dec', 'pm_ra', 'pm_dec']:
        assert_quantity_allclose(getattr(sc, attrnm), getattr(sc_cpy, attrnm))

    sc_newdiff = SkyCoord(sc, differential_type=SphericalCosLatDifferential)
    reprepr = sc.represent_as(SphericalRepresentation, SphericalCosLatDifferential)
    assert_quantity_allclose(sc_newdiff.pm_ra_cosdec,
                             reprepr.differentials['s'].d_lon_coslat)


def test_creation_cartesian():
    rep = CartesianRepresentation([10, 0., 0.]*u.pc)
    dif = CartesianDifferential([0, 100, 0.]*u.pc/u.Myr)
    rep = rep.with_differentials(dif)
    c = SkyCoord(rep)

    sdif = dif.represent_as(SphericalCosLatDifferential, rep)
    assert_quantity_allclose(c.pm_ra_cosdec, sdif.d_lon_coslat)


def test_useful_error_missing():
    sc_nod = SkyCoord(ICRS(1*u.deg, 2*u.deg))
    try:
        sc_nod.l
    except AttributeError as e:
        # this is double-checking the *normal* behavior
        msg_l = e.args[0]

    try:
        sc_nod.pm_dec
    except Exception as e:
        msg_pm_dec = e.args[0]

    assert "has no attribute" in msg_l
    assert "has no associated differentials" in msg_pm_dec


# ----------------------Operations on SkyCoords w/ velocities-------------------

# define some fixtures to get baseline coordinates to try operations with
@pytest.fixture(scope="module", params=[(False, False),
                                        (True, False),
                                        (False, True),
                                        (True, True)])
def sc(request):
    incldist, inclrv = request.param

    args = [1*u.deg, 2*u.deg]
    kwargs = dict(pm_dec=1*u.mas/u.yr, pm_ra_cosdec=2*u.mas/u.yr)
    if incldist:
        kwargs['distance'] = 213.4*u.pc
    if inclrv:
        kwargs['radial_velocity'] = 61*u.km/u.s

    return SkyCoord(*args, **kwargs)

@pytest.fixture(scope="module")
def scmany():
    return SkyCoord(ICRS(ra=[1]*100*u.deg, dec=[2]*100*u.deg,
                     pm_ra_cosdec=np.random.randn(100)*u.mas/u.yr,
                     pm_dec=np.random.randn(100)*u.mas/u.yr,))

@pytest.fixture(scope="module")
def sc_for_sep():
    return SkyCoord(1*u.deg, 2*u.deg,
                    pm_dec=1*u.mas/u.yr, pm_ra_cosdec=2*u.mas/u.yr)


def test_separation(sc, sc_for_sep):
    sc.separation(sc_for_sep)


def test_accessors(sc, scmany):
    sc.data.differentials['s']
    sph = sc.spherical
    gal = sc.galactic

    if (sc.data.get_name().startswith('unit') and not
        sc.data.differentials['s'].get_name().startswith('unit')):
        # this xfail can be eliminated when issue #7028 is resolved
        pytest.xfail('.velocity fails if there is an RV but not distance')
    sc.velocity

    assert isinstance(sph, SphericalRepresentation)
    assert gal.data.differentials is not None

    scmany[0]
    sph = scmany.spherical
    gal = scmany.galactic

    assert isinstance(sph, SphericalRepresentation)
    assert gal.data.differentials is not None

def test_transforms(sc):
    trans = sc.transform_to('galactic')
    assert isinstance(trans.frame, Galactic)


def test_transforms_diff(sc):
    # note that arguably this *should* fail for the no-distance cases: 3D
    # information is necessary to truly solve this, hence the xfail
    if not sc.distance.unit.is_equivalent(u.m):
        pytest.xfail('Should fail for no-distance cases')
    else:
        trans = sc.transform_to(PrecessedGeocentric(equinox='B1975'))
        assert isinstance(trans.frame, PrecessedGeocentric)


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_matching(sc, scmany):
    # just check that it works and yields something
    idx, d2d, d3d = sc.match_to_catalog_sky(scmany)


def test_position_angle(sc, sc_for_sep):
    sc.position_angle(sc_for_sep)


def test_constellations(sc):
    const = sc.get_constellation()
    assert const == 'Pisces'


def test_separation_3d_with_differentials():
    c1 = SkyCoord(ra=138*u.deg, dec=-17*u.deg, distance=100*u.pc,
                  pm_ra_cosdec=5*u.mas/u.yr,
                  pm_dec=-7*u.mas/u.yr,
                  radial_velocity=160*u.km/u.s)
    c2 = SkyCoord(ra=138*u.deg, dec=-17*u.deg, distance=105*u.pc,
                  pm_ra_cosdec=15*u.mas/u.yr,
                  pm_dec=-74*u.mas/u.yr,
                  radial_velocity=-60*u.km/u.s)

    sep = c1.separation_3d(c2)
    assert_quantity_allclose(sep, 5*u.pc)

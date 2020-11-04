# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Accuracy tests for Ecliptic coordinate systems.
"""

import numpy as np

import pytest

from astropy.units import allclose as quantity_allclose
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.builtin_frames import FK5, ICRS, GCRS, GeocentricMeanEcliptic, BarycentricMeanEcliptic, HeliocentricMeanEcliptic, GeocentricTrueEcliptic, BarycentricTrueEcliptic, HeliocentricTrueEcliptic, HeliocentricEclipticIAU76
from astropy.coordinates.solar_system import get_body_barycentric_posvel
from astropy.constants import R_sun, R_earth
from astropy.time import Time


def test_against_pytpm_doc_example():
    """
    Check that Astropy's Ecliptic systems give answers consistent with pyTPM

    Currently this is only testing against the example given in the pytpm docs
    """
    fk5_in = SkyCoord('12h22m54.899s', '15d49m20.57s', frame=FK5(equinox='J2000'))
    pytpm_out = BarycentricMeanEcliptic(lon=178.78256462*u.deg,
                                        lat=16.7597002513*u.deg,
                                        equinox='J2000')
    astropy_out = fk5_in.transform_to(pytpm_out)

    assert pytpm_out.separation(astropy_out) < (1*u.arcsec)


def test_ecliptic_heliobary():
    """
    Check that the ecliptic transformations for heliocentric and barycentric
    at least more or less make sense
    """
    icrs = ICRS(1*u.deg, 2*u.deg, distance=1.5*R_sun)

    bary = icrs.transform_to(BarycentricMeanEcliptic())
    helio = icrs.transform_to(HeliocentricMeanEcliptic())

    # make sure there's a sizable distance shift - in 3d hundreds of km, but
    # this is 1D so we allow it to be somewhat smaller
    assert np.abs(bary.distance - helio.distance) > 1*u.km

    # now make something that's got the location of helio but in bary's frame.
    # this is a convenience to allow `separation` to work as expected
    helio_in_bary_frame = bary.realize_frame(helio.cartesian)
    assert bary.separation(helio_in_bary_frame) > 1*u.arcmin


@pytest.mark.parametrize(('trueframe', 'meanframe'),
                        [(BarycentricTrueEcliptic, BarycentricMeanEcliptic),
                         (HeliocentricTrueEcliptic, HeliocentricMeanEcliptic),
                         (GeocentricTrueEcliptic, GeocentricMeanEcliptic),
                         (HeliocentricEclipticIAU76, HeliocentricMeanEcliptic)])
def test_ecliptic_roundtrips(trueframe, meanframe):
    """
    Check that the various ecliptic transformations at least roundtrip
    """
    icrs = ICRS(1*u.deg, 2*u.deg, distance=1.5*R_sun)

    truecoo = icrs.transform_to(trueframe())
    meancoo = truecoo.transform_to(meanframe())
    truecoo2 = meancoo.transform_to(trueframe())

    assert not quantity_allclose(truecoo.cartesian.xyz, meancoo.cartesian.xyz)
    assert quantity_allclose(truecoo.cartesian.xyz, truecoo2.cartesian.xyz)


@pytest.mark.parametrize(('trueframe', 'meanframe'),
                        [(BarycentricTrueEcliptic, BarycentricMeanEcliptic),
                         (HeliocentricTrueEcliptic, HeliocentricMeanEcliptic),
                         (GeocentricTrueEcliptic, GeocentricMeanEcliptic)])
def test_ecliptic_true_mean_preserve_latitude(trueframe, meanframe):
    """
    Check that the ecliptic true/mean transformations preserve latitude
    """
    truecoo = trueframe(90*u.deg, 0*u.deg, distance=1*u.AU)
    meancoo = truecoo.transform_to(meanframe())

    assert not quantity_allclose(truecoo.lon, meancoo.lon)
    assert quantity_allclose(truecoo.lat, meancoo.lat, atol=1e-10*u.arcsec)


@pytest.mark.parametrize('frame',
                         [HeliocentricMeanEcliptic,
                          HeliocentricTrueEcliptic,
                          HeliocentricEclipticIAU76])
def test_helioecliptic_induced_velocity(frame):
    # Create a coordinate with zero speed in ICRS
    time = Time('2021-01-01')
    icrs = ICRS(ra=1*u.deg, dec=2*u.deg, distance=3*u.AU,
                pm_ra_cosdec=0*u.deg/u.s, pm_dec=0*u.deg/u.s, radial_velocity=0*u.m/u.s)

    # Transforming to a helioecliptic frame should give an induced speed equal to the Sun's speed
    transformed = icrs.transform_to(frame(obstime=time))
    _, vel = get_body_barycentric_posvel('sun', time)
    assert quantity_allclose(transformed.velocity.norm(), vel.norm())

    # Transforming back to ICRS should get back to zero speed
    back = transformed.transform_to(ICRS())
    assert quantity_allclose(back.velocity.norm(), 0*u.m/u.s, atol=1e-10*u.m/u.s)


def test_ecl_geo():
    """
    Check that the geocentric version at least gets well away from GCRS.  For a
    true "accuracy" test we need a comparison dataset that is similar to the
    geocentric/GCRS comparison we want to do here.  Contributions welcome!
    """
    gcrs = GCRS(10*u.deg, 20*u.deg, distance=1.5*R_earth)
    gecl = gcrs.transform_to(GeocentricMeanEcliptic())

    assert quantity_allclose(gecl.distance, gcrs.distance)


def test_arraytransforms():
    """
    Test that transforms to/from ecliptic coordinates work on array coordinates
    (not testing for accuracy.)
    """
    ra = np.ones((4, ), dtype=float) * u.deg
    dec = 2*np.ones((4, ), dtype=float) * u.deg
    distance = np.ones((4, ), dtype=float) * u.au

    test_icrs = ICRS(ra=ra, dec=dec, distance=distance)
    test_gcrs = GCRS(test_icrs.data)

    bary_arr = test_icrs.transform_to(BarycentricMeanEcliptic())
    assert bary_arr.shape == ra.shape

    helio_arr = test_icrs.transform_to(HeliocentricMeanEcliptic())
    assert helio_arr.shape == ra.shape

    geo_arr = test_gcrs.transform_to(GeocentricMeanEcliptic())
    assert geo_arr.shape == ra.shape

    # now check that we also can go back the other way without shape problems
    bary_icrs = bary_arr.transform_to(ICRS())
    assert bary_icrs.shape == test_icrs.shape

    helio_icrs = helio_arr.transform_to(ICRS())
    assert helio_icrs.shape == test_icrs.shape

    geo_gcrs = geo_arr.transform_to(GCRS())
    assert geo_gcrs.shape == test_gcrs.shape


def test_roundtrip_scalar():
    icrs = ICRS(ra=1*u.deg, dec=2*u.deg, distance=3*u.au)
    gcrs = GCRS(icrs.cartesian)

    bary = icrs.transform_to(BarycentricMeanEcliptic())
    helio = icrs.transform_to(HeliocentricMeanEcliptic())
    geo = gcrs.transform_to(GeocentricMeanEcliptic())

    bary_icrs = bary.transform_to(ICRS())
    helio_icrs = helio.transform_to(ICRS())
    geo_gcrs = geo.transform_to(GCRS())

    assert quantity_allclose(bary_icrs.cartesian.xyz, icrs.cartesian.xyz)
    assert quantity_allclose(helio_icrs.cartesian.xyz, icrs.cartesian.xyz)
    assert quantity_allclose(geo_gcrs.cartesian.xyz, gcrs.cartesian.xyz)

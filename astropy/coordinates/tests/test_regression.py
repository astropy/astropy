# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Regression tests for coordinates-related bugs that don't have an obvious other
place to live
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ... import units as u
from .. import (AltAz, EarthLocation, SkyCoord, get_sun, ICRS, CIRS, ITRS,
                GeocentricTrueEcliptic, Longitude, Latitude, GCRS, HCRS,
                get_moon, FK4, FK4NoETerms)
from ..sites import get_builtin_sites
from ...time import Time
from ...utils import iers

from ...tests.helper import pytest, assert_quantity_allclose, catch_warnings, quantity_allclose
from .test_matching import HAS_SCIPY, OLDER_SCIPY


def test_regression_5085():
    """
    PR #5085 was put in place to fix the following issue.

    Issue: https://github.com/astropy/astropy/issues/5069
    At root was the transformation of Ecliptic coordinates with
    non-scalar times.
    """
    times = Time(["2015-08-28 03:30", "2015-09-05 10:30", "2015-09-15 18:35"])
    latitudes = Latitude([3.9807075, -5.00733806, 1.69539491]*u.deg)
    longitudes = Longitude([311.79678613,  72.86626741, 199.58698226]*u.deg)
    distances = u.Quantity([0.00243266, 0.0025424, 0.00271296]*u.au)
    coo = GeocentricTrueEcliptic(lat=latitudes,
                                 lon=longitudes,
                                 distance=distances, equinox=times)
    # expected result
    ras = Longitude([310.50095387, 314.67109863, 319.56507471]*u.deg)
    decs = Latitude([-18.25190707, -17.1556641, -15.71616651]*u.deg)
    distances = u.Quantity([1.78309902, 1.710874, 1.61326648]*u.au)
    expected_result = GCRS(ra=ras, dec=decs,
                           distance=distances, obstime="J2000").cartesian.xyz
    actual_result = coo.transform_to(GCRS(obstime="J2000")).cartesian.xyz
    assert_quantity_allclose(expected_result, actual_result)


def test_regression_3920():
    """
    Issue: https://github.com/astropy/astropy/issues/3920
    """
    loc = EarthLocation.from_geodetic(0*u.deg, 0*u.deg, 0)
    time = Time('2010-1-1')

    aa = AltAz(location=loc, obstime=time)
    sc = SkyCoord(10*u.deg, 3*u.deg)
    assert sc.transform_to(aa).shape == tuple()
    # That part makes sense: the input is a scalar so the output is too

    sc2 = SkyCoord(10*u.deg, 3*u.deg, 1*u.AU)
    assert sc2.transform_to(aa).shape == tuple()
    # in 3920 that assert fails, because the shape is (1,)

    # check that the same behavior occurs even if transform is from low-level classes
    icoo = ICRS(sc.data)
    icoo2 = ICRS(sc2.data)
    assert icoo.transform_to(aa).shape == tuple()
    assert icoo2.transform_to(aa).shape == tuple()


def test_regression_3938():
    """
    Issue: https://github.com/astropy/astropy/issues/3938
    """
    # Set up list of targets - we don't use `from_name` here to avoid
    # remote_data requirements, but it does the same thing
    # vega = SkyCoord.from_name('Vega')
    vega = SkyCoord(279.23473479*u.deg, 38.78368896*u.deg)
    # capella = SkyCoord.from_name('Capella')
    capella = SkyCoord(79.17232794*u.deg, 45.99799147*u.deg)
    # sirius = SkyCoord.from_name('Sirius')
    sirius = SkyCoord(101.28715533*u.deg, -16.71611586*u.deg)
    targets = [vega, capella, sirius]

    # Feed list of targets into SkyCoord
    combined_coords = SkyCoord(targets)

    # Set up AltAz frame
    time = Time('2012-01-01 00:00:00')
    location = EarthLocation('10d', '45d', 0)
    aa = AltAz(location=location, obstime=time)

    combined_coords.transform_to(aa)
    # in 3938 the above yields ``UnitConversionError: '' (dimensionless) and 'pc' (length) are not convertible``


def test_regression_3998():
    """
    Issue: https://github.com/astropy/astropy/issues/3998
    """
    time = Time('2012-01-01 00:00:00')
    assert time.isscalar

    sun = get_sun(time)
    assert sun.isscalar
    # in 3998, the above yields False - `sun` is a length-1 vector

    assert sun.obstime is time


def test_regression_4033():
    """
    Issue: https://github.com/astropy/astropy/issues/4033
    """
    # alb = SkyCoord.from_name('Albireo')
    alb = SkyCoord(292.68033548*u.deg, 27.95968007*u.deg)
    alb_wdist = SkyCoord(alb, distance=133*u.pc)

    # de = SkyCoord.from_name('Deneb')
    de = SkyCoord(310.35797975*u.deg, 45.28033881*u.deg)
    de_wdist = SkyCoord(de, distance=802*u.pc)

    aa = AltAz(location=EarthLocation(lat=45*u.deg, lon=0*u.deg), obstime='2010-1-1')
    deaa = de.transform_to(aa)
    albaa = alb.transform_to(aa)
    alb_wdistaa = alb_wdist.transform_to(aa)
    de_wdistaa = de_wdist.transform_to(aa)

    # these work fine
    sepnod = deaa.separation(albaa)
    sepwd = deaa.separation(alb_wdistaa)
    assert_quantity_allclose(sepnod, 22.2862*u.deg, rtol=1e-6)
    assert_quantity_allclose(sepwd, 22.2862*u.deg, rtol=1e-6)
    # parallax should be present when distance added
    assert np.abs(sepnod - sepwd) > 1*u.marcsec

    # in 4033, the following fail with a recursion error
    assert_quantity_allclose(de_wdistaa.separation(alb_wdistaa), 22.2862*u.deg, rtol=1e-3)
    assert_quantity_allclose(alb_wdistaa.separation(deaa), 22.2862*u.deg, rtol=1e-3)


@pytest.mark.skipif(not HAS_SCIPY, reason='No Scipy')
@pytest.mark.skipif(OLDER_SCIPY, reason='Scipy too old')
def test_regression_4082():
    """
    Issue: https://github.com/astropy/astropy/issues/4082
    """
    from .. import search_around_sky, search_around_3d
    cat = SkyCoord([10.076, 10.00455], [18.54746, 18.54896], unit='deg')
    search_around_sky(cat[0:1], cat, seplimit=u.arcsec * 60, storekdtree=False)
    # in the issue, this raises a TypeError

    # also check 3d for good measure, although it's not really affected by this bug directly
    cat3d = SkyCoord([10.076, 10.00455]*u.deg, [18.54746, 18.54896]*u.deg, distance=[0.1, 1.5]*u.kpc)
    search_around_3d(cat3d[0:1], cat3d, 1*u.kpc, storekdtree=False)


def test_regression_4210():
    """
    Issue: https://github.com/astropy/astropy/issues/4210
    Related PR with actual change: https://github.com/astropy/astropy/pull/4211
    """
    crd = SkyCoord(0*u.deg, 0*u.deg, distance=1*u.AU)
    ecl = crd.geocentrictrueecliptic
    # bug was that "lambda", which at the time was the name of the geocentric
    # ecliptic longitude, is a reserved keyword. So this just makes sure the
    # new name is are all valid
    ecl.lon

    # and for good measure, check the other ecliptic systems are all the same
    # names for their attributes
    from ..builtin_frames import ecliptic
    for frame_name in ecliptic.__all__:
        eclcls = getattr(ecliptic, frame_name)
        eclobj = eclcls(1*u.deg, 2*u.deg, 3*u.AU)

        eclobj.lat
        eclobj.lon
        eclobj.distance


def test_regression_futuretimes_4302():
    """
    Checks that an error is not raised for future times not covered by IERS
    tables (at least in a simple transform like CIRS->ITRS that simply requires
    the UTC<->UT1 conversion).

    Relevant comment: https://github.com/astropy/astropy/pull/4302#discussion_r44836531
    """
    from ...utils.exceptions import AstropyWarning

    # this is an ugly hack to get the warning to show up even if it has already
    # appeared
    from ..builtin_frames import utils
    if hasattr(utils, '__warningregistry__'):
        utils.__warningregistry__.clear()

    with catch_warnings() as found_warnings:
        future_time = Time('2511-5-1')
        c = CIRS(1*u.deg, 2*u.deg, obstime=future_time)
        c.transform_to(ITRS(obstime=future_time))

    if not isinstance(iers.IERS_Auto.iers_table, iers.IERS_Auto):
        saw_iers_warnings = False
        for w in found_warnings:
            if issubclass(w.category, AstropyWarning):
                if '(some) times are outside of range covered by IERS table' in str(w.message):
                    saw_iers_warnings = True
                    break
        assert saw_iers_warnings, 'Never saw IERS warning'


def test_regression_4996():
    # this part is the actual regression test
    deltat = np.linspace(-12, 12, 1000)*u.hour
    times = Time('2012-7-13 00:00:00') + deltat
    suncoo = get_sun(times)
    assert suncoo.shape == (len(times),)

    # and this is an additional test to make sure more complex arrays work
    times2 = Time('2012-7-13 00:00:00') + deltat.reshape(10, 20, 5)
    suncoo2 = get_sun(times2)
    assert suncoo2.shape == times2.shape

    # this is intentionally not allclose - they should be *exactly* the same
    assert np.all(suncoo.ra.ravel() == suncoo2.ra.ravel())


def test_regression_4293():
    """Really just an extra test on FK4 no e, after finding that the units
    were not always taken correctly.  This test is against explicitly doing
    the transformations on pp170 of Explanatory Supplement to the Astronomical
    Almanac (Seidelmann, 2005).

    See https://github.com/astropy/astropy/pull/4293#issuecomment-234973086
    """
    # Check all over sky, but avoiding poles (note that FK4 did not ignore
    # e terms within 10∘ of the poles...  see p170 of explan.supp.).
    ra, dec = np.meshgrid(np.arange(0, 359, 45), np.arange(-80, 81, 40))
    fk4 = FK4(ra.ravel() * u.deg, dec.ravel() * u.deg)

    Dc = -0.065838*u.arcsec
    Dd = +0.335299*u.arcsec
    # Dc * tan(obliquity), as given on p.170
    Dctano = -0.028553*u.arcsec

    fk4noe_dec = (fk4.dec - (Dd*np.cos(fk4.ra) -
                             Dc*np.sin(fk4.ra))*np.sin(fk4.dec) -
                  Dctano*np.cos(fk4.dec))
    fk4noe_ra = fk4.ra - (Dc*np.cos(fk4.ra) +
                          Dd*np.sin(fk4.ra)) / np.cos(fk4.dec)

    fk4noe = fk4.transform_to(FK4NoETerms)
    # Tolerance here just set to how well the coordinates match, which is much
    # better than the claimed accuracy of <1 mas for this first-order in
    # v_earth/c approximation.
    # Interestingly, if one divides by np.cos(fk4noe_dec) in the ra correction,
    # the match becomes good to 2 μas.
    assert_quantity_allclose(fk4noe.ra, fk4noe_ra, atol=11.*u.uas, rtol=0)
    assert_quantity_allclose(fk4noe.dec, fk4noe_dec, atol=3.*u.uas, rtol=0)


def test_regression_4926():
    times = Time('2010-01-1') + np.arange(20)*u.day
    green = get_builtin_sites()['greenwich']
    # this is the regression test
    moon = get_moon(times, green)

    # this is an additional test to make sure the GCRS->ICRS transform works for complex shapes
    moon.transform_to(ICRS())

    # and some others to increase coverage of transforms
    moon.transform_to(HCRS(obstime="J2000"))
    moon.transform_to(HCRS(obstime=times))


def test_regression_5209():
    "check that distances are not lost on SkyCoord init"
    time = Time('2015-01-01')
    moon = get_moon(time)
    new_coord = SkyCoord([moon])
    assert_quantity_allclose(new_coord[0].distance, moon.distance)


def test_regression_5133():
    N = 1000
    np.random.seed(12345)
    lon = np.random.uniform(-10, 10, N) * u.deg
    lat = np.random.uniform(50, 52, N) * u.deg
    alt = np.random.uniform(0, 10., N) * u.km

    time = Time('2010-1-1')

    objects = EarthLocation.from_geodetic(lon, lat, height=alt)
    itrs_coo = objects.get_itrs(time)

    homes = [EarthLocation.from_geodetic(lon=-1 * u.deg, lat=52 * u.deg, height=h)
             for h in (0, 1000, 10000)*u.km]

    altaz_frames = [AltAz(obstime=time, location=h) for h in homes]
    altaz_coos = [itrs_coo.transform_to(f) for f in altaz_frames]

    # they should all be different
    for coo in altaz_coos[1:]:
        assert not quantity_allclose(coo.az, coo.az[0])
        assert not quantity_allclose(coo.alt, coo.alt[0])


def test_itrs_vals_5133():
    time = Time('2010-1-1')
    el = EarthLocation.from_geodetic(lon=20*u.deg, lat=45*u.deg, height=0*u.km)

    lons = [20, 30, 20]*u.deg
    lats = [44, 45, 45]*u.deg
    alts = [0, 0, 10]*u.km
    coos = [EarthLocation.from_geodetic(lon, lat, height=alt).get_itrs(time)
            for lon, lat, alt in zip(lons, lats, alts)]

    aaf = AltAz(obstime=time, location=el)
    aacs = [coo.transform_to(aaf) for coo in coos]

    assert all([coo.isscalar for coo in aacs])

    # the ~1 arcsec tolerance is b/c aberration makes it not exact
    assert_quantity_allclose(aacs[0].az, 180*u.deg, atol=1*u.arcsec)
    assert aacs[0].alt < 0*u.deg
    assert aacs[0].distance > 50*u.km

    # it should *not* actually be 90 degrees, b/c constant latitude is not
    # straight east anywhere except the equator... but should be close-ish
    assert_quantity_allclose(aacs[1].az, 90*u.deg, atol=5*u.deg)
    assert aacs[1].alt < 0*u.deg
    assert aacs[1].distance > 50*u.km

    assert_quantity_allclose(aacs[2].alt, 90*u.deg, atol=1*u.arcsec)
    assert_quantity_allclose(aacs[2].distance, 10*u.km)


def test_regression_simple_5133():
    t = Time('J2010')
    obj = EarthLocation(-1*u.deg, 52*u.deg, height=[100., 0.]*u.km)
    home = EarthLocation(-1*u.deg, 52*u.deg, height=10.*u.km)
    aa = obj.get_itrs(t).transform_to(AltAz(obstime=t, location=home))

    # az is more-or-less undefined for straight up or down
    assert_quantity_allclose(aa.alt, [90, -90]*u.deg, rtol=1e-5)
    assert_quantity_allclose(aa.distance, [90, 10]*u.km)


def test_regression_5884():
    # ensure comparison of frames with non-scalar attributes works
    times = Time(["2015-08-28 03:30", "2015-08-28 12:00", "2015-09-05 10:30", "2015-09-15 18:35"])
    coo1 = SkyCoord(40*u.deg, 50*u.deg, obstime=times)
    coo2 = SkyCoord(41*u.deg, 52*u.deg, obstime=times)
    # check this does not raise a ValueError
    coo1.is_equivalent_frame(coo2)
    # let's do N-D whilst we're at it
    times = times.reshape((2,2))
    coo1 = SkyCoord(40*u.deg, 50*u.deg, obstime=times)
    coo2 = SkyCoord(41*u.deg, 52*u.deg, obstime=times)
    coo1.is_equivalent_frame(coo2)


def test_regression_5889_5890():
    # ensure we can represent all Representations and transform to ND frames
    greenwich = EarthLocation(
        *u.Quantity([3980608.90246817, -102.47522911,  4966861.27310067],
        unit=u.m))
    times = Time("2017-03-20T12:00:00") + np.linspace(-2,2,3)*u.hour
    moon = get_moon(times, location=greenwich)
    targets = SkyCoord([350.7*u.deg, 260.7*u.deg], [18.4*u.deg, 22.4*u.deg])
    targs2d = targets[:, np.newaxis]
    targs2d.transform_to(moon)

# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Regression tests for coordinates-related bugs that don't have an obvious other
place to live
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest
import numpy as np

from ...extern import six

from ... import units as u
from .. import (AltAz, EarthLocation, SkyCoord, get_sun, ICRS, CIRS, ITRS,
                GeocentricTrueEcliptic, Longitude, Latitude, GCRS, HCRS,
                get_moon, FK4, FK4NoETerms, BaseCoordinateFrame,
                QuantityAttribute, SphericalRepresentation,
                UnitSphericalRepresentation, CartesianRepresentation)
from ..sites import get_builtin_sites
from ...time import Time
from ...utils import iers
from ...table import Table

from ...tests.helper import assert_quantity_allclose, catch_warnings, quantity_allclose
from .test_matching import HAS_SCIPY, OLDER_SCIPY

try:
    import yaml  # pylint: disable=W0611
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


def test_regression_5085():
    """
    PR #5085 was put in place to fix the following issue.

    Issue: https://github.com/astropy/astropy/issues/5069
    At root was the transformation of Ecliptic coordinates with
    non-scalar times.
    """
    times = Time(["2015-08-28 03:30", "2015-09-05 10:30", "2015-09-15 18:35"])
    latitudes = Latitude([3.9807075, -5.00733806, 1.69539491]*u.deg)
    longitudes = Longitude([311.79678613, 72.86626741, 199.58698226]*u.deg)
    distances = u.Quantity([0.00243266, 0.0025424, 0.00271296]*u.au)
    coo = GeocentricTrueEcliptic(lat=latitudes,
                                 lon=longitudes,
                                 distance=distances, equinox=times)
    # expected result
    ras = Longitude([310.50095400, 314.67109920, 319.56507428]*u.deg)
    decs = Latitude([-18.25190443, -17.1556676, -15.71616522]*u.deg)
    distances = u.Quantity([1.78309901, 1.710874, 1.61326649]*u.au)
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


def test_regression_5743():
    sc = SkyCoord([5, 10], [20, 30], unit=u.deg,
                  obstime=['2017-01-01T00:00', '2017-01-01T00:10'])
    assert sc[0].obstime.shape == tuple()


def test_regression_5889_5890():
    # ensure we can represent all Representations and transform to ND frames
    greenwich = EarthLocation(
        *u.Quantity([3980608.90246817, -102.47522911, 4966861.27310067],
        unit=u.m))
    times = Time("2017-03-20T12:00:00") + np.linspace(-2, 2, 3)*u.hour
    moon = get_moon(times, location=greenwich)
    targets = SkyCoord([350.7*u.deg, 260.7*u.deg], [18.4*u.deg, 22.4*u.deg])
    targs2d = targets[:, np.newaxis]
    targs2d.transform_to(moon)


def test_regression_6236():
    # sunpy changes its representation upon initialisation of a frame,
    # including via `realize_frame`. Ensure this works.
    class MyFrame(BaseCoordinateFrame):
        default_representation = CartesianRepresentation
        my_attr = QuantityAttribute(default=0, unit=u.m)

    class MySpecialFrame(MyFrame):
        def __init__(self, *args, **kwargs):
            _rep_kwarg = kwargs.get('representation', None)
            super(MyFrame, self).__init__(*args, **kwargs)
            if not _rep_kwarg:
                self.representation = self.default_representation
                self._data = self.data.represent_as(self.representation)

    rep1 = UnitSphericalRepresentation([0., 1]*u.deg, [2., 3.]*u.deg)
    rep2 = SphericalRepresentation([10., 11]*u.deg, [12., 13.]*u.deg,
                                   [14., 15.]*u.kpc)
    mf1 = MyFrame(rep1, my_attr=1.*u.km)
    mf2 = mf1.realize_frame(rep2)
    # Normally, data is stored as is, but the representation gets set to a
    # default, even if a different representation instance was passed in.
    # realize_frame should do the same. Just in case, check attrs are passed.
    assert mf1.data is rep1
    assert mf2.data is rep2
    assert mf1.representation is CartesianRepresentation
    assert mf2.representation is CartesianRepresentation
    assert mf2.my_attr == mf1.my_attr
    # It should be independent of whether I set the reprensentation explicitly
    mf3 = MyFrame(rep1, my_attr=1.*u.km, representation='unitspherical')
    mf4 = mf3.realize_frame(rep2)
    assert mf3.data is rep1
    assert mf4.data is rep2
    assert mf3.representation is UnitSphericalRepresentation
    assert mf4.representation is CartesianRepresentation
    assert mf4.my_attr == mf3.my_attr
    # This should be enough to help sunpy, but just to be sure, a test
    # even closer to what is done there, i.e., transform the representation.
    msf1 = MySpecialFrame(rep1, my_attr=1.*u.km)
    msf2 = msf1.realize_frame(rep2)
    assert msf1.data is not rep1  # Gets transformed to Cartesian.
    assert msf2.data is not rep2
    assert type(msf1.data) is CartesianRepresentation
    assert type(msf2.data) is CartesianRepresentation
    assert msf1.representation is CartesianRepresentation
    assert msf2.representation is CartesianRepresentation
    assert msf2.my_attr == msf1.my_attr
    # And finally a test where the input is not transformed.
    msf3 = MySpecialFrame(rep1, my_attr=1.*u.km,
                          representation='unitspherical')
    msf4 = msf3.realize_frame(rep2)
    assert msf3.data is rep1
    assert msf4.data is not rep2
    assert msf3.representation is UnitSphericalRepresentation
    assert msf4.representation is CartesianRepresentation
    assert msf4.my_attr == msf3.my_attr


@pytest.mark.skipif(not HAS_SCIPY, reason='No Scipy')
@pytest.mark.skipif(OLDER_SCIPY, reason='Scipy too old')
def test_regression_6347():
    sc1 = SkyCoord([1, 2]*u.deg, [3, 4]*u.deg)
    sc2 = SkyCoord([1.1, 2.1]*u.deg, [3.1, 4.1]*u.deg)
    sc0 = sc1[:0]

    idx1_10, idx2_10, d2d_10, d3d_10 = sc1.search_around_sky(sc2, 10*u.arcmin)
    idx1_1, idx2_1, d2d_1, d3d_1 = sc1.search_around_sky(sc2, 1*u.arcmin)
    idx1_0, idx2_0, d2d_0, d3d_0 = sc0.search_around_sky(sc2, 10*u.arcmin)

    assert len(d2d_10) == 2

    assert len(d2d_0) == 0
    assert type(d2d_0) is type(d2d_10)

    assert len(d2d_1) == 0
    assert type(d2d_1) is type(d2d_10)


@pytest.mark.skipif(not HAS_SCIPY, reason='No Scipy')
@pytest.mark.skipif(OLDER_SCIPY, reason='Scipy too old')
def test_regression_6347_3d():
    sc1 = SkyCoord([1, 2]*u.deg, [3, 4]*u.deg, [5, 6]*u.kpc)
    sc2 = SkyCoord([1, 2]*u.deg, [3, 4]*u.deg, [5.1, 6.1]*u.kpc)
    sc0 = sc1[:0]

    idx1_10, idx2_10, d2d_10, d3d_10 = sc1.search_around_3d(sc2, 500*u.pc)
    idx1_1, idx2_1, d2d_1, d3d_1 = sc1.search_around_3d(sc2, 50*u.pc)
    idx1_0, idx2_0, d2d_0, d3d_0 = sc0.search_around_3d(sc2, 500*u.pc)

    assert len(d2d_10) > 0

    assert len(d2d_0) == 0
    assert type(d2d_0) is type(d2d_10)

    assert len(d2d_1) == 0
    assert type(d2d_1) is type(d2d_10)

def test_regression_6300():
    """Check that importing old frame attribute names from astropy.coordinates
    still works. See comments at end of #6300
    """
    from ...utils.exceptions import AstropyDeprecationWarning
    from .. import CartesianRepresentation
    from .. import (TimeFrameAttribute, QuantityFrameAttribute,
                    CartesianRepresentationFrameAttribute)

    with catch_warnings() as found_warnings:
        attr = TimeFrameAttribute(default=Time("J2000"))

        for w in found_warnings:
            if issubclass(w.category, AstropyDeprecationWarning):
                break
        else:
            assert False, "Deprecation warning not raised"

    with catch_warnings() as found_warnings:
        attr = QuantityFrameAttribute(default=5*u.km)

        for w in found_warnings:
            if issubclass(w.category, AstropyDeprecationWarning):
                break
        else:
            assert False, "Deprecation warning not raised"

    with catch_warnings() as found_warnings:
        attr = CartesianRepresentationFrameAttribute(
            default=CartesianRepresentation([5,6,7]*u.kpc))

        for w in found_warnings:
            if issubclass(w.category, AstropyDeprecationWarning):
                break
        else:
            assert False, "Deprecation warning not raised"


def test_gcrs_itrs_cartesian_repr():
    # issue 6436: transformation failed if coordinate representation was
    # Cartesian
    gcrs = GCRS(CartesianRepresentation((859.07256, -4137.20368,  5295.56871),
                                        unit='km'), representation='cartesian')
    gcrs.transform_to(ITRS)


@pytest.mark.skipif('not HAS_YAML')
def test_regression_6446():
    # this succeeds even before 6446:
    sc1 = SkyCoord([1, 2], [3, 4], unit='deg')
    t1 = Table([sc1])
    sio1 = six.StringIO()
    t1.write(sio1, format='ascii.ecsv')

    # but this fails due to the 6446 bug
    c1 = SkyCoord(1, 3, unit='deg')
    c2 = SkyCoord(2, 4, unit='deg')
    sc2 = SkyCoord([c1, c2])
    t2 = Table([sc2])
    sio2 = six.StringIO()
    t2.write(sio2, format='ascii.ecsv')

    assert sio1.getvalue() == sio2.getvalue()


def test_regression_6448():
    """
    This tests the more narrow problem reported in 6446 that 6448 is meant to
    fix. `test_regression_6446` also covers this, but this test is provided
    so that this is still tested even if YAML isn't installed.
    """
    sc1 = SkyCoord([1, 2], [3, 4], unit='deg')
    # this should always succeed even prior to 6448
    assert sc1.galcen_v_sun is None

    c1 = SkyCoord(1, 3, unit='deg')
    c2 = SkyCoord(2, 4, unit='deg')
    sc2 = SkyCoord([c1, c2])
    # without 6448 this fails
    assert sc2.galcen_v_sun is None


def test_regression_6597():
    frame_name = 'galactic'
    c1 = SkyCoord(1, 3, unit='deg', frame=frame_name)
    c2 = SkyCoord(2, 4, unit='deg', frame=frame_name)
    sc1 = SkyCoord([c1, c2])

    assert sc1.frame.name == frame_name


def test_regression_6597_2():
    """
    This tests the more subtle flaw that #6597 indirectly uncovered: that even
    in the case that the frames are ra/dec, they still might be the wrong *kind*
    """
    frame = FK4(equinox='J1949')
    c1 = SkyCoord(1, 3, unit='deg', frame=frame)
    c2 = SkyCoord(2, 4, unit='deg', frame=frame)
    sc1 = SkyCoord([c1, c2])

    assert sc1.frame.name == frame.name


def test_regression_6697():
    """
    Test for regression of a bug in get_gcrs_posvel that introduced errors at the 1m/s level.

    Comparison data is derived from calculation in PINT
    https://github.com/nanograv/PINT/blob/master/pint/erfautils.py
    """
    pint_vels = CartesianRepresentation(*(348.63632871, -212.31704928, -0.60154936), unit=u.m/u.s)
    location = EarthLocation(*(5327448.9957829, -1718665.73869569,  3051566.90295403), unit=u.m)
    t = Time(2458036.161966612, format='jd', scale='utc')
    obsgeopos, obsgeovel = location.get_gcrs_posvel(t)
    delta = (obsgeovel-pint_vels).norm()
    assert delta < 1*u.cm/u.s

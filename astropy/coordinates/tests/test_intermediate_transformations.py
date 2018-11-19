# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Accuracy tests for GCRS coordinate transformations, primarily to/from AltAz.

"""

import pytest
import numpy as np

from ... import units as u
from ...tests.helper import (assert_quantity_allclose as assert_allclose)
from ...time import Time
from .. import (EarthLocation, get_sun, ICRS, GCRS, CIRS, ITRS, AltAz,
                PrecessedGeocentric, CartesianRepresentation, SkyCoord,
                SphericalRepresentation, UnitSphericalRepresentation,
                HCRS, HeliocentricTrueEcliptic)


from ..._erfa import epv00

from .utils import randomly_sample_sphere
from ..builtin_frames.utils import get_jd12
from .. import solar_system_ephemeris
from ...units import allclose

try:
    import jplephem  # pylint: disable=W0611
except ImportError:
    HAS_JPLEPHEM = False
else:
    HAS_JPLEPHEM = True


def test_icrs_cirs():
    """
    Check a few cases of ICRS<->CIRS for consistency.

    Also includes the CIRS<->CIRS transforms at different times, as those go
    through ICRS
    """
    ra, dec, dist = randomly_sample_sphere(200)
    inod = ICRS(ra=ra, dec=dec)
    iwd = ICRS(ra=ra, dec=dec, distance=dist*u.pc)

    cframe1 = CIRS()
    cirsnod = inod.transform_to(cframe1)  # uses the default time
    # first do a round-tripping test
    inod2 = cirsnod.transform_to(ICRS)
    assert_allclose(inod.ra, inod2.ra)
    assert_allclose(inod.dec, inod2.dec)

    # now check that a different time yields different answers
    cframe2 = CIRS(obstime=Time('J2005', scale='utc'))
    cirsnod2 = inod.transform_to(cframe2)
    assert not allclose(cirsnod.ra, cirsnod2.ra, rtol=1e-8)
    assert not allclose(cirsnod.dec, cirsnod2.dec, rtol=1e-8)

    # parallax effects should be included, so with and w/o distance should be different
    cirswd = iwd.transform_to(cframe1)
    assert not allclose(cirswd.ra, cirsnod.ra, rtol=1e-8)
    assert not allclose(cirswd.dec, cirsnod.dec, rtol=1e-8)
    # and the distance should transform at least somehow
    assert not allclose(cirswd.distance, iwd.distance, rtol=1e-8)

    # now check that the cirs self-transform works as expected
    cirsnod3 = cirsnod.transform_to(cframe1)  # should be a no-op
    assert_allclose(cirsnod.ra, cirsnod3.ra)
    assert_allclose(cirsnod.dec, cirsnod3.dec)

    cirsnod4 = cirsnod.transform_to(cframe2)  # should be different
    assert not allclose(cirsnod4.ra, cirsnod.ra, rtol=1e-8)
    assert not allclose(cirsnod4.dec, cirsnod.dec, rtol=1e-8)

    cirsnod5 = cirsnod4.transform_to(cframe1)  # should be back to the same
    assert_allclose(cirsnod.ra, cirsnod5.ra)
    assert_allclose(cirsnod.dec, cirsnod5.dec)


ra, dec, dist = randomly_sample_sphere(200)
icrs_coords = [ICRS(ra=ra, dec=dec), ICRS(ra=ra, dec=dec, distance=dist*u.pc)]
gcrs_frames = [GCRS(), GCRS(obstime=Time('J2005', scale='utc'))]


@pytest.mark.parametrize('icoo', icrs_coords)
def test_icrs_gcrs(icoo):
    """
    Check ICRS<->GCRS for consistency
    """
    gcrscoo = icoo.transform_to(gcrs_frames[0])  # uses the default time
    # first do a round-tripping test
    icoo2 = gcrscoo.transform_to(ICRS)
    assert_allclose(icoo.distance, icoo2.distance)
    assert_allclose(icoo.ra, icoo2.ra)
    assert_allclose(icoo.dec, icoo2.dec)
    assert isinstance(icoo2.data, icoo.data.__class__)

    # now check that a different time yields different answers
    gcrscoo2 = icoo.transform_to(gcrs_frames[1])
    assert not allclose(gcrscoo.ra, gcrscoo2.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert not allclose(gcrscoo.dec, gcrscoo2.dec, rtol=1e-8, atol=1e-10*u.deg)

    # now check that the cirs self-transform works as expected
    gcrscoo3 = gcrscoo.transform_to(gcrs_frames[0])  # should be a no-op
    assert_allclose(gcrscoo.ra, gcrscoo3.ra)
    assert_allclose(gcrscoo.dec, gcrscoo3.dec)

    gcrscoo4 = gcrscoo.transform_to(gcrs_frames[1])  # should be different
    assert not allclose(gcrscoo4.ra, gcrscoo.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert not allclose(gcrscoo4.dec, gcrscoo.dec, rtol=1e-8, atol=1e-10*u.deg)

    gcrscoo5 = gcrscoo4.transform_to(gcrs_frames[0])  # should be back to the same
    assert_allclose(gcrscoo.ra, gcrscoo5.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert_allclose(gcrscoo.dec, gcrscoo5.dec, rtol=1e-8, atol=1e-10*u.deg)

    # also make sure that a GCRS with a different geoloc/geovel gets a different answer
    # roughly a moon-like frame
    gframe3 = GCRS(obsgeoloc=[385000., 0, 0]*u.km, obsgeovel=[1, 0, 0]*u.km/u.s)
    gcrscoo6 = icoo.transform_to(gframe3)  # should be different
    assert not allclose(gcrscoo.ra, gcrscoo6.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert not allclose(gcrscoo.dec, gcrscoo6.dec, rtol=1e-8, atol=1e-10*u.deg)
    icooviag3 = gcrscoo6.transform_to(ICRS)  # and now back to the original
    assert_allclose(icoo.ra, icooviag3.ra)
    assert_allclose(icoo.dec, icooviag3.dec)


@pytest.mark.parametrize('gframe', gcrs_frames)
def test_icrs_gcrs_dist_diff(gframe):
    """
    Check that with and without distance give different ICRS<->GCRS answers
    """
    gcrsnod = icrs_coords[0].transform_to(gframe)
    gcrswd = icrs_coords[1].transform_to(gframe)

    # parallax effects should be included, so with and w/o distance should be different
    assert not allclose(gcrswd.ra, gcrsnod.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert not allclose(gcrswd.dec, gcrsnod.dec, rtol=1e-8, atol=1e-10*u.deg)
    # and the distance should transform at least somehow
    assert not allclose(gcrswd.distance, icrs_coords[1].distance, rtol=1e-8,
                        atol=1e-10*u.pc)


def test_cirs_to_altaz():
    """
    Check the basic CIRS<->AltAz transforms.  More thorough checks implicitly
    happen in `test_iau_fullstack`
    """
    from .. import EarthLocation

    ra, dec, dist = randomly_sample_sphere(200)
    cirs = CIRS(ra=ra, dec=dec, obstime='J2000')
    crepr = SphericalRepresentation(lon=ra, lat=dec, distance=dist)
    cirscart = CIRS(crepr, obstime=cirs.obstime, representation_type=CartesianRepresentation)

    loc = EarthLocation(lat=0*u.deg, lon=0*u.deg, height=0*u.m)
    altazframe = AltAz(location=loc, obstime=Time('J2005'))

    cirs2 = cirs.transform_to(altazframe).transform_to(cirs)
    cirs3 = cirscart.transform_to(altazframe).transform_to(cirs)

    # check round-tripping
    assert_allclose(cirs.ra, cirs2.ra)
    assert_allclose(cirs.dec, cirs2.dec)
    assert_allclose(cirs.ra, cirs3.ra)
    assert_allclose(cirs.dec, cirs3.dec)


def test_gcrs_itrs():
    """
    Check basic GCRS<->ITRS transforms for round-tripping.
    """
    ra, dec, _ = randomly_sample_sphere(200)
    gcrs = GCRS(ra=ra, dec=dec, obstime='J2000')
    gcrs6 = GCRS(ra=ra, dec=dec, obstime='J2006')

    gcrs2 = gcrs.transform_to(ITRS).transform_to(gcrs)
    gcrs6_2 = gcrs6.transform_to(ITRS).transform_to(gcrs)

    assert_allclose(gcrs.ra, gcrs2.ra)
    assert_allclose(gcrs.dec, gcrs2.dec)
    assert not allclose(gcrs.ra, gcrs6_2.ra)
    assert not allclose(gcrs.dec, gcrs6_2.dec)

    # also try with the cartesian representation
    gcrsc = gcrs.realize_frame(gcrs.data)
    gcrsc.representation_type = CartesianRepresentation
    gcrsc2 = gcrsc.transform_to(ITRS).transform_to(gcrsc)
    assert_allclose(gcrsc.spherical.lon.deg, gcrsc2.ra.deg)
    assert_allclose(gcrsc.spherical.lat, gcrsc2.dec)


def test_cirs_itrs():
    """
    Check basic CIRS<->ITRS transforms for round-tripping.
    """
    ra, dec, _ = randomly_sample_sphere(200)
    cirs = CIRS(ra=ra, dec=dec, obstime='J2000')
    cirs6 = CIRS(ra=ra, dec=dec, obstime='J2006')

    cirs2 = cirs.transform_to(ITRS).transform_to(cirs)
    cirs6_2 = cirs6.transform_to(ITRS).transform_to(cirs)  # different obstime

    # just check round-tripping
    assert_allclose(cirs.ra, cirs2.ra)
    assert_allclose(cirs.dec, cirs2.dec)
    assert not allclose(cirs.ra, cirs6_2.ra)
    assert not allclose(cirs.dec, cirs6_2.dec)


def test_gcrs_cirs():
    """
    Check GCRS<->CIRS transforms for round-tripping.  More complicated than the
    above two because it's multi-hop
    """
    ra, dec, _ = randomly_sample_sphere(200)
    gcrs = GCRS(ra=ra, dec=dec, obstime='J2000')
    gcrs6 = GCRS(ra=ra, dec=dec, obstime='J2006')

    gcrs2 = gcrs.transform_to(CIRS).transform_to(gcrs)
    gcrs6_2 = gcrs6.transform_to(CIRS).transform_to(gcrs)

    assert_allclose(gcrs.ra, gcrs2.ra)
    assert_allclose(gcrs.dec, gcrs2.dec)
    assert not allclose(gcrs.ra, gcrs6_2.ra)
    assert not allclose(gcrs.dec, gcrs6_2.dec)

    # now try explicit intermediate pathways and ensure they're all consistent
    gcrs3 = gcrs.transform_to(ITRS).transform_to(CIRS).transform_to(ITRS).transform_to(gcrs)
    assert_allclose(gcrs.ra, gcrs3.ra)
    assert_allclose(gcrs.dec, gcrs3.dec)

    gcrs4 = gcrs.transform_to(ICRS).transform_to(CIRS).transform_to(ICRS).transform_to(gcrs)
    assert_allclose(gcrs.ra, gcrs4.ra)
    assert_allclose(gcrs.dec, gcrs4.dec)


def test_gcrs_altaz():
    """
    Check GCRS<->AltAz transforms for round-tripping.  Has multiple paths
    """
    from .. import EarthLocation

    ra, dec, _ = randomly_sample_sphere(1)
    gcrs = GCRS(ra=ra[0], dec=dec[0], obstime='J2000')

    # check array times sure N-d arrays work
    times = Time(np.linspace(2456293.25, 2456657.25, 51) * u.day,
                 format='jd', scale='utc')

    loc = EarthLocation(lon=10 * u.deg, lat=80. * u.deg)
    aaframe = AltAz(obstime=times, location=loc)

    aa1 = gcrs.transform_to(aaframe)
    aa2 = gcrs.transform_to(ICRS).transform_to(CIRS).transform_to(aaframe)
    aa3 = gcrs.transform_to(ITRS).transform_to(CIRS).transform_to(aaframe)

    # make sure they're all consistent
    assert_allclose(aa1.alt, aa2.alt)
    assert_allclose(aa1.az, aa2.az)
    assert_allclose(aa1.alt, aa3.alt)
    assert_allclose(aa1.az, aa3.az)


def test_precessed_geocentric():
    assert PrecessedGeocentric().equinox.jd == Time('J2000', scale='utc').jd

    gcrs_coo = GCRS(180*u.deg, 2*u.deg, distance=10000*u.km)
    pgeo_coo = gcrs_coo.transform_to(PrecessedGeocentric)
    assert np.abs(gcrs_coo.ra - pgeo_coo.ra) > 10*u.marcsec
    assert np.abs(gcrs_coo.dec - pgeo_coo.dec) > 10*u.marcsec
    assert_allclose(gcrs_coo.distance, pgeo_coo.distance)

    gcrs_roundtrip = pgeo_coo.transform_to(GCRS)
    assert_allclose(gcrs_coo.ra, gcrs_roundtrip.ra)
    assert_allclose(gcrs_coo.dec, gcrs_roundtrip.dec)
    assert_allclose(gcrs_coo.distance, gcrs_roundtrip.distance)

    pgeo_coo2 = gcrs_coo.transform_to(PrecessedGeocentric(equinox='B1850'))
    assert np.abs(gcrs_coo.ra - pgeo_coo2.ra) > 1.5*u.deg
    assert np.abs(gcrs_coo.dec - pgeo_coo2.dec) > 0.5*u.deg
    assert_allclose(gcrs_coo.distance, pgeo_coo2.distance)

    gcrs2_roundtrip = pgeo_coo2.transform_to(GCRS)
    assert_allclose(gcrs_coo.ra, gcrs2_roundtrip.ra)
    assert_allclose(gcrs_coo.dec, gcrs2_roundtrip.dec)
    assert_allclose(gcrs_coo.distance, gcrs2_roundtrip.distance)


# shared by parametrized tests below.  Some use the whole AltAz, others use just obstime
totest_frames = [AltAz(location=EarthLocation(-90*u.deg, 65*u.deg),
                       obstime=Time('J2000')),  # J2000 is often a default so this might work when others don't
                 AltAz(location=EarthLocation(120*u.deg, -35*u.deg),
                       obstime=Time('J2000')),
                 AltAz(location=EarthLocation(-90*u.deg, 65*u.deg),
                       obstime=Time('2014-01-01 00:00:00')),
                 AltAz(location=EarthLocation(-90*u.deg, 65*u.deg),
                       obstime=Time('2014-08-01 08:00:00')),
                 AltAz(location=EarthLocation(120*u.deg, -35*u.deg),
                       obstime=Time('2014-01-01 00:00:00'))
                ]
MOONDIST = 385000*u.km  # approximate moon semi-major orbit axis of moon
MOONDIST_CART = CartesianRepresentation(3**-0.5*MOONDIST, 3**-0.5*MOONDIST, 3**-0.5*MOONDIST)
EARTHECC = 0.017 + 0.005  # roughly earth orbital eccentricity, but with an added tolerance


@pytest.mark.parametrize('testframe', totest_frames)
def test_gcrs_altaz_sunish(testframe):
    """
    Sanity-check that the sun is at a reasonable distance from any altaz
    """
    sun = get_sun(testframe.obstime)

    assert sun.frame.name == 'gcrs'

    # the .to(u.au) is not necessary, it just makes the asserts on failure more readable
    assert (EARTHECC - 1)*u.au < sun.distance.to(u.au) < (EARTHECC + 1)*u.au

    sunaa = sun.transform_to(testframe)
    assert (EARTHECC - 1)*u.au < sunaa.distance.to(u.au) < (EARTHECC + 1)*u.au


@pytest.mark.parametrize('testframe', totest_frames)
def test_gcrs_altaz_moonish(testframe):
    """
    Sanity-check that an object resembling the moon goes to the right place with
    a GCRS->AltAz transformation
    """
    moon = GCRS(MOONDIST_CART, obstime=testframe.obstime)

    moonaa = moon.transform_to(testframe)

    # now check that the distance change is similar to earth radius
    assert 1000*u.km < np.abs(moonaa.distance - moon.distance).to(u.au) < 7000*u.km

    # now check that it round-trips
    moon2 = moonaa.transform_to(moon)
    assert_allclose(moon.cartesian.xyz, moon2.cartesian.xyz)

    # also should add checks that the alt/az are different for different earth locations


@pytest.mark.parametrize('testframe', totest_frames)
def test_gcrs_altaz_bothroutes(testframe):
    """
    Repeat of both the moonish and sunish tests above to make sure the two
    routes through the coordinate graph are consistent with each other
    """
    sun = get_sun(testframe.obstime)
    sunaa_viaicrs = sun.transform_to(ICRS).transform_to(testframe)
    sunaa_viaitrs = sun.transform_to(ITRS(obstime=testframe.obstime)).transform_to(testframe)

    moon = GCRS(MOONDIST_CART, obstime=testframe.obstime)
    moonaa_viaicrs = moon.transform_to(ICRS).transform_to(testframe)
    moonaa_viaitrs = moon.transform_to(ITRS(obstime=testframe.obstime)).transform_to(testframe)

    assert_allclose(sunaa_viaicrs.cartesian.xyz, sunaa_viaitrs.cartesian.xyz)
    assert_allclose(moonaa_viaicrs.cartesian.xyz, moonaa_viaitrs.cartesian.xyz)


@pytest.mark.parametrize('testframe', totest_frames)
def test_cirs_altaz_moonish(testframe):
    """
    Sanity-check that an object resembling the moon goes to the right place with
    a CIRS<->AltAz transformation
    """
    moon = CIRS(MOONDIST_CART, obstime=testframe.obstime)

    moonaa = moon.transform_to(testframe)
    assert 1000*u.km < np.abs(moonaa.distance - moon.distance).to(u.km) < 7000*u.km

    # now check that it round-trips
    moon2 = moonaa.transform_to(moon)
    assert_allclose(moon.cartesian.xyz, moon2.cartesian.xyz)


@pytest.mark.parametrize('testframe', totest_frames)
def test_cirs_altaz_nodist(testframe):
    """
    Check that a UnitSphericalRepresentation coordinate round-trips for the
    CIRS<->AltAz transformation.
    """
    coo0 = CIRS(UnitSphericalRepresentation(10*u.deg, 20*u.deg), obstime=testframe.obstime)

    # check that it round-trips
    coo1 = coo0.transform_to(testframe).transform_to(coo0)
    assert_allclose(coo0.cartesian.xyz, coo1.cartesian.xyz)


@pytest.mark.parametrize('testframe', totest_frames)
def test_cirs_icrs_moonish(testframe):
    """
    check that something like the moon goes to about the right distance from the
    ICRS origin when starting from CIRS
    """
    moonish = CIRS(MOONDIST_CART, obstime=testframe.obstime)
    moonicrs = moonish.transform_to(ICRS)

    assert 0.97*u.au < moonicrs.distance < 1.03*u.au


@pytest.mark.parametrize('testframe', totest_frames)
def test_gcrs_icrs_moonish(testframe):
    """
    check that something like the moon goes to about the right distance from the
    ICRS origin when starting from GCRS
    """
    moonish = GCRS(MOONDIST_CART, obstime=testframe.obstime)
    moonicrs = moonish.transform_to(ICRS)

    assert 0.97*u.au < moonicrs.distance < 1.03*u.au


@pytest.mark.parametrize('testframe', totest_frames)
def test_icrs_gcrscirs_sunish(testframe):
    """
    check that the ICRS barycenter goes to about the right distance from various
    ~geocentric frames (other than testframe)
    """
    # slight offset to avoid divide-by-zero errors
    icrs = ICRS(0*u.deg, 0*u.deg, distance=10*u.km)

    gcrs = icrs.transform_to(GCRS(obstime=testframe.obstime))
    assert (EARTHECC - 1)*u.au < gcrs.distance.to(u.au) < (EARTHECC + 1)*u.au

    cirs = icrs.transform_to(CIRS(obstime=testframe.obstime))
    assert (EARTHECC - 1)*u.au < cirs.distance.to(u.au) < (EARTHECC + 1)*u.au

    itrs = icrs.transform_to(ITRS(obstime=testframe.obstime))
    assert (EARTHECC - 1)*u.au < itrs.spherical.distance.to(u.au) < (EARTHECC + 1)*u.au


@pytest.mark.parametrize('testframe', totest_frames)
def test_icrs_altaz_moonish(testframe):
    """
    Check that something expressed in *ICRS* as being moon-like goes to the
    right AltAz distance
    """
    # we use epv00 instead of get_sun because get_sun includes aberration
    earth_pv_helio, earth_pv_bary = epv00(*get_jd12(testframe.obstime, 'tdb'))
    earth_icrs_xyz = earth_pv_bary[0]*u.au
    moonoffset = [0, 0, MOONDIST.value]*MOONDIST.unit
    moonish_icrs = ICRS(CartesianRepresentation(earth_icrs_xyz + moonoffset))
    moonaa = moonish_icrs.transform_to(testframe)

    # now check that the distance change is similar to earth radius
    assert 1000*u.km < np.abs(moonaa.distance - MOONDIST).to(u.au) < 7000*u.km


def test_gcrs_self_transform_closeby():
    """
    Tests GCRS self transform for objects which are nearby and thus
    have reasonable parallax.

    Moon positions were originally created using JPL DE432s ephemeris.

    The two lunar positions (one geocentric, one at a defined location)
    are created via a transformation from ICRS to two different GCRS frames.

    We test that the GCRS-GCRS self transform can correctly map one GCRS
    frame onto the other.
    """
    t = Time("2014-12-25T07:00")
    moon_geocentric = SkyCoord(GCRS(318.10579159*u.deg,
                                    -11.65281165*u.deg,
                                    365042.64880308*u.km, obstime=t))

    # this is the location of the Moon as seen from La Palma
    obsgeoloc = [-5592982.59658935, -63054.1948592, 3059763.90102216]*u.m
    obsgeovel = [4.59798494, -407.84677071, 0.]*u.m/u.s
    moon_lapalma = SkyCoord(GCRS(318.7048445*u.deg,
                                 -11.98761996*u.deg,
                                 369722.8231031*u.km,
                                 obstime=t,
                                 obsgeoloc=obsgeoloc,
                                 obsgeovel=obsgeovel))

    transformed = moon_geocentric.transform_to(moon_lapalma.frame)
    delta = transformed.separation_3d(moon_lapalma)
    assert_allclose(delta, 0.0*u.m, atol=1*u.m)


@pytest.mark.remote_data
@pytest.mark.skipif('not HAS_JPLEPHEM')
def test_ephemerides():
    """
    We test that using different ephemerides gives very similar results
    for transformations
    """
    t = Time("2014-12-25T07:00")
    moon = SkyCoord(GCRS(318.10579159*u.deg,
                         -11.65281165*u.deg,
                         365042.64880308*u.km, obstime=t))

    icrs_frame = ICRS()
    hcrs_frame = HCRS(obstime=t)
    ecl_frame = HeliocentricTrueEcliptic(equinox=t)
    cirs_frame = CIRS(obstime=t)

    moon_icrs_builtin = moon.transform_to(icrs_frame)
    moon_hcrs_builtin = moon.transform_to(hcrs_frame)
    moon_helioecl_builtin = moon.transform_to(ecl_frame)
    moon_cirs_builtin = moon.transform_to(cirs_frame)

    with solar_system_ephemeris.set('jpl'):
        moon_icrs_jpl = moon.transform_to(icrs_frame)
        moon_hcrs_jpl = moon.transform_to(hcrs_frame)
        moon_helioecl_jpl = moon.transform_to(ecl_frame)
        moon_cirs_jpl = moon.transform_to(cirs_frame)

    # most transformations should differ by an amount which is
    # non-zero but of order milliarcsecs
    sep_icrs = moon_icrs_builtin.separation(moon_icrs_jpl)
    sep_hcrs = moon_hcrs_builtin.separation(moon_hcrs_jpl)
    sep_helioecl = moon_helioecl_builtin.separation(moon_helioecl_jpl)
    sep_cirs = moon_cirs_builtin.separation(moon_cirs_jpl)

    assert_allclose([sep_icrs, sep_hcrs, sep_helioecl], 0.0*u.deg, atol=10*u.mas)
    assert all(sep > 10*u.microarcsecond for sep in (sep_icrs, sep_hcrs, sep_helioecl))

    # CIRS should be the same
    assert_allclose(sep_cirs, 0.0*u.deg, atol=1*u.microarcsecond)

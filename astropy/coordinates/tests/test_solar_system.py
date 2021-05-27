import os

import pytest
import numpy as np
from urllib.error import HTTPError

from astropy.time import Time
from astropy import units as u
from astropy.constants import c
from astropy.coordinates.builtin_frames import GCRS, TETE
from astropy.coordinates.earth import EarthLocation
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.coordinates.representation import CartesianRepresentation, UnitSphericalRepresentation
from astropy.coordinates.solar_system import (get_body, get_moon, BODY_NAME_TO_KERNEL_SPEC,
                                              _get_apparent_body_position, solar_system_ephemeris,
                                              get_body_barycentric, get_body_barycentric_posvel)
from astropy.coordinates.funcs import get_sun
from astropy.tests.helper import assert_quantity_allclose
from astropy.units import allclose as quantity_allclose
from astropy.utils.data import download_file
from astropy.utils.compat.optional_deps import (HAS_JPLEPHEM,  # noqa
                                                HAS_SKYFIELD)

if HAS_SKYFIELD:
    from skyfield.api import Loader, Topos

de432s_separation_tolerance_planets = 5*u.arcsec
de432s_separation_tolerance_moon = 5*u.arcsec
de432s_distance_tolerance = 20*u.km

skyfield_angular_separation_tolerance = 1*u.arcsec
skyfield_separation_tolerance = 10*u.km


@pytest.mark.remote_data
@pytest.mark.skipif('not HAS_SKYFIELD')
def test_positions_skyfield(tmpdir):
    """
    Test positions against those generated by skyfield.
    """
    load = Loader(tmpdir)
    t = Time('1980-03-25 00:00')
    location = None

    # skyfield ephemeris
    try:
        planets = load('de421.bsp')
        ts = load.timescale()
    except OSError as e:
        if os.environ.get('CI', False) and 'timed out' in str(e):
            pytest.xfail('Timed out in CI')
        else:
            raise

    mercury, jupiter, moon = planets['mercury'], planets['jupiter barycenter'], planets['moon']
    earth = planets['earth']

    skyfield_t = ts.from_astropy(t)

    if location is not None:
        earth = earth+Topos(latitude_degrees=location.lat.to_value(u.deg),
                            longitude_degrees=location.lon.to_value(u.deg),
                            elevation_m=location.height.to_value(u.m))

    skyfield_mercury = earth.at(skyfield_t).observe(mercury).apparent()
    skyfield_jupiter = earth.at(skyfield_t).observe(jupiter).apparent()
    skyfield_moon = earth.at(skyfield_t).observe(moon).apparent()

    if location is not None:
        frame = TETE(obstime=t, location=location)
    else:
        frame = TETE(obstime=t)

    ra, dec, dist = skyfield_mercury.radec(epoch='date')
    skyfield_mercury = SkyCoord(ra.to(u.deg), dec.to(u.deg), distance=dist.to(u.km),
                                frame=frame)
    ra, dec, dist = skyfield_jupiter.radec(epoch='date')
    skyfield_jupiter = SkyCoord(ra.to(u.deg), dec.to(u.deg), distance=dist.to(u.km),
                                frame=frame)
    ra, dec, dist = skyfield_moon.radec(epoch='date')
    skyfield_moon = SkyCoord(ra.to(u.deg), dec.to(u.deg), distance=dist.to(u.km),
                             frame=frame)

    # planet positions w.r.t true equator and equinox
    moon_astropy = get_moon(t, location, ephemeris='de430').transform_to(frame)
    mercury_astropy = get_body('mercury', t, location, ephemeris='de430').transform_to(frame)
    jupiter_astropy = get_body('jupiter', t, location, ephemeris='de430').transform_to(frame)

    assert (moon_astropy.separation(skyfield_moon) <
            skyfield_angular_separation_tolerance)
    assert (moon_astropy.separation_3d(skyfield_moon) < skyfield_separation_tolerance)

    assert (jupiter_astropy.separation(skyfield_jupiter) <
            skyfield_angular_separation_tolerance)
    assert (jupiter_astropy.separation_3d(skyfield_jupiter) <
            skyfield_separation_tolerance)

    assert (mercury_astropy.separation(skyfield_mercury) <
            skyfield_angular_separation_tolerance)
    assert (mercury_astropy.separation_3d(skyfield_mercury) <
            skyfield_separation_tolerance)


class TestPositionsGeocentric:
    """
    Test positions against those generated by JPL Horizons accessed on
    2016-03-28, with refraction turned on.
    """

    def setup(self):
        self.t = Time('1980-03-25 00:00')
        self.apparent_frame = TETE(obstime=self.t)
        # Results returned by JPL Horizons web interface
        self.horizons = {
            'mercury': SkyCoord(ra='22h41m47.78s', dec='-08d29m32.0s',
                                distance=c*6.323037*u.min, frame=self.apparent_frame),
            'moon': SkyCoord(ra='07h32m02.62s', dec='+18d34m05.0s',
                             distance=c*0.021921*u.min, frame=self.apparent_frame),
            'jupiter': SkyCoord(ra='10h17m12.82s', dec='+12d02m57.0s',
                                distance=c*37.694557*u.min, frame=self.apparent_frame),
            'sun': SkyCoord(ra='00h16m31.00s', dec='+01d47m16.9s',
                            distance=c*8.294858*u.min, frame=self.apparent_frame)}

    @pytest.mark.parametrize(('body', 'sep_tol', 'dist_tol'),
                             (('mercury', 7.*u.arcsec, 1000*u.km),
                              ('jupiter', 78.*u.arcsec, 76000*u.km),
                              ('moon', 20.*u.arcsec, 80*u.km),
                              ('sun', 5.*u.arcsec, 11.*u.km)))
    def test_erfa_planet(self, body, sep_tol, dist_tol):
        """Test predictions using erfa/plan94.

        Accuracies are maximum deviations listed in erfa/plan94.c, for Jupiter and
        Mercury, and that quoted in Meeus "Astronomical Algorithms" (1998) for the Moon.
        """
        astropy = get_body(body, self.t, ephemeris='builtin')
        horizons = self.horizons[body]

        # convert to true equator and equinox
        astropy = astropy.transform_to(self.apparent_frame)

        # Assert sky coordinates are close.
        assert astropy.separation(horizons) < sep_tol

        # Assert distances are close.
        assert_quantity_allclose(astropy.distance, horizons.distance,
                                 atol=dist_tol)

    @pytest.mark.remote_data
    @pytest.mark.skipif('not HAS_JPLEPHEM')
    @pytest.mark.parametrize('body', ('mercury', 'jupiter', 'sun'))
    def test_de432s_planet(self, body):
        astropy = get_body(body, self.t, ephemeris='de432s')
        horizons = self.horizons[body]

        # convert to true equator and equinox
        astropy = astropy.transform_to(self.apparent_frame)

        # Assert sky coordinates are close.
        assert (astropy.separation(horizons) <
                de432s_separation_tolerance_planets)

        # Assert distances are close.
        assert_quantity_allclose(astropy.distance, horizons.distance,
                                 atol=de432s_distance_tolerance)

    @pytest.mark.remote_data
    @pytest.mark.skipif('not HAS_JPLEPHEM')
    def test_de432s_moon(self):
        astropy = get_moon(self.t, ephemeris='de432s')
        horizons = self.horizons['moon']

        # convert to true equator and equinox
        astropy = astropy.transform_to(self.apparent_frame)

        # Assert sky coordinates are close.
        assert (astropy.separation(horizons) <
                de432s_separation_tolerance_moon)

        # Assert distances are close.
        assert_quantity_allclose(astropy.distance, horizons.distance,
                                 atol=de432s_distance_tolerance)


class TestPositionKittPeak:
    """
    Test positions against those generated by JPL Horizons accessed on
    2016-03-28, with refraction turned on.
    """

    def setup(self):
        kitt_peak = EarthLocation.from_geodetic(lon=-111.6*u.deg,
                                                lat=31.963333333333342*u.deg,
                                                height=2120*u.m)
        self.t = Time('2014-09-25T00:00', location=kitt_peak)
        self.apparent_frame = TETE(obstime=self.t, location=kitt_peak)
        # Results returned by JPL Horizons web interface
        self.horizons = {
            'mercury': SkyCoord(ra='13h38m58.50s', dec='-13d34m42.6s',
                                distance=c*7.699020*u.min, frame=self.apparent_frame),
            'moon': SkyCoord(ra='12h33m12.85s', dec='-05d17m54.4s',
                             distance=c*0.022054*u.min, frame=self.apparent_frame),
            'jupiter': SkyCoord(ra='09h09m55.55s', dec='+16d51m57.8s',
                                distance=c*49.244937*u.min, frame=self.apparent_frame)}

    @pytest.mark.parametrize(('body', 'sep_tol', 'dist_tol'),
                             (('mercury', 7.*u.arcsec, 500*u.km),
                              ('jupiter', 78.*u.arcsec, 82000*u.km)))
    def test_erfa_planet(self, body, sep_tol, dist_tol):
        """Test predictions using erfa/plan94.

        Accuracies are maximum deviations listed in erfa/plan94.c.
        """
        # Add uncertainty in position of Earth
        dist_tol = dist_tol + 1300 * u.km

        astropy = get_body(body, self.t, ephemeris='builtin')
        horizons = self.horizons[body]

        # convert to true equator and equinox
        astropy = astropy.transform_to(self.apparent_frame)

        # Assert sky coordinates are close.
        assert astropy.separation(horizons) < sep_tol

        # Assert distances are close.
        assert_quantity_allclose(astropy.distance, horizons.distance,
                                 atol=dist_tol)

    @pytest.mark.remote_data
    @pytest.mark.skipif('not HAS_JPLEPHEM')
    @pytest.mark.parametrize('body', ('mercury', 'jupiter'))
    def test_de432s_planet(self, body):
        astropy = get_body(body, self.t, ephemeris='de432s')
        horizons = self.horizons[body]

        # convert to true equator and equinox
        astropy = astropy.transform_to(self.apparent_frame)

        # Assert sky coordinates are close.
        assert (astropy.separation(horizons) <
                de432s_separation_tolerance_planets)

        # Assert distances are close.
        assert_quantity_allclose(astropy.distance, horizons.distance,
                                 atol=de432s_distance_tolerance)

    @pytest.mark.remote_data
    @pytest.mark.skipif('not HAS_JPLEPHEM')
    def test_de432s_moon(self):
        astropy = get_moon(self.t, ephemeris='de432s')
        horizons = self.horizons['moon']

        # convert to true equator and equinox
        astropy = astropy.transform_to(self.apparent_frame)

        # Assert sky coordinates are close.
        assert (astropy.separation(horizons) <
                de432s_separation_tolerance_moon)

        # Assert distances are close.
        assert_quantity_allclose(astropy.distance, horizons.distance,
                                 atol=de432s_distance_tolerance)

    @pytest.mark.remote_data
    @pytest.mark.skipif('not HAS_JPLEPHEM')
    @pytest.mark.parametrize('bodyname', ('mercury', 'jupiter'))
    def test_custom_kernel_spec_body(self, bodyname):
        """
        Checks that giving a kernel specifier instead of a body name works
        """
        coord_by_name = get_body(bodyname, self.t, ephemeris='de432s')
        kspec = BODY_NAME_TO_KERNEL_SPEC[bodyname]
        coord_by_kspec = get_body(kspec, self.t, ephemeris='de432s')

        assert_quantity_allclose(coord_by_name.ra, coord_by_kspec.ra)
        assert_quantity_allclose(coord_by_name.dec, coord_by_kspec.dec)
        assert_quantity_allclose(coord_by_name.distance, coord_by_kspec.distance)


@pytest.mark.remote_data
def test_horizons_consistency_with_precision():
    """
    A test to compare at high precision against output of JPL horizons.

    Tests ephemerides, and conversions from ICRS to GCRS to TETE. We are aiming for
    better than 2 milli-arcsecond precision.

    We use the Moon since it is nearby, and moves fast in the sky so we are
    testing for parallax, proper handling of light deflection and aberration.
    """
    # JPL Horizon values for 2020_04_06 00:00 to 23:00 in 1 hour steps
    # JPL Horizons has a known offset (frame bias) of 51.02 mas in RA. We correct that here
    ra_apparent_horizons = [
        170.167332531, 170.560688674, 170.923834838, 171.271663481, 171.620188972, 171.985340827,
        172.381766539, 172.821772139, 173.314502650, 173.865422398, 174.476108551, 175.144332386,
        175.864375310, 176.627519827, 177.422655853, 178.236955730, 179.056584831, 179.867427392,
        180.655815385, 181.409252074, 182.117113814, 182.771311578, 183.366872837, 183.902395443
    ] * u.deg + 51.02376467 * u.mas
    dec_apparent_horizons = [
        10.269112037, 10.058820647, 9.837152044, 9.603724551, 9.358956528, 9.104012390, 8.840674927,
        8.571162442, 8.297917326, 8.023394488, 7.749873882, 7.479312991, 7.213246666, 6.952732614,
        6.698336823, 6.450150213, 6.207828142, 5.970645962, 5.737565957, 5.507313851, 5.278462034,
        5.049521497, 4.819038911, 4.585696512
    ] * u.deg
    with solar_system_ephemeris.set('de430'):
        loc = EarthLocation.from_geodetic(-67.787260*u.deg, -22.959748*u.deg, 5186*u.m)
        times = Time('2020-04-06 00:00') + np.arange(0, 24, 1)*u.hour
        astropy = get_body('moon', times, loc)

        apparent_frame = TETE(obstime=times, location=loc)
        astropy = astropy.transform_to(apparent_frame)
        usrepr = UnitSphericalRepresentation(ra_apparent_horizons, dec_apparent_horizons)
        horizons = apparent_frame.realize_frame(usrepr)
    assert_quantity_allclose(astropy.separation(horizons), 0*u.mas, atol=1.5*u.mas)


@pytest.mark.remote_data
@pytest.mark.skipif('not HAS_JPLEPHEM')
@pytest.mark.parametrize('time', (Time('1960-01-12 00:00'),
                                  Time('1980-03-25 00:00'),
                                  Time('2010-10-13 00:00')))
def test_get_sun_consistency(time):
    """
    Test that the sun from JPL and the builtin get_sun match
    """
    sun_jpl_gcrs = get_body('sun', time, ephemeris='de432s')
    builtin_get_sun = get_sun(time)
    sep = builtin_get_sun.separation(sun_jpl_gcrs)
    assert sep < 0.1*u.arcsec


def test_get_moon_nonscalar_regression():
    """
    Test that the builtin ephemeris works with non-scalar times.

    See Issue #5069.
    """
    times = Time(["2015-08-28 03:30", "2015-09-05 10:30"])
    # the following line will raise an Exception if the bug recurs.
    get_moon(times, ephemeris='builtin')


def test_barycentric_pos_posvel_same():
    # Check that the two routines give identical results.
    ep1 = get_body_barycentric('earth', Time('2016-03-20T12:30:00'))
    ep2, _ = get_body_barycentric_posvel('earth', Time('2016-03-20T12:30:00'))
    assert np.all(ep1.xyz == ep2.xyz)


def test_earth_barycentric_velocity_rough():
    # Check that a time near the equinox gives roughly the right result.
    ep, ev = get_body_barycentric_posvel('earth', Time('2016-03-20T12:30:00'))
    assert_quantity_allclose(ep.xyz, [-1., 0., 0.]*u.AU, atol=0.01*u.AU)
    expected = u.Quantity([0.*u.one,
                           np.cos(23.5*u.deg),
                           np.sin(23.5*u.deg)]) * -30. * u.km / u.s
    assert_quantity_allclose(ev.xyz, expected, atol=1.*u.km/u.s)


def test_earth_barycentric_velocity_multi_d():
    # Might as well test it with a multidimensional array too.
    t = Time('2016-03-20T12:30:00') + np.arange(8.).reshape(2, 2, 2) * u.yr / 2.
    ep, ev = get_body_barycentric_posvel('earth', t)
    # note: assert_quantity_allclose doesn't like the shape mismatch.
    # this is a problem with np.testing.assert_allclose.
    assert quantity_allclose(ep.get_xyz(xyz_axis=-1),
                             [[-1., 0., 0.], [+1., 0., 0.]]*u.AU,
                             atol=0.06*u.AU)
    expected = u.Quantity([0.*u.one,
                           np.cos(23.5*u.deg),
                           np.sin(23.5*u.deg)]) * ([[-30.], [30.]] * u.km / u.s)
    assert quantity_allclose(ev.get_xyz(xyz_axis=-1), expected,
                             atol=2.*u.km/u.s)


@pytest.mark.remote_data
@pytest.mark.skipif('not HAS_JPLEPHEM')
@pytest.mark.parametrize(('body', 'pos_tol', 'vel_tol'),
                         (('mercury', 1000.*u.km, 1.*u.km/u.s),
                          ('jupiter', 100000.*u.km, 2.*u.km/u.s),
                          ('earth', 10*u.km, 10*u.mm/u.s),
                          ('moon', 18*u.km, 50*u.mm/u.s)))
def test_barycentric_velocity_consistency(body, pos_tol, vel_tol):
    # Tolerances are about 1.5 times the rms listed for plan94 and epv00,
    # except for Mercury (which nominally is 334 km rms), and the Moon
    # (which nominally is 6 km rms).
    t = Time('2016-03-20T12:30:00')
    ep, ev = get_body_barycentric_posvel(body, t, ephemeris='builtin')
    dp, dv = get_body_barycentric_posvel(body, t, ephemeris='de432s')
    assert_quantity_allclose(ep.xyz, dp.xyz, atol=pos_tol)
    assert_quantity_allclose(ev.xyz, dv.xyz, atol=vel_tol)
    # Might as well test it with a multidimensional array too.
    t = Time('2016-03-20T12:30:00') + np.arange(8.).reshape(2, 2, 2) * u.yr / 2.
    ep, ev = get_body_barycentric_posvel(body, t, ephemeris='builtin')
    dp, dv = get_body_barycentric_posvel(body, t, ephemeris='de432s')
    assert_quantity_allclose(ep.xyz, dp.xyz, atol=pos_tol)
    assert_quantity_allclose(ev.xyz, dv.xyz, atol=vel_tol)


@pytest.mark.remote_data
@pytest.mark.skipif('not HAS_JPLEPHEM')
@pytest.mark.parametrize('time', (Time('1960-01-12 00:00'),
                                  Time('1980-03-25 00:00'),
                                  Time('2010-10-13 00:00')))
def test_url_or_file_ephemeris(time):
    # URL for ephemeris de432s used for testing:
    url = 'http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp'

    # Pass the ephemeris directly as a URL.
    coord_by_url = get_body('earth', time, ephemeris=url)

    # Translate the URL to the cached location on the filesystem.
    # Since we just used the url above, it should already have been downloaded.
    filepath = download_file(url, cache=True)

    # Get the coordinates using the file path directly:
    coord_by_filepath = get_body('earth', time, ephemeris=filepath)

    # Using the URL or filepath should give exactly the same results:
    assert_quantity_allclose(coord_by_url.ra, coord_by_filepath.ra)
    assert_quantity_allclose(coord_by_url.dec, coord_by_filepath.dec)
    assert_quantity_allclose(coord_by_url.distance, coord_by_filepath.distance)


@pytest.mark.remote_data
@pytest.mark.skipif('not HAS_JPLEPHEM')
def test_url_ephemeris_wrong_input():
    # Try loading a non-existing URL:
    time = Time('1960-01-12 00:00')
    with pytest.raises(HTTPError):
        get_body('earth', time, ephemeris='http://data.astropy.org/path/to/nonexisting/file.bsp')


@pytest.mark.skipif('not HAS_JPLEPHEM')
def test_file_ephemeris_wrong_input():
    time = Time('1960-01-12 00:00')
    # Try loading a non-existing file:
    with pytest.raises(ValueError):
        get_body('earth', time, ephemeris='/path/to/nonexisting/file.bsp')

    # This test currently leaves the file open. To fix this issue, an
    # upstream fix is required in jplephem package. This test is removed
    # until that issue is fixed.
    # Try loading a file that does exist, but is not an ephemeris file:
    # with pytest.raises(ValueError):
    #     get_body('earth', time, ephemeris=__file__)


def test_regression_10271():
    t = Time(58973.534052125986, format='mjd')
    # GCRS position of ALMA at this time
    obs_p = CartesianRepresentation(5724535.74068625, -1311071.58985697, -2492738.93017009, u.m)

    geocentre = CartesianRepresentation(0, 0, 0, u.m)
    icrs_sun_from_alma = _get_apparent_body_position('sun', t, 'builtin', obs_p)
    icrs_sun_from_geocentre = _get_apparent_body_position('sun', t, 'builtin', geocentre)

    difference = (icrs_sun_from_alma - icrs_sun_from_geocentre).norm()
    assert_quantity_allclose(difference, 0.13046941*u.m, atol=1*u.mm)

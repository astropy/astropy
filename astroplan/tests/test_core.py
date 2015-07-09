from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, Latitude, Longitude, SkyCoord)
import astropy.units as u
from astropy.time import Time
from astropy.tests.helper import remote_data
import numpy as np
from numpy.testing import assert_allclose
import pytz
import datetime
import unittest

from ..core import FixedTarget, Observer
from ..exceptions import TargetAlwaysUpWarning, TargetNeverUpWarning

def test_Observer_constructor_location():
    '''
    Show that location defined by latitude/longitude/elevation is parsed
    identically to passing in an `~astropy.coordinates.EarthLocation` directly.
    '''

    lat = '+19:00:00'
    lon = '-155:00:00'
    elevation = 0.0 * u.m
    location = EarthLocation.from_geodetic(lon, lat, elevation)

    environment_kwargs = dict(pressure=1*u.bar, relative_humidity=0.1,
                              temperature=10*u.deg_C)

    obs1 = Observer(name='Observatory',
                    latitude=lat,
                    longitude=lon,
                    elevation=elevation,
                    **environment_kwargs)

    obs2 = Observer(name='Observatory',
                    location=location,
                    **environment_kwargs)

    assert obs1.location == obs2.location, ('using latitude/longitude/'
                                            'elevation keywords gave a '
                                            'different answer from passing in '
                                            'an EarthLocation directly')

@remote_data
def test_FixedTarget_from_name():
    '''
    Check that resolving target names with the `SkyCoord.from_name` constructor
    to produce a `FixedTarget` accurately resolves the coordinates of Polaris.
    '''

    # Resolve coordinates with SkyCoord.from_name classmethod
    polaris_from_name = FixedTarget.from_name('Polaris')
    polaris_from_name = FixedTarget.from_name('Polaris', name='Target 1')
    # Coordinates grabbed from SIMBAD
    polaris_from_SIMBAD = SkyCoord('02h31m49.09456s', '+89d15m50.7923s')

    # Make sure separation is small
    assert polaris_from_name.coord.separation(polaris_from_SIMBAD) < 1*u.arcsec

def test_Observer_altaz():
    '''
    Check that the altitude/azimuth computed by `Observer.altaz` is similar
    to the result from PyEphem when pressure = 0 (no atmosphere) for Vega at
    2000-01-01 12:00:00 UTC.
    '''
    # Define the test case
    latitude = '00:00:00'
    longitude = '00:00:00'
    elevation = 0 # [m]
    pressure = 0  * u.bar # no atmosphere
    time = Time('2000-01-01 12:00:00')
    vega_coords = SkyCoord('18h36m56.33635s', '+38d47m01.2802s')

    # Calculate altitude/azimuth with astroplan
    location = EarthLocation.from_geodetic(longitude, latitude, elevation*u.m)
    astroplan_obs = Observer(name='Observatory', location=location,
                             pressure=pressure*u.bar)
    astroplan_vega = FixedTarget(vega_coords)
    altaz = astroplan_obs.altaz(time, astroplan_vega)
    astroplan_altitude = altaz.alt
    astroplan_azimuth = altaz.az

    # Calculate altitude/azimuth with PyEphem, like so with print_pyephem_altaz
    pyephem_altitude = Latitude('51.198848716510874 deg')
    pyephem_azimuth = Longitude('358.4676707379987 deg')

    # Assert that altitudes/azimuths are within 30 arcsec - this is a wide
    # tolerance because the IERS tables used by astroplan may offset astroplan's
    # positions due to leap seconds.
    tolerance = (30*u.arcsec).to('deg').value
    assert_allclose(pyephem_altitude.value, astroplan_altitude.value,
                    atol=tolerance)
    assert_allclose(pyephem_azimuth.value, astroplan_azimuth.value,
                    atol=tolerance)

    # Check that alt/az without target returns AltAz frame
    from astropy.coordinates import AltAz
    assert isinstance(astroplan_obs.altaz(time), AltAz)

def print_pyephem_altaz(latitude, longitude, elevation, time, pressure,
                      target_coords):
    '''
    Run PyEphem to compute the altitude/azimuth of a target at specified time
    and observatory, for comparison with astroplan calucation tested in
    `test_Observer_altaz`.
    '''
    import ephem
    pyephem_obs = ephem.Observer()
    pyephem_obs.lat = latitude
    pyephem_obs.lon = longitude
    pyephem_obs.elevation = elevation
    pyephem_obs.date = time.datetime
    pyephem_obs.pressure = pressure
    pyephem_target = ephem.FixedBody()
    pyephem_target._ra = ephem.degrees(np.radians(target_coords.ra.value))
    pyephem_target._dec = ephem.degrees(np.radians(target_coords.dec.value))
    pyephem_target.compute(pyephem_obs)
    pyephem_altitude = Latitude(np.degrees(pyephem_target.alt)*u.degree)
    pyephem_azimuth = Longitude(np.degrees(pyephem_target.az)*u.degree)
    print(pyephem_altitude, pyephem_azimuth)

def test_Observer_timezone_parser():
    lat = '+19:00:00'
    lon = '-155:00:00'
    elevation = 0.0 * u.m
    location = EarthLocation.from_geodetic(lon, lat, elevation)

    obs1 = Observer(name='Observatory', location=location,
                    timezone=pytz.timezone('UTC'))
    obs2 = Observer(name='Observatory', location=location, timezone='UTC')
    obs3 = Observer(name='Observatory', location=location)

    assert obs1.timezone == obs2.timezone, ('Accept both strings to pass to '
                                            'the pytz.timezone() constructor '
                                            'and instances of pytz.timezone')

    assert obs2.timezone == obs3.timezone, ('Default timezone should be UTC')

def test_FixedTarget_ra_dec():
    '''
    Confirm that FixedTarget.ra and FixedTarget.dec are the same as the
    right ascension and declination stored in the FixedTarget.coord variable -
    which is a SkyCoord
    '''

    vega_coords = SkyCoord('18h36m56.33635s', '+38d47m01.2802s')
    vega = FixedTarget(vega_coords, name='Vega')
    assert vega.coord == vega_coords, 'Store coordinates directly'
    assert vega.coord.ra == vega_coords.ra == vega.ra, ('Retrieve RA from '
                                                        'SkyCoord')
    assert vega.coord.dec == vega_coords.dec == vega.dec, ('Retrieve Dec from '
                                                           'SkyCoord')

def test_parallactic_angle():
    '''
    Compute parallactic angle for targets at hour angle = {3, 19} for
    at observer at IRTF using the online SpeX calculator and PyEphem
    '''
    # Set up position for IRTF
    lat = 19.826218*u.deg
    lon = -155.471999*u.deg
    elevation = 4160.0 * u.m
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2015-01-01 00:00:00')
    LST = time.sidereal_time('mean', longitude=lon)
    desired_HA_1 = 3*u.hourangle
    desired_HA_2 = 19*u.hourangle # = -5*u.hourangle
    target1 = SkyCoord(LST - desired_HA_1, -30*u.degree)
    target2 = SkyCoord(LST - desired_HA_2, -30*u.degree)

    obs = Observer(location=location)
    q1 = obs.parallactic_angle(time, target1)
    q2 = obs.parallactic_angle(time, target2)

    # Get values from PyEphem for comparison from print_pyephem_parallactic_angle()
    pyephem_q1 = 46.54610060782033*u.deg
    pyephem_q2 = -65.51818282032019*u.deg

    assert_allclose(q1.to(u.degree).value, pyephem_q1, atol=1)
    assert_allclose(q2.to(u.degree).value, pyephem_q2, atol=1)

    # Get SpeX parallactic angle calculator values for comparison from
    # http://irtfweb.ifa.hawaii.edu/cgi-bin/spex/parangle.cgi to produce

    SpeX_q1 = 46.7237968 # deg
    SpeX_q2 = -65.428924 # deg

    assert_allclose(q1.to(u.degree).value, SpeX_q1, atol=0.1)
    assert_allclose(q2.to(u.degree).value, SpeX_q2, atol=0.1)

def print_pyephem_parallactic_angle():
    lat = 19.826218*u.deg
    lon = -155.471999*u.deg
    time = Time('2015-01-01 00:00:00')
    LST = time.sidereal_time('mean', longitude=lon)
    desired_HA_1 = 3*u.hourangle
    desired_HA_2 = 19*u.hourangle # = -5*u.hourangle

    import ephem
    obs = ephem.Observer()
    obs.lat = '19:49:34.3848'
    obs.lon = '-155:28:19.1964'
    obs.date = time.datetime
    pyephem_target1 = ephem.FixedBody()
    pyephem_target1._ra = ephem.degrees((LST - desired_HA_1).to(u.rad).value)
    pyephem_target1._dec = ephem.degrees((-30*u.deg).to(u.rad).value)
    pyephem_target1.compute(obs)
    pyephem_q1 = (float(pyephem_target1.parallactic_angle())*u.rad).to(u.deg)

    pyephem_target2 = ephem.FixedBody()
    pyephem_target2._ra = ephem.degrees((LST - desired_HA_2).to(u.rad).value)
    pyephem_target2._dec = ephem.degrees((-30*u.deg).to(u.rad).value)
    pyephem_target2.compute(obs)
    pyephem_q2 = (float(pyephem_target2.parallactic_angle())*u.rad).to(u.deg)
    print(pyephem_q1, pyephem_q2)

    assert (obs.astropy_to_local_time(obs.local_to_astropy_time(dt)).replace(
            tzinfo=None) == dt)

def test_sunrise_sunset_equator():
    '''
    Check that time of sunrise/set for an observer on the equator is
    consistent with PyEphem results (for no atmosphere/pressure=0)
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0 * u.bar
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    obs = Observer(location=location, pressure=pressure)
    astroplan_next_sunrise = obs.sunrise(time, which='next').datetime
    astroplan_next_sunset = obs.sunset(time, which='next').datetime

    astroplan_prev_sunrise = obs.sunrise(time, which='previous').datetime
    astroplan_prev_sunset = obs.sunset(time, which='previous').datetime

    # Run print_pyephem_sunrise_sunset() to compute analogous
    # result from PyEphem:
    pyephem_next_sunrise = datetime.datetime(2000, 1, 2, 6, 3, 39, 150790)
    pyephem_next_sunset = datetime.datetime(2000, 1, 1, 18, 3, 23, 676686)
    pyephem_prev_sunrise = datetime.datetime(2000, 1, 1, 6, 3, 10, 720052)
    pyephem_prev_sunset = datetime.datetime(1999, 12, 31, 18, 2, 55, 100786)

    # Typical difference in this example between PyEphem and astroplan
    # with an atmosphere is <2 min
    threshold_minutes = 8
    assert (abs(pyephem_next_sunrise - astroplan_next_sunrise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_next_sunset - astroplan_next_sunset) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_sunrise - astroplan_prev_sunrise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_sunset - astroplan_prev_sunset) <
            datetime.timedelta(minutes=threshold_minutes))

def print_pyephem_sunrise_sunset():
    '''
    To run:

    python -c 'from astroplan.tests.test_core import print_pyephem_sunrise_sunset as f; f()'
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    time = Time('2000-01-01 12:00:00')

    import ephem
    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.elevation = elevation
    obs.date = time.datetime
    obs.pressure = pressure
    next_sunrise = obs.next_rising(ephem.Sun(), use_center=True)
    next_sunset = obs.next_setting(ephem.Sun(), use_center=True)
    prev_sunrise = obs.previous_rising(ephem.Sun(), use_center=True)
    prev_sunset = obs.previous_setting(ephem.Sun(), use_center=True)

    print(map(repr, [next_sunrise.datetime(), next_sunset.datetime(),
                     prev_sunrise.datetime(), prev_sunset.datetime()]))

def test_vega_rise_set_equator():
    '''
    Check that time of rise/set of Vega for an observer on the equator is
    consistent with PyEphem results (for no atmosphere/pressure=0)
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0 * u.bar
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    vega_ra, vega_dec = (279.23473479*u.degree, 38.78368896*u.degree)
    vega = SkyCoord(vega_ra, vega_dec)

    obs = Observer(location=location, pressure=pressure)
    astroplan_next_rise = obs.calc_rise(time, vega, which='next').datetime
    astroplan_next_set = obs.calc_set(time, vega, which='next').datetime

    astroplan_prev_rise = obs.calc_rise(time, vega, which='previous').datetime
    astroplan_prev_set = obs.calc_set(time, vega, which='previous').datetime

    astroplan_nearest_rise = obs.calc_rise(time, vega, which='nearest').datetime
    astroplan_nearest_set = obs.calc_set(time, vega, which='nearest').datetime

    # Run print_pyephem_vega_rise_set() to compute analogous
    # result from PyEphem:
    pyephem_next_rise = datetime.datetime(2000, 1, 2, 5, 52, 8, 257401)
    pyephem_next_set = datetime.datetime(2000, 1, 1, 17, 54, 6, 211705)
    pyephem_prev_rise = datetime.datetime(2000, 1, 1, 5, 56, 4, 165852)
    pyephem_prev_set = datetime.datetime(1999, 12, 31, 17, 58, 2, 120088)

    # Typical difference in this example between PyEphem and astroplan
    # with an atmosphere is <2 min
    threshold_minutes = 8
    assert (abs(pyephem_next_rise - astroplan_next_rise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_next_set - astroplan_next_set) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_rise - astroplan_prev_rise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_set - astroplan_prev_set) <
            datetime.timedelta(minutes=threshold_minutes))

    # Check that the 'nearest' option selects the nearest rise/set
    assert astroplan_nearest_rise == astroplan_prev_rise
    assert astroplan_nearest_set == astroplan_next_set

def print_pyephem_vega_rise_set():
    '''
    To run:

    python -c 'from astroplan.tests.test_core import print_pyephem_vega_rise_set as f; f()'
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    time = Time('2000-01-01 12:00:00')
    vega_ra, vega_dec = (279.23473479*u.degree, 38.78368896*u.degree)
    vega = SkyCoord(vega_ra, vega_dec)

    import ephem
    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.elevation = elevation
    obs.date = time.datetime
    obs.pressure = pressure
    target = ephem.FixedBody()
    target._ra = ephem.degrees(np.radians(vega.ra.value))
    target._dec = ephem.degrees(np.radians(vega.dec.value))
    target.compute(obs)
    next_rising = obs.next_rising(target).datetime()
    next_setting = obs.next_setting(target).datetime()
    prev_rising = obs.previous_rising(target).datetime()
    prev_setting = obs.previous_setting(target).datetime()

    print(map(repr, [next_rising, next_setting, prev_rising, prev_setting]))

def test_sunrise_sunset_equator_civil_twilight():
    '''
    Check that time of sunrise/set for an observer on the equator is
    consistent with PyEphem results (for no atmosphere/pressure=0)
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0 * u.bar
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    obs = Observer(location=location, pressure=pressure)
    # Manually impose horizon equivalent to civil twilight
    horizon = -6*u.degree
    astroplan_next_sunrise = obs.sunrise(time, which='next',
                                         horizon=horizon).datetime
    astroplan_next_sunset = obs.sunset(time, which='next',
                                       horizon=horizon).datetime

    astroplan_prev_sunrise = obs.sunrise(time, which='previous',
                                         horizon=horizon).datetime
    astroplan_prev_sunset = obs.sunset(time, which='previous',
                                       horizon=horizon).datetime

    # Run print_pyephem_sunrise_sunset_equator_civil_twilight() to compute
    # analogous result from PyEphem:
    pyephem_next_rise = datetime.datetime(2000, 1, 2, 5, 37, 34, 83328)
    pyephem_next_set = datetime.datetime(2000, 1, 1, 18, 29, 29, 195908)
    pyephem_prev_rise = datetime.datetime(2000, 1, 1, 5, 37, 4, 701708)
    pyephem_prev_set = datetime.datetime(1999, 12, 31, 18, 29, 1, 530987)


    threshold_minutes = 8
    assert (abs(pyephem_next_rise - astroplan_next_sunrise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_next_set - astroplan_next_sunset) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_rise - astroplan_prev_sunrise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_set - astroplan_prev_sunset) <
            datetime.timedelta(minutes=threshold_minutes))

def print_pyephem_sunrise_sunset_equator_civil_twilight():
    '''
    Calculate next sunrise and sunset with PyEphem for an observer
    on the equator.

    To run:
    python -c 'from astroplan.tests.test_core import print_pyephem_sunrise_sunset_equator_civil_twilight as f; f()'
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    time = Time('2000-01-01 12:00:00')

    import ephem
    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.elevation = elevation
    obs.date = time.datetime
    obs.pressure = pressure
    obs.horizon = '-06:00:00'
    next_sunrise = obs.next_rising(ephem.Sun(), use_center=True)
    next_sunset = obs.next_setting(ephem.Sun(), use_center=True)
    prev_sunrise = obs.previous_rising(ephem.Sun(), use_center=True)
    prev_sunset = obs.previous_setting(ephem.Sun(), use_center=True)

    pyephem_time_to_datetime_str = lambda t: repr(t.datetime())
    print(map(pyephem_time_to_datetime_str, [next_sunrise, next_sunset,
                                             prev_sunrise, prev_sunset]))

def test_twilight_convenience_funcs():
    '''
    Check that the convenience functions for evening
    astronomical/nautical/civil twilight correspond to their
    PyEphem equivalents
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0 * u.bar
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    obs = Observer(location=location, pressure=pressure)
    # Compute morning twilights with astroplan
    astroplan_morning_civil = obs.morning_civil(time, which='previous').datetime
    astroplan_morning_nautical = obs.morning_nautical(time,
                                                      which='previous').datetime
    astroplan_morning_astro = obs.morning_astronomical(time,
                                                       which='previous').datetime
    # Compute evening twilights with astroplan
    astroplan_evening_civil = obs.evening_civil(time, which='next').datetime
    astroplan_evening_nautical = obs.evening_nautical(time,
                                                      which='next').datetime
    astroplan_evening_astro = obs.evening_astronomical(time,
                                                       which='next').datetime

    # Compute morning and evening twilights with PyEphem from
    # the function print_pyephem_twilight_convenience_funcs()
    pyephem_morning_civil, pyephem_morning_nautical, pyephem_morning_astronomical, = (
        datetime.datetime(2000, 1, 1, 5, 37, 4, 701708),
        datetime.datetime(2000, 1, 1, 5, 10, 55, 450939),
        datetime.datetime(2000, 1, 1, 4, 44, 39, 415865))

    pyephem_evening_civil, pyephem_evening_nautical, pyephem_evening_astronomical = (
        datetime.datetime(2000, 1, 1, 18, 29, 29, 195908),
        datetime.datetime(2000, 1, 1, 18, 55, 37, 864882),
        datetime.datetime(2000, 1, 1, 19, 21, 53, 213768))

    threshold_minutes = 8
    # Compare morning twilights
    assert (abs(astroplan_morning_civil - pyephem_morning_civil) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_morning_nautical - pyephem_morning_nautical) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_morning_astro - pyephem_morning_astronomical) <
            datetime.timedelta(minutes=threshold_minutes))

    # Compare evening twilights
    assert (abs(astroplan_evening_civil - pyephem_evening_civil) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_evening_nautical - pyephem_evening_nautical) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_evening_astro - pyephem_evening_astronomical) <
            datetime.timedelta(minutes=threshold_minutes))

def print_pyephem_twilight_convenience_funcs():
    '''
    To run:
    python -c 'from astroplan.tests.test_core import print_pyephem_twilight_convenience_funcs as f; f()'
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    time = Time('2000-01-01 12:00:00')

    import ephem
    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.elevation = elevation
    obs.date = time.datetime
    obs.pressure = pressure

    # Morning twilights
    obs.horizon = '-06:00:00'
    morning_civil = obs.previous_rising(ephem.Sun(), use_center=True)
    obs.horizon = '-12:00:00'
    morning_nautical = obs.previous_rising(ephem.Sun(), use_center=True)
    obs.horizon = '-18:00:00'
    morning_astronomical = obs.previous_rising(ephem.Sun(), use_center=True)

    # Evening twilights
    obs.horizon = '-06:00:00'
    evening_civil = obs.next_setting(ephem.Sun(), use_center=True)
    obs.horizon = '-12:00:00'
    evening_nautical = obs.next_setting(ephem.Sun(), use_center=True)
    obs.horizon = '-18:00:00'
    evening_astronomical = obs.next_setting(ephem.Sun(), use_center=True)

    pyephem_time_to_datetime_str = lambda t: repr(t.datetime())
    print(map(pyephem_time_to_datetime_str, [morning_civil, morning_nautical,
                                             morning_astronomical,
                                             evening_civil, evening_nautical,
                                             evening_astronomical]))

def test_solar_transit():
    '''
    Test that astroplan's solar transit/antitransit (which are noon and
    midnight) agree with PyEphem's
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0 * u.bar
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    from astropy.coordinates import get_sun
    obs = Observer(location=location, pressure=pressure)

    # Compute next/previous noon/midnight using generic calc_transit methods
    astroplan_next_transit = obs.calc_meridian_transit(time, get_sun(time),
                                                       which='next').datetime
    astroplan_next_antitransit = obs.calc_meridian_antitransit(time,
                                                               get_sun(time),
                                                               which='next').datetime
    astroplan_prev_transit = obs.calc_meridian_transit(time, get_sun(time),
                                                       which='previous').datetime
    astroplan_prev_antitransit = obs.calc_meridian_antitransit(time,
                                                               get_sun(time),
                                                               which='previous').datetime

    astroplan_nearest_transit = obs.calc_meridian_transit(time, get_sun(time),
                                                          which='nearest').datetime
    astroplan_nearest_antitransit = obs.calc_meridian_antitransit(time, get_sun(time),
                                                                  which='nearest').datetime

    # Computed in print_pyephem_solar_transit_noon()
    pyephem_next_transit = datetime.datetime(2000, 1, 1, 12, 3, 17, 207300)
    pyephem_next_antitransit = datetime.datetime(2000, 1, 2, 0, 3, 31, 423333)
    pyephem_prev_transit = datetime.datetime(1999, 12, 31, 12, 2, 48, 562755)
    pyephem_prev_antitransit = datetime.datetime(2000, 1, 1, 0, 3, 2, 918943)

    threshold_minutes = 8
    assert (abs(astroplan_next_transit - pyephem_next_transit) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_next_antitransit - pyephem_next_antitransit) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_prev_transit - pyephem_prev_transit) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_prev_antitransit - pyephem_prev_antitransit) <
            datetime.timedelta(minutes=threshold_minutes))

    # Check nearest
    assert astroplan_next_transit == astroplan_nearest_transit
    assert astroplan_nearest_antitransit == astroplan_prev_antitransit

def test_solar_transit_convenience_methods():
    '''
    Test that astroplan's noon and midnight convenience methods agree with
    PyEphem's solar transit/antitransit time.
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0 * u.bar
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    from astropy.coordinates import get_sun
    obs = Observer(location=location, pressure=pressure)

    # Compute next/previous noon/midnight using generic calc_transit methods
    astroplan_next_noon = obs.noon(time, which='next').datetime
    astroplan_next_midnight = obs.midnight(time, which='next').datetime
    astroplan_prev_noon = obs.noon(time, which='previous').datetime
    astroplan_prev_midnight = obs.midnight(time, which='previous').datetime

    # Computed in print_pyephem_solar_transit_noon()
    pyephem_next_transit = datetime.datetime(2000, 1, 1, 12, 3, 17, 207300)
    pyephem_next_antitransit = datetime.datetime(2000, 1, 2, 0, 3, 31, 423333)
    pyephem_prev_transit = datetime.datetime(1999, 12, 31, 12, 2, 48, 562755)
    pyephem_prev_antitransit = datetime.datetime(2000, 1, 1, 0, 3, 2, 918943)

    threshold_minutes = 8
    assert (abs(astroplan_next_noon - pyephem_next_transit) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_next_midnight - pyephem_next_antitransit) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_prev_noon - pyephem_prev_transit) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(astroplan_prev_midnight - pyephem_prev_antitransit) <
            datetime.timedelta(minutes=threshold_minutes))

def print_pyephem_solar_transit_noon():
    '''
    Calculate next sunrise and sunset with PyEphem for an observer
    on the equator.

    To run:
    python -c 'from astroplan.tests.test_core import print_pyephem_transit_noon as f; f()'
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    time = Time('2000-01-01 12:00:00')

    import ephem
    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.elevation = elevation
    obs.date = time.datetime
    obs.pressure = pressure
    next_transit = obs.next_transit(ephem.Sun())
    next_antitransit = obs.next_antitransit(ephem.Sun())
    prev_transit = obs.previous_transit(ephem.Sun())
    prev_antitransit = obs.previous_antitransit(ephem.Sun())

    pyephem_time_to_datetime_str = lambda t: repr(t.datetime())
    print(map(pyephem_time_to_datetime_str, [next_transit, next_antitransit,
                                             prev_transit, prev_antitransit]))

def test_can_see():
    '''
    Test that Polaris is/isn't observable from north/south pole
    '''
    elevation = 0.0 * u.m
    pressure = 0 * u.bar
    north = EarthLocation.from_geodetic('00:00:00',
                                        '90:00:00', elevation)
    south = EarthLocation.from_geodetic('00:00:00',
                                        '-90:00:00', elevation)
    time = Time('2000-01-01 12:00:00')
    polaris = SkyCoord(37.95456067*u.degree, 89.26410897*u.degree)
    north_pole = Observer(location=north, pressure=pressure)
    south_pole = Observer(location=south, pressure=pressure)
    assert north_pole.can_see(time, polaris)
    assert not south_pole.can_see(time, polaris)

def test_string_times():
    '''
    Test that strings passed to time argument get successfully
    passed to Time constructor. Analogous test to test_vega_rise_set_equator(),
    just with a string for a time.
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0 * u.bar
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = '2000-01-01 12:00:00'
    vega_ra, vega_dec = (279.23473479*u.degree, 38.78368896*u.degree)
    vega = SkyCoord(vega_ra, vega_dec)

    obs = Observer(location=location, pressure=pressure)
    astroplan_next_rise = obs.calc_rise(time, vega, which='next').datetime
    astroplan_next_set = obs.calc_set(time, vega, which='next').datetime

    astroplan_prev_rise = obs.calc_rise(time, vega, which='previous').datetime
    astroplan_prev_set = obs.calc_set(time, vega, which='previous').datetime

    astroplan_nearest_rise = obs.calc_rise(time, vega, which='nearest').datetime
    astroplan_nearest_set = obs.calc_set(time, vega, which='nearest').datetime

    # Run print_pyephem_vega_rise_set() to compute analogous
    # result from PyEphem:
    pyephem_next_rise = datetime.datetime(2000, 1, 2, 5, 52, 8, 257401)
    pyephem_next_set = datetime.datetime(2000, 1, 1, 17, 54, 6, 211705)
    pyephem_prev_rise = datetime.datetime(2000, 1, 1, 5, 56, 4, 165852)
    pyephem_prev_set = datetime.datetime(1999, 12, 31, 17, 58, 2, 120088)

    # Typical difference in this example between PyEphem and astroplan
    # with an atmosphere is <2 min
    threshold_minutes = 8
    assert (abs(pyephem_next_rise - astroplan_next_rise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_next_set - astroplan_next_set) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_rise - astroplan_prev_rise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_set - astroplan_prev_set) <
            datetime.timedelta(minutes=threshold_minutes))

    # Check that the 'nearest' option selects the nearest rise/set
    assert astroplan_nearest_rise == astroplan_prev_rise
    assert astroplan_nearest_set == astroplan_next_set

def test_TargetAlwaysUpWarning(recwarn):
    lat = '90:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    polaris = SkyCoord(37.95456067*u.degree, 89.26410897*u.degree)

    obs = Observer(location=location)
    no_time = obs.calc_rise(time, polaris, which='next')

    w = recwarn.pop(TargetAlwaysUpWarning)
    assert issubclass(w.category, TargetAlwaysUpWarning)
    assert np.isnan(no_time)

def test_TargetNeverUpWarning(recwarn):
    lat = '-90:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    polaris = SkyCoord(37.95456067*u.degree, 89.26410897*u.degree)

    obs = Observer(location=location)
    no_time = obs.calc_rise(time, polaris, which='next')

    w = recwarn.pop(TargetNeverUpWarning)
    assert issubclass(w.category, TargetNeverUpWarning)
    assert np.isnan(no_time)

def test_timezone_convenience_methods():
    location = EarthLocation(-74.0*u.deg, 40.7*u.deg, 0*u.m)
    obs = Observer(location=location,timezone=pytz.timezone('US/Eastern'))
    t = Time(57100.3, format='mjd')
    assert (obs.astropy_time_to_datetime(t).hour == 3)

    dt = datetime.datetime(2015, 3, 19, 3, 12)
    assert (obs.datetime_to_astropy_time(dt).datetime ==
            datetime.datetime(2015, 3, 19, 7, 12))

    assert (obs.astropy_time_to_datetime(obs.datetime_to_astropy_time(dt)).replace(
            tzinfo=None) == dt)

    # Test ndarray of times:
    times = t + np.linspace(0, 24, 10)*u.hour
    times_dt_ndarray = times.datetime
    assert all((obs.datetime_to_astropy_time(times_dt_ndarray)).jd ==
               (times + 4*u.hour).jd)

    # Test list of times:
    times_dt_list = list(times.datetime)
    assert all((obs.datetime_to_astropy_time(times_dt_list)).jd ==
               (times + 4*u.hour).jd)

    dts = obs.astropy_time_to_datetime(times)
    naive_dts = list(map(lambda t: t.replace(tzinfo=None), dts))
    assert all(naive_dts == times_dt_ndarray - datetime.timedelta(hours=4))

class TestExceptions(unittest.TestCase):
    def test_rise_set_transit_which(self):
        lat = '00:00:00'
        lon = '00:00:00'
        elevation = 0.0 * u.m
        location = EarthLocation.from_geodetic(lon, lat, elevation)
        time = Time('2000-01-01 12:00:00')
        vega_coords = SkyCoord('18h36m56.33635s', '+38d47m01.2802s')

        obs = Observer(location=location)

        with self.assertRaises(ValueError):
            obs.calc_rise(time, vega_coords, which='oops').datetime

        with self.assertRaises(ValueError):
            obs.calc_set(time, vega_coords, which='oops').datetime

        with self.assertRaises(ValueError):
            obs.calc_meridian_transit(time, vega_coords, which='oops').datetime

        with self.assertRaises(ValueError):
            obs.calc_meridian_antitransit(time, vega_coords, which='oops').datetime

    def test_FixedTarget_duck_typing(self):
        with self.assertRaises(TypeError):
            FixedTarget(['00:00:00', '00:00:00'], name='VE')

    def test_Observer_init(self):
        with self.assertRaises(TypeError):
            Observer(location='Greenwich')

        with self.assertRaises(TypeError):
            Observer(location=EarthLocation(0, 0, 0), timezone=-6)

    def test_Observer_altaz(self):
        with self.assertRaises(TypeError):
            obs = Observer(location=EarthLocation(0, 0, 0))
            obs.altaz(Time('2000-01-01 00:00:00'), ['00:00:00','00:00:00'])

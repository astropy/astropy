# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the direct ITRS to observed transformations.

These tests verify the new ITRS<->AltAz and ITRS<->HADec transformations that
stay within the ITRS coordinate system, which is appropriate for nearby objects
like satellites, aircraft, and ground features.
"""

import numpy as np
import pytest

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation, ITRS, AltAz, HADec, CartesianRepresentation
)


def test_itrs_to_altaz_straight_up():
    """
    Test ITRS to AltAz for an object directly overhead.

    This is a basic sanity check - an object at a position directly above
    the observer in ITRS should appear at altitude 90 degrees.
    """
    # Define observer location
    lon = 0 * u.deg
    lat = 45 * u.deg
    height = 0 * u.m
    location = EarthLocation(lon=lon, lat=lat, height=height)

    # Create an ITRS position 100 km directly above the observer
    # First get the observer's ITRS position
    obstime = Time('2020-01-01T00:00:00')
    observer_itrs = location.get_itrs(obstime)

    # Calculate the local vertical direction (pointing up)
    # In ITRS, this is along the radial direction from Earth center
    vertical_direction = observer_itrs.cartesian.xyz / np.linalg.norm(observer_itrs.cartesian.xyz.value)

    # Place object 100 km above observer
    object_itrs = ITRS(
        CartesianRepresentation(
            observer_itrs.cartesian.xyz + vertical_direction * 100 * u.km
        ),
        obstime=obstime
    )

    # Transform to AltAz
    altaz_frame = AltAz(obstime=obstime, location=location)
    object_altaz = object_itrs.transform_to(altaz_frame)

    # Should be at altitude 90 degrees (straight up)
    assert_allclose(object_altaz.alt, 90 * u.deg, atol=1e-10 * u.deg)
    # Distance should be 100 km
    assert_allclose(object_altaz.distance, 100 * u.km, rtol=1e-10)


def test_itrs_to_hadec_straight_up():
    """
    Test ITRS to HADec for an object directly overhead.

    An object directly overhead should have hour angle 0 and declination
    equal to the observer's latitude.
    """
    # Define observer location
    lon = 0 * u.deg
    lat = 45 * u.deg
    height = 0 * u.m
    location = EarthLocation(lon=lon, lat=lat, height=height)

    # Create an ITRS position 100 km directly above the observer
    obstime = Time('2020-01-01T00:00:00')
    observer_itrs = location.get_itrs(obstime)

    # Calculate the local vertical direction
    vertical_direction = observer_itrs.cartesian.xyz / np.linalg.norm(observer_itrs.cartesian.xyz.value)

    # Place object 100 km above observer
    object_itrs = ITRS(
        CartesianRepresentation(
            observer_itrs.cartesian.xyz + vertical_direction * 100 * u.km
        ),
        obstime=obstime
    )

    # Transform to HADec
    hadec_frame = HADec(obstime=obstime, location=location)
    object_hadec = object_itrs.transform_to(hadec_frame)

    # Should have hour angle 0
    assert_allclose(object_hadec.ha, 0 * u.hourangle, atol=1e-10 * u.hourangle)
    # Declination should equal observer's latitude
    assert_allclose(object_hadec.dec, lat, atol=1e-10 * u.deg)
    # Distance should be 100 km
    assert_allclose(object_hadec.distance, 100 * u.km, rtol=1e-10)


def test_itrs_to_altaz_on_horizon():
    """
    Test ITRS to AltAz for an object on the horizon.

    An object at the same height as the observer but displaced horizontally
    should appear on or near the horizon.
    """
    # Define observer location
    lon = 0 * u.deg
    lat = 45 * u.deg
    height = 0 * u.m
    location = EarthLocation(lon=lon, lat=lat, height=height)

    obstime = Time('2020-01-01T00:00:00')
    observer_itrs = location.get_itrs(obstime)

    # Create a position at the same radius from Earth center but displaced
    # in the direction of north (along the +z direction in ITRS, then rotated)
    # We'll place it 100 km north along the surface
    # This is approximate - for testing purposes

    # For simplicity, use a position directly north at horizon level
    # Get a point on the horizon by going perpendicular to the vertical
    obs_x, obs_y, obs_z = observer_itrs.cartesian.xyz.to(u.m).value

    # Calculate tangent direction (north)
    # In ITRS, the north direction at a point is approximately along increasing z
    # but we need to be tangent to the sphere
    # Use a simpler approach: offset in x direction (east) for testing

    # Place object 100 km to the east at same distance from Earth center
    # This should be near the horizon
    from astropy.coordinates.matrix_utilities import rotation_matrix

    # Rotate the observer position slightly around z-axis (longitude)
    # Small angle to stay near horizon
    angle = 1 * u.deg
    rot_mat = rotation_matrix(angle, 'z')
    new_pos = rot_mat @ observer_itrs.cartesian.xyz

    object_itrs = ITRS(CartesianRepresentation(new_pos), obstime=obstime)

    # Transform to AltAz
    altaz_frame = AltAz(obstime=obstime, location=location)
    object_altaz = object_itrs.transform_to(altaz_frame)

    # Should be near the horizon (altitude close to 0, but not exactly due to Earth curvature)
    # The exact value depends on the geometry, but should be relatively small
    assert object_altaz.alt < 10 * u.deg
    assert object_altaz.alt > -10 * u.deg


def test_itrs_altaz_roundtrip():
    """
    Test round-trip transformation ITRS -> AltAz -> ITRS.

    The coordinates should remain unchanged after a round trip.
    """
    # Define observer location
    location = EarthLocation(lon=10 * u.deg, lat=45 * u.deg, height=0 * u.m)
    obstime = Time('2020-01-01T00:00:00')

    # Create an ITRS position (a satellite 400 km above a point)
    observer_itrs = location.get_itrs(obstime)
    vertical_direction = observer_itrs.cartesian.xyz / np.linalg.norm(observer_itrs.cartesian.xyz.value)
    object_itrs = ITRS(
        CartesianRepresentation(
            observer_itrs.cartesian.xyz + vertical_direction * 400 * u.km
        ),
        obstime=obstime
    )

    # Transform to AltAz and back
    altaz_frame = AltAz(obstime=obstime, location=location)
    object_altaz = object_itrs.transform_to(altaz_frame)
    object_itrs_2 = object_altaz.transform_to(ITRS(obstime=obstime))

    # Check round-trip accuracy
    assert_allclose(object_itrs.cartesian.xyz, object_itrs_2.cartesian.xyz, rtol=1e-10)


def test_itrs_hadec_roundtrip():
    """
    Test round-trip transformation ITRS -> HADec -> ITRS.

    The coordinates should remain unchanged after a round trip.
    """
    # Define observer location
    location = EarthLocation(lon=10 * u.deg, lat=45 * u.deg, height=0 * u.m)
    obstime = Time('2020-01-01T00:00:00')

    # Create an ITRS position
    observer_itrs = location.get_itrs(obstime)
    vertical_direction = observer_itrs.cartesian.xyz / np.linalg.norm(observer_itrs.cartesian.xyz.value)
    object_itrs = ITRS(
        CartesianRepresentation(
            observer_itrs.cartesian.xyz + vertical_direction * 400 * u.km
        ),
        obstime=obstime
    )

    # Transform to HADec and back
    hadec_frame = HADec(obstime=obstime, location=location)
    object_hadec = object_itrs.transform_to(hadec_frame)
    object_itrs_2 = object_hadec.transform_to(ITRS(obstime=obstime))

    # Check round-trip accuracy
    assert_allclose(object_itrs.cartesian.xyz, object_itrs_2.cartesian.xyz, rtol=1e-10)


def test_itrs_to_altaz_different_locations():
    """
    Test that the same ITRS position gives different AltAz coordinates
    for different observer locations.
    """
    # Fixed ITRS position (e.g., a geostationary satellite)
    obstime = Time('2020-01-01T00:00:00')
    # Geostationary altitude is ~35,786 km above equator
    sat_itrs = ITRS(
        CartesianRepresentation(x=42164 * u.km, y=0 * u.km, z=0 * u.km),
        obstime=obstime
    )

    # Two different observer locations
    location1 = EarthLocation(lon=0 * u.deg, lat=0 * u.deg, height=0 * u.m)
    location2 = EarthLocation(lon=90 * u.deg, lat=0 * u.deg, height=0 * u.m)

    # Transform to AltAz for both locations
    altaz1 = sat_itrs.transform_to(AltAz(obstime=obstime, location=location1))
    altaz2 = sat_itrs.transform_to(AltAz(obstime=obstime, location=location2))

    # The coordinates should be different
    assert not np.allclose(altaz1.alt.value, altaz2.alt.value, rtol=0.01)
    assert not np.allclose(altaz1.az.value, altaz2.az.value, rtol=0.01)


def test_itrs_obstime_handling():
    """
    Test that ITRS coordinates are treated as time-invariant.

    When transforming from ITRS to AltAz with different obstimes, the ITRS
    position should be treated as fixed to the Earth, not moving with the
    solar system barycenter.
    """
    # Define a fixed location on Earth
    location = EarthLocation(lon=0 * u.deg, lat=45 * u.deg, height=0 * u.m)

    # Create an ITRS position at one time
    obstime1 = Time('2020-01-01T00:00:00')
    observer_itrs = location.get_itrs(obstime1)
    vertical_direction = observer_itrs.cartesian.xyz / np.linalg.norm(observer_itrs.cartesian.xyz.value)
    object_itrs = ITRS(
        CartesianRepresentation(
            observer_itrs.cartesian.xyz + vertical_direction * 100 * u.km
        ),
        obstime=obstime1
    )

    # Transform to AltAz at a much later time
    obstime2 = Time('2020-01-02T00:00:00')
    altaz_frame = AltAz(obstime=obstime2, location=location)
    object_altaz = object_itrs.transform_to(altaz_frame)

    # The object should still appear near zenith (within reason, accounting for
    # Earth's rotation). If ITRS were incorrectly transformed through the SSB,
    # the position would be wildly off.
    # Note: Due to Earth's rotation, the altitude won't be exactly 90 degrees,
    # but it should still be high in the sky.
    # This is more of a sanity check that we haven't moved millions of km away.
    assert object_altaz.alt > 0 * u.deg


def test_altaz_hadec_consistency():
    """
    Test that ITRS->AltAz and ITRS->HADec give consistent results.

    For an object straight overhead, both should agree that it's at the zenith.
    """
    # Define observer location
    location = EarthLocation(lon=0 * u.deg, lat=45 * u.deg, height=0 * u.m)
    obstime = Time('2020-01-01T00:00:00')

    # Create an ITRS position directly overhead
    observer_itrs = location.get_itrs(obstime)
    vertical_direction = observer_itrs.cartesian.xyz / np.linalg.norm(observer_itrs.cartesian.xyz.value)
    object_itrs = ITRS(
        CartesianRepresentation(
            observer_itrs.cartesian.xyz + vertical_direction * 100 * u.km
        ),
        obstime=obstime
    )

    # Transform to both AltAz and HADec
    altaz = object_itrs.transform_to(AltAz(obstime=obstime, location=location))
    hadec = object_itrs.transform_to(HADec(obstime=obstime, location=location))

    # AltAz should show altitude 90 degrees
    assert_allclose(altaz.alt, 90 * u.deg, atol=1e-10 * u.deg)

    # HADec should show HA=0 and dec=latitude
    assert_allclose(hadec.ha, 0 * u.hourangle, atol=1e-10 * u.hourangle)
    assert_allclose(hadec.dec, 45 * u.deg, atol=1e-10 * u.deg)


def test_satellite_tracking_scenario():
    """
    Test a realistic satellite tracking scenario.

    This simulates tracking a satellite in low Earth orbit from a ground station.
    """
    # Ground station location
    location = EarthLocation(lon=-104.0 * u.deg, lat=40.0 * u.deg, height=1650 * u.m)
    obstime = Time('2020-06-01T12:00:00')

    # Simulate a satellite in LEO at ~500 km altitude
    # Place it at a specific ITRS position
    # For this test, we'll use an arbitrary but reasonable position
    sat_itrs = ITRS(
        CartesianRepresentation(x=5000 * u.km, y=3000 * u.km, z=3000 * u.km),
        obstime=obstime
    )

    # Transform to AltAz to get pointing information
    altaz = sat_itrs.transform_to(AltAz(obstime=obstime, location=location))

    # Basic sanity checks
    assert altaz.alt > -90 * u.deg
    assert altaz.alt < 90 * u.deg
    assert altaz.az >= 0 * u.deg
    assert altaz.az < 360 * u.deg
    assert altaz.distance > 0 * u.km

    # Check that round-trip works
    sat_itrs_2 = altaz.transform_to(ITRS(obstime=obstime))
    assert_allclose(sat_itrs.cartesian.xyz, sat_itrs_2.cartesian.xyz, rtol=1e-10)


def test_multiple_objects():
    """
    Test transforming multiple objects at once using arrays.
    """
    # Define observer location
    location = EarthLocation(lon=0 * u.deg, lat=45 * u.deg, height=0 * u.m)
    obstime = Time('2020-01-01T00:00:00')

    # Create multiple ITRS positions (e.g., a constellation of satellites)
    observer_itrs = location.get_itrs(obstime)
    vertical_direction = observer_itrs.cartesian.xyz / np.linalg.norm(observer_itrs.cartesian.xyz.value)

    # Create 5 objects at different altitudes above the observer
    altitudes = np.array([100, 200, 300, 400, 500]) * u.km
    positions = observer_itrs.cartesian.xyz[:, np.newaxis] + vertical_direction[:, np.newaxis] * altitudes

    objects_itrs = ITRS(CartesianRepresentation(positions), obstime=obstime)

    # Transform to AltAz
    altaz_frame = AltAz(obstime=obstime, location=location)
    objects_altaz = objects_itrs.transform_to(altaz_frame)

    # All should be at altitude 90 degrees
    assert_allclose(objects_altaz.alt, 90 * u.deg, atol=1e-10 * u.deg)

    # Distances should match the altitudes
    assert_allclose(objects_altaz.distance, altitudes, rtol=1e-10)

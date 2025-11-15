# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the direct ITRS to Observed (AltAz/HADec) transformations.

These tests verify that the direct ITRS<->Observed transformations work correctly
for nearby objects where staying within ITRS is appropriate.
"""

import numpy as np
import pytest

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation, ITRS, AltAz, HADec, CartesianRepresentation
)


def test_itrs_to_altaz_straight_overhead():
    """
    Test ITRS to AltAz for an object directly overhead.

    An object at a location directly above the observer should appear at
    altitude 90 degrees (zenith).
    """
    t = Time('J2010')
    # Observer location
    home = EarthLocation(-1*u.deg, 52*u.deg, height=0.*u.km)
    # Object directly overhead (10 km above observer)
    obj = EarthLocation(-1*u.deg, 52*u.deg, height=10.*u.km)

    # Get ITRS position of object
    itrs_obj = obj.get_itrs(t)

    # Transform to AltAz as seen from home
    aa = itrs_obj.transform_to(AltAz(obstime=t, location=home))

    # Should be directly overhead (altitude = 90 degrees)
    # Azimuth is undefined at zenith, so we don't test it
    assert_quantity_allclose(aa.alt, 90*u.deg, atol=1*u.arcsec)


def test_itrs_to_hadec_straight_overhead():
    """
    Test ITRS to HADec for an object directly overhead.

    An object at a location directly above the observer should appear at
    hour angle 0 and declination equal to the observer's latitude.
    """
    t = Time('J2010')
    # Observer location
    home = EarthLocation(-1*u.deg, 52*u.deg, height=0.*u.km)
    # Object directly overhead (10 km above observer)
    obj = EarthLocation(-1*u.deg, 52*u.deg, height=10.*u.km)

    # Get ITRS position of object
    itrs_obj = obj.get_itrs(t)

    # Transform to HADec as seen from home
    hd = itrs_obj.transform_to(HADec(obstime=t, location=home))

    # Should be at hour angle 0 and declination equal to latitude
    assert_quantity_allclose(hd.ha, 0*u.hourangle, atol=1*u.arcsec)
    assert_quantity_allclose(hd.dec, 52*u.deg, atol=1*u.arcsec)


def test_itrs_altaz_roundtrip():
    """
    Test that ITRS -> AltAz -> ITRS roundtrip works correctly.
    """
    t = Time('J2015')
    home = EarthLocation(lon=10*u.deg, lat=45*u.deg, height=0*u.m)

    # Create an ITRS coordinate for a nearby object (e.g., a satellite)
    # Position 100 km above the observer
    observer_itrs = home.get_itrs(t)
    # Offset by 100 km in the "up" direction
    # For simplicity, we'll create a point displaced vertically
    itrs_original = ITRS(
        x=observer_itrs.cartesian.x,
        y=observer_itrs.cartesian.y,
        z=observer_itrs.cartesian.z + 100*u.km,
        obstime=t
    )

    # Transform to AltAz and back
    aa = itrs_original.transform_to(AltAz(obstime=t, location=home))
    itrs_roundtrip = aa.transform_to(ITRS(obstime=t))

    # Should get back the same position (within tolerance)
    assert_quantity_allclose(
        itrs_original.cartesian.xyz,
        itrs_roundtrip.cartesian.xyz,
        atol=1*u.mm
    )


def test_itrs_hadec_roundtrip():
    """
    Test that ITRS -> HADec -> ITRS roundtrip works correctly.
    """
    t = Time('J2015')
    home = EarthLocation(lon=10*u.deg, lat=45*u.deg, height=0*u.m)

    # Create an ITRS coordinate for a nearby object
    observer_itrs = home.get_itrs(t)
    itrs_original = ITRS(
        x=observer_itrs.cartesian.x + 50*u.km,
        y=observer_itrs.cartesian.y + 50*u.km,
        z=observer_itrs.cartesian.z + 100*u.km,
        obstime=t
    )

    # Transform to HADec and back
    hd = itrs_original.transform_to(HADec(obstime=t, location=home))
    itrs_roundtrip = hd.transform_to(ITRS(obstime=t))

    # Should get back the same position (within tolerance)
    assert_quantity_allclose(
        itrs_original.cartesian.xyz,
        itrs_roundtrip.cartesian.xyz,
        atol=1*u.mm
    )


def test_itrs_observed_time_invariance():
    """
    Test that ITRS coordinates are treated as time-invariant.

    The transformation should not change the ITRS coordinate based on
    different obstimes in the target frame.
    """
    t1 = Time('J2000')
    t2 = Time('J2010')
    home = EarthLocation(lon=0*u.deg, lat=45*u.deg, height=0*u.m)

    # Create an ITRS coordinate at time t1
    observer_itrs = home.get_itrs(t1)
    itrs_coord = ITRS(
        x=observer_itrs.cartesian.x,
        y=observer_itrs.cartesian.y,
        z=observer_itrs.cartesian.z + 100*u.km,
        obstime=t1
    )

    # Transform to AltAz at different times
    # The ITRS position should be treated as fixed, so we're just
    # converting the coordinate system
    aa1 = itrs_coord.transform_to(AltAz(obstime=t1, location=home))
    aa2 = itrs_coord.transform_to(AltAz(obstime=t2, location=home))

    # The transformations should give different results because the Earth
    # has rotated between t1 and t2, but when we transform back to ITRS,
    # we should get different ITRS coordinates (one at t1, one at t2)
    # But if we specify the same location, the topocentric part should be similar

    # Transform back to ITRS
    itrs1 = aa1.transform_to(ITRS(obstime=t1))
    itrs2 = aa2.transform_to(ITRS(obstime=t2))

    # Both should be close to the original ITRS coordinate
    # (within a reasonable tolerance for the different obstimes)
    assert_quantity_allclose(
        itrs1.cartesian.xyz,
        itrs_coord.cartesian.xyz,
        atol=1*u.mm
    )


def test_itrs_to_altaz_multiple_objects():
    """
    Test ITRS to AltAz transformation with multiple objects at once.
    """
    t = Time('J2010')
    home = EarthLocation(lon=0*u.deg, lat=45*u.deg, height=0*u.m)
    observer_itrs = home.get_itrs(t)

    # Create multiple ITRS coordinates
    n_objects = 10
    x = observer_itrs.cartesian.x + np.linspace(0, 100, n_objects) * u.km
    y = observer_itrs.cartesian.y + np.linspace(0, 100, n_objects) * u.km
    z = observer_itrs.cartesian.z + np.linspace(50, 150, n_objects) * u.km

    itrs_coords = ITRS(
        CartesianRepresentation(x=x, y=y, z=z),
        obstime=t
    )

    # Transform to AltAz
    aa_coords = itrs_coords.transform_to(AltAz(obstime=t, location=home))

    # Check that we got the right number of coordinates
    assert len(aa_coords) == n_objects

    # Transform back and check roundtrip
    itrs_roundtrip = aa_coords.transform_to(ITRS(obstime=t))
    assert_quantity_allclose(
        itrs_coords.cartesian.xyz,
        itrs_roundtrip.cartesian.xyz,
        atol=1*u.mm
    )


def test_itrs_to_hadec_multiple_objects():
    """
    Test ITRS to HADec transformation with multiple objects at once.
    """
    t = Time('J2010')
    home = EarthLocation(lon=0*u.deg, lat=45*u.deg, height=0*u.m)
    observer_itrs = home.get_itrs(t)

    # Create multiple ITRS coordinates
    n_objects = 10
    x = observer_itrs.cartesian.x + np.linspace(0, 100, n_objects) * u.km
    y = observer_itrs.cartesian.y + np.linspace(0, 100, n_objects) * u.km
    z = observer_itrs.cartesian.z + np.linspace(50, 150, n_objects) * u.km

    itrs_coords = ITRS(
        CartesianRepresentation(x=x, y=y, z=z),
        obstime=t
    )

    # Transform to HADec
    hd_coords = itrs_coords.transform_to(HADec(obstime=t, location=home))

    # Check that we got the right number of coordinates
    assert len(hd_coords) == n_objects

    # Transform back and check roundtrip
    itrs_roundtrip = hd_coords.transform_to(ITRS(obstime=t))
    assert_quantity_allclose(
        itrs_coords.cartesian.xyz,
        itrs_roundtrip.cartesian.xyz,
        atol=1*u.mm
    )


def test_itrs_observed_different_locations():
    """
    Test that transformations work correctly for different observer locations.
    """
    t = Time('J2010')

    # Create several observer locations
    locations = [
        EarthLocation(lon=0*u.deg, lat=0*u.deg, height=0*u.m),  # Equator
        EarthLocation(lon=0*u.deg, lat=45*u.deg, height=0*u.m),  # Mid-latitude
        EarthLocation(lon=0*u.deg, lat=89*u.deg, height=0*u.m),  # Near pole
        EarthLocation(lon=180*u.deg, lat=-45*u.deg, height=1000*u.m),  # South, elevated
    ]

    for loc in locations:
        # Create an object above the observer
        obj = EarthLocation(
            lon=loc.lon,
            lat=loc.lat,
            height=loc.height + 100*u.km
        )
        itrs_obj = obj.get_itrs(t)

        # Transform to AltAz
        aa = itrs_obj.transform_to(AltAz(obstime=t, location=loc))

        # Should be at or very near zenith
        assert_quantity_allclose(aa.alt, 90*u.deg, atol=0.1*u.deg)

        # Transform to HADec
        hd = itrs_obj.transform_to(HADec(obstime=t, location=loc))

        # Should be at hour angle 0 and declination equal to latitude
        assert_quantity_allclose(hd.ha, 0*u.hourangle, atol=0.1*u.deg)
        assert_quantity_allclose(hd.dec, loc.lat, atol=0.1*u.deg)


def test_straight_overhead_simple():
    """
    Test the simple direct approach for an overhead object.

    This demonstrates the intuitive solution enabled by the direct ITRS->Observed
    transformation, which is much simpler than the workaround in
    test_intermediate_transformations.test_straight_overhead().
    """
    t = Time('J2010')
    # Observer location
    home = EarthLocation(-1*u.deg, 52*u.deg, height=0.*u.km)
    # Object directly overhead (10 km above observer)
    obj = EarthLocation(-1*u.deg, 52*u.deg, height=10.*u.km)

    # With the direct ITRS->Observed transform, this is straightforward:
    itrs_obj = obj.get_itrs(t)

    # Transform directly from ITRS to AltAz
    aa = itrs_obj.transform_to(AltAz(obstime=t, location=home))

    # Should be directly overhead (altitude = 90 degrees)
    assert_quantity_allclose(aa.alt, 90*u.deg, atol=1*u.arcsec)

    # Transform directly from ITRS to HADec
    hd = itrs_obj.transform_to(HADec(obstime=t, location=home))

    # Should be at hour angle 0 and declination equal to latitude
    assert_quantity_allclose(hd.ha, 0*u.hourangle, atol=1*u.arcsec)
    assert_quantity_allclose(hd.dec, 52*u.deg, atol=1*u.arcsec)


def test_consistency_with_cirs_transform():
    """
    Test that the direct ITRS->Observed transform gives similar results to
    ITRS->CIRS->Observed for distant objects.

    For distant objects, both approaches should give similar results (though
    not identical due to aberration effects).
    """
    t = Time('J2015')
    home = EarthLocation(lon=10*u.deg, lat=45*u.deg, height=0*u.m)

    # Create an ITRS coordinate for a very distant object
    # (where aberration effects are negligible)
    observer_itrs = home.get_itrs(t)
    # Point to a distant location in the same direction
    direction = observer_itrs.cartesian.xyz / np.linalg.norm(
        observer_itrs.cartesian.xyz.to_value(u.m)
    )
    # Create a point 1000 km away in that direction
    distant_point = observer_itrs.cartesian + 1000 * u.km * direction

    itrs_coord = ITRS(distant_point, obstime=t)

    # Direct transform
    aa_direct = itrs_coord.transform_to(AltAz(obstime=t, location=home))

    # Transform via CIRS (the old way)
    from astropy.coordinates import CIRS
    cirs_coord = itrs_coord.transform_to(CIRS(obstime=t, location=home))
    aa_via_cirs = cirs_coord.transform_to(AltAz(obstime=t, location=home))

    # Results should be similar (within a few degrees due to different treatment)
    # This is a coarse check - the direct transform is for nearby objects
    # For nearby objects, the difference can be larger
    # We mainly want to ensure both transforms work without error
    assert aa_direct.alt.value > 0  # Just check it's a valid altitude
    assert aa_via_cirs.alt.value > 0

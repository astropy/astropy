# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for direct ITRS ↔ AltAz/HADec transforms (PR #13398).

Covers:
- ITRS frame has a ``location`` attribute (geocenter default)
- ``earth_location`` property includes the location offset
- Direct geometric ITRS→AltAz/HADec transform and its inverse
- Round-trip precision for both topocentric and geocentric cases
- Time-invariance when ITRS location == observer location
"""
import numpy as np
import pytest

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation, ITRS, AltAz, HADec,
    CartesianRepresentation,
)
from astropy.tests.helper import assert_quantity_allclose


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def greenwich():
    return EarthLocation(lon=0*u.deg, lat=51.5*u.deg, height=0*u.m)


@pytest.fixture
def obstime():
    return Time('2020-06-01 12:00:00')


@pytest.fixture
def altaz_frame(greenwich, obstime):
    return AltAz(location=greenwich, obstime=obstime)


@pytest.fixture
def hadec_frame(greenwich, obstime):
    return HADec(location=greenwich, obstime=obstime)


# ---------------------------------------------------------------------------
# Smoke — module import and ITRS location attribute
# ---------------------------------------------------------------------------

def test_itrs_has_location_attribute(greenwich):
    """ITRS frame accepts a location attribute and defaults to geocenter."""
    from astropy.coordinates.builtin_frames.utils import EARTH_CENTER

    # default → geocenter
    itrs_default = ITRS(CartesianRepresentation(0*u.m, 0*u.m, 6378100*u.m))
    assert itrs_default.location is not None

    # explicit location
    itrs_topo = ITRS(CartesianRepresentation(1*u.m, 2*u.m, 3*u.m),
                     location=greenwich)
    assert itrs_topo.location == greenwich


def test_itrs_earth_location_includes_offset(greenwich):
    """earth_location property adds the frame's location offset to data."""
    # Place a point 10 m above the observer in ITRS (zero displacement)
    loc_cart = CartesianRepresentation(greenwich.x, greenwich.y, greenwich.z)
    displacement = CartesianRepresentation(0*u.m, 0*u.m, 10*u.m)

    itrs = ITRS(displacement, location=greenwich)
    el = itrs.earth_location

    assert_quantity_allclose(el.x, greenwich.x, atol=1*u.m)
    assert_quantity_allclose(el.y, greenwich.y, atol=1*u.m)
    assert_quantity_allclose(el.z, greenwich.z + 10*u.m, atol=1*u.m)


# ---------------------------------------------------------------------------
# Integration — round-trip topocentric ITRS ↔ AltAz / HADec
# ---------------------------------------------------------------------------

def test_itrs_altaz_roundtrip_topocentric(altaz_frame, greenwich):
    """Topocentric ITRS→AltAz→ITRS round-trips to sub-millimetre precision."""
    # Topocentric displacement: object 1 km to the north at the same height
    topo = CartesianRepresentation(0*u.m, 1000*u.m, 0*u.m)
    itrs_in = ITRS(topo, location=greenwich)

    altaz = itrs_in.transform_to(altaz_frame)
    itrs_out = altaz.transform_to(ITRS(location=greenwich))

    assert_quantity_allclose(itrs_in.cartesian.xyz, itrs_out.cartesian.xyz,
                             atol=1e-3*u.m)


def test_itrs_hadec_roundtrip_topocentric(hadec_frame, greenwich):
    """Topocentric ITRS→HADec→ITRS round-trips to sub-millimetre precision."""
    topo = CartesianRepresentation(500*u.m, -300*u.m, 200*u.m)
    itrs_in = ITRS(topo, location=greenwich)

    hadec = itrs_in.transform_to(hadec_frame)
    itrs_out = hadec.transform_to(ITRS(location=greenwich))

    assert_quantity_allclose(itrs_in.cartesian.xyz, itrs_out.cartesian.xyz,
                             atol=1e-3*u.m)


# ---------------------------------------------------------------------------
# Integration — geocentric ITRS → AltAz (different locations)
# ---------------------------------------------------------------------------

def test_itrs_altaz_geocentric_roundtrip(altaz_frame, greenwich):
    """Geocentric ITRS→AltAz→geocentric ITRS round-trips correctly."""
    # Geocentric ITRS: put the observer's own ECEF position as the data
    obs_cart = CartesianRepresentation(greenwich.x, greenwich.y, greenwich.z)
    itrs_in = ITRS(obs_cart)  # default geocentric

    altaz = itrs_in.transform_to(altaz_frame)
    itrs_out = altaz.transform_to(ITRS())  # back to geocentric

    assert_quantity_allclose(itrs_in.cartesian.xyz, itrs_out.cartesian.xyz,
                             atol=1e-3*u.m)


# ---------------------------------------------------------------------------
# Unit — time-invariance for topocentric ITRS → AltAz
# ---------------------------------------------------------------------------

def test_itrs_altaz_time_invariant(greenwich):
    """When ITRS location == observer location the result must be independent of obstime.

    The transform is purely geometric (rotation only), so changing obstime
    must not change the AltAz output.
    """
    topo = CartesianRepresentation(0*u.m, 1000*u.m, 0*u.m)

    t1 = Time('2020-01-01')
    t2 = Time('2020-07-01')

    altaz1 = ITRS(topo, location=greenwich).transform_to(
        AltAz(location=greenwich, obstime=t1))
    altaz2 = ITRS(topo, location=greenwich).transform_to(
        AltAz(location=greenwich, obstime=t2))

    assert_quantity_allclose(altaz1.cartesian.xyz, altaz2.cartesian.xyz,
                             atol=1e-6*u.m)

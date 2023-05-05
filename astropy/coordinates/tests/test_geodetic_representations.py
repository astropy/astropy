# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test geodetic representations"""
import pytest

from astropy import units as u
from astropy.coordinates.geodetic import (
    BaseGeodeticRepresentation,
    GRS80GeodeticRepresentation,
    WGS72GeodeticRepresentation,
    WGS84GeodeticRepresentation,
)

from astropy.coordinates.representation import (
    CartesianRepresentation,
    REPRESENTATION_CLASSES,
)
from astropy.units import allclose as quantity_allclose
from astropy.units.tests.test_quantity_erfa_ufuncs import vvd

# Preserve the original REPRESENTATION_CLASSES dict so that importing
#   the test file doesn't add a persistent test subclass (LogDRepresentation)
from astropy.coordinates.tests.test_representation import (
    setup_function,
    teardown_function,
)


def test_custombasegeodetic_error():
    # Test incomplete initialization.

    msg = "module 'erfa' has no attribute 'foo'"
    with pytest.raises(AttributeError, match=msg):

        class CustomGeodeticEllipsoidError(BaseGeodeticRepresentation):
            _ellipsoid = "foo"

    assert "customgeodeticellipsoiderror" not in REPRESENTATION_CLASSES

    msg = "CustomGeodeticError requires '_ellipsoid' or '_equatorial_radius' and '_flattening'."
    with pytest.raises(AttributeError, match=msg):

        class CustomGeodeticError(BaseGeodeticRepresentation):
            _flattening = 0.075 * u.dimensionless_unscaled

    assert "customgeodeticerror" not in REPRESENTATION_CLASSES


def test_custombasegeodetic_roundtrip():
    # Test array-valued input in the process.

    class CustomGeodetic(BaseGeodeticRepresentation):
        _flattening = 0.075 * u.dimensionless_unscaled
        _equatorial_radius = 3000000.0 * u.m

    # Test cartesian initialization round-trip
    s1c = CartesianRepresentation(
        x=[1, 3000.0] * u.km, y=[7000.0, 4.0] * u.km, z=[5.0, 6000.0] * u.km
    )

    s2c = CustomGeodetic.from_representation(s1c)

    s3c = CartesianRepresentation.from_representation(s2c)

    assert quantity_allclose(s1c.x, s3c.x)
    assert quantity_allclose(s1c.y, s3c.y)
    assert quantity_allclose(s1c.z, s3c.z)

    # Test geodetic initialization round-trip
    s1g = CustomGeodetic(
        lon=[0.8, 1.3] * u.radian,
        lat=[0.3, 0.98] * u.radian,
        height=[100.0, 367.0] * u.m,
    )

    s2g = CartesianRepresentation.from_representation(s1g)

    s3g = CustomGeodetic.from_representation(s2g)

    assert quantity_allclose(s1g.lon, s3g.lon)
    assert quantity_allclose(s1g.lat, s3g.lat)
    assert quantity_allclose(s1g.height, s3g.height)


def test_wgs84geodetic_roundtrip():
    # Test array-valued input in the process.

    # Test cartesian initialization round-trip
    s1c = CartesianRepresentation(
        x=[1, 3000.0] * u.km, y=[7000.0, 4.0] * u.km, z=[5.0, 6000.0] * u.km
    )

    s2c = WGS84GeodeticRepresentation.from_representation(s1c)

    s3c = CartesianRepresentation.from_representation(s2c)

    assert quantity_allclose(s1c.x, s3c.x)
    assert quantity_allclose(s1c.y, s3c.y)
    assert quantity_allclose(s1c.z, s3c.z)

    # Test geodetic initialization round-trip
    s1g = WGS84GeodeticRepresentation(
        lon=[0.8, 1.3] * u.radian,
        lat=[0.3, 0.98] * u.radian,
        height=[100.0, 367.0] * u.m,
    )

    s2g = CartesianRepresentation.from_representation(s1g)

    s3g = WGS84GeodeticRepresentation.from_representation(s2g)

    assert quantity_allclose(s1g.lon, s3g.lon)
    assert quantity_allclose(s1g.lat, s3g.lat)
    assert quantity_allclose(s1g.height, s3g.height)


def test_geocentric_to_geodetic():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gc2gd"""
    # Here, test the chain.  Direct conversion from Cartesian to
    # various Geodetic representations is done indirectly in test_earth.
    x, y, z = (2e6, 3e6, 5.244e6)

    status = 0  # help for copy & paste of vvd

    gc = CartesianRepresentation(x, y, z, u.m)
    gd = WGS84GeodeticRepresentation.from_cartesian(gc)
    e, p, h = gd.lon.to(u.radian), gd.lat.to(u.radian), gd.height.to(u.m)
    vvd(e, 0.9827937232473290680, 1e-14, "eraGc2gd", "e1", status)
    vvd(p, 0.97160184819075459, 1e-14, "eraGc2gd", "p1", status)
    vvd(h, 331.4172461426059892, 1e-8, "eraGc2gd", "h1", status)

    gd = gd.represent_as(GRS80GeodeticRepresentation)
    e, p, h = gd.lon.to(u.radian), gd.lat.to(u.radian), gd.height.to(u.m)
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e2", status)
    vvd(p, 0.97160184820607853, 1e-14, "eraGc2gd", "p2", status)
    vvd(h, 331.41731754844348, 1e-8, "eraGc2gd", "h2", status)

    gd = gd.represent_as(WGS72GeodeticRepresentation)
    e, p, h = gd.lon.to(u.radian), gd.lat.to(u.radian), gd.height.to(u.m)
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e3", status)
    vvd(p, 0.97160181811015119, 1e-14, "eraGc2gd", "p3", status)
    vvd(h, 333.27707261303181, 1e-8, "eraGc2gd", "h3", status)


def test_geodetic_to_geocentric():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gd2gc"""
    # These tests are also done implicitly in test_earth.py.
    e = 3.1 * u.rad
    p = -0.5 * u.rad
    h = 2500.0 * u.m

    status = 0  # help for copy & paste of vvd

    gd = WGS84GeodeticRepresentation(e, p, h)
    xyz = gd.to_cartesian().get_xyz()
    vvd(xyz[0], -5599000.5577049947, 1e-7, "eraGd2gc", "0/1", status)
    vvd(xyz[1], 233011.67223479203, 1e-7, "eraGd2gc", "1/1", status)
    vvd(xyz[2], -3040909.4706983363, 1e-7, "eraGd2gc", "2/1", status)

    gd = GRS80GeodeticRepresentation(e, p, h)
    xyz = gd.to_cartesian().get_xyz()
    vvd(xyz[0], -5599000.5577260984, 1e-7, "eraGd2gc", "0/2", status)
    vvd(xyz[1], 233011.6722356703, 1e-7, "eraGd2gc", "1/2", status)
    vvd(xyz[2], -3040909.4706095476, 1e-7, "eraGd2gc", "2/2", status)

    gd = WGS72GeodeticRepresentation(e, p, h)
    xyz = gd.to_cartesian().get_xyz()
    vvd(xyz[0], -5598998.7626301490, 1e-7, "eraGd2gc", "0/3", status)
    vvd(xyz[1], 233011.5975297822, 1e-7, "eraGd2gc", "1/3", status)
    vvd(xyz[2], -3040908.6861467111, 1e-7, "eraGd2gc", "2/3", status)


def test_default_height_is_zero():
    gd = WGS84GeodeticRepresentation(10 * u.deg, 20 * u.deg)
    assert gd.lon == 10 * u.deg
    assert gd.lat == 20 * u.deg
    assert gd.height == 0 * u.m


def test_non_angle_error():
    with pytest.raises(u.UnitTypeError):
        WGS84GeodeticRepresentation(20 * u.m, 20 * u.deg, 20 * u.m)


def test_non_length_error():
    with pytest.raises(u.UnitTypeError, match="units of length"):
        WGS84GeodeticRepresentation(10 * u.deg, 20 * u.deg, 30)

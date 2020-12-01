# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test geodetic representations"""

from astropy.coordinates.representation import CartesianRepresentation

from astropy.coordinates.earth import WGS84GeodeticRepresentation
from astropy.units import allclose as quantity_allclose
from astropy import units as u


def test_cartesian_wgs84geodetic_roundtrip():

    s1 = CartesianRepresentation(x=[1, 3000.] * u.km,
                                 y=[7000., 4.] * u.km,
                                 z=[5., 6000.] * u.km)

    s2 = WGS84GeodeticRepresentation.from_representation(s1)

    s3 = CartesianRepresentation.from_representation(s2)

    s4 = WGS84GeodeticRepresentation.from_representation(s3)

    assert quantity_allclose(s1.x, s3.x)
    assert quantity_allclose(s1.y, s3.y)
    assert quantity_allclose(s1.z, s3.z)

    assert quantity_allclose(s2.lon, s4.lon)
    assert quantity_allclose(s2.lat, s4.lat)
    assert quantity_allclose(s2.height, s4.height)


def test_geocentric_to_geodetic():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gc2gd"""
    gc = CartesianRepresentation(2e6, 3e6, 5.244e6, u.m)
    gd = WGS84GeodeticRepresentation.from_cartesian(gc)

    assert quantity_allclose(gd.lon, 0.98279372324732907 * u.rad, atol=1e-14 * u.rad)
    assert quantity_allclose(gd.lat, 0.97160184820607853 * u.rad, atol=1e-14 * u.rad)
    assert quantity_allclose(gd.height, 331.41731754844348 * u.m, rtol=1e-5, atol=1e-8 * u.m)


def test_geodetic_to_geocentric():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gd2gc"""
    gd = WGS84GeodeticRepresentation(3.1 * u.rad, -0.5 * u.rad, 2500.0 * u.m)
    gc = gd.to_cartesian()

    assert quantity_allclose(gc.x, -5599000.5577049947 * u.m, atol=1e-7 * u.m)
    assert quantity_allclose(gc.y, 233011.67223479203 * u.m, atol=1e-7 * u.m)
    assert quantity_allclose(gc.z, -3040909.4706983363 * u.m, atol=1e-7 * u.m)

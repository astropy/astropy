# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

asdf = pytest.importorskip('asdf')
from asdf.tests.helpers import assert_roundtrip_tree

from astropy import units as u
from astropy.coordinates.angles import Longitude, Latitude
from astropy.coordinates.earth import EarthLocation, ELLIPSOIDS


@pytest.fixture
def position():

    lon = Longitude([0., 45., 90., 135., 180., -180, -90, -45], u.deg,
                    wrap_angle=180*u.deg)
    lat = Latitude([+0., 30., 60., +90., -90., -60., -30., 0.], u.deg)
    h = u.Quantity([0.1, 0.5, 1.0, -0.5, -1.0, +4.2, -11., -.1], u.m)

    return lon, lat, h


def test_earthlocation_quantity(tmpdir):

    location = EarthLocation(lat=34.4900*u.deg, lon=-104.221800*u.deg,
                             height=40*u.km)

    tree = dict(location=location)
    assert_roundtrip_tree(tree, tmpdir)


def test_earthlocation(position, tmpdir):

    x, y, z = EarthLocation.from_geodetic(*position).to_geocentric()
    geocentric = EarthLocation(x, y, z)

    tree = dict(location=geocentric)
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.parametrize('ellipsoid', ELLIPSOIDS)
def test_earthlocation_geodetic(position, ellipsoid, tmpdir):

    location = EarthLocation.from_geodetic(*position, ellipsoid=ellipsoid)

    tree = dict(location=location)
    assert_roundtrip_tree(tree, tmpdir)


def test_earthlocation_site(tmpdir):
    orig_sites = getattr(EarthLocation, '_site_registry', None)
    try:
        EarthLocation._get_site_registry(force_builtin=True)
        rog = EarthLocation.of_site('greenwich')
        tree = dict(location=rog)
        assert_roundtrip_tree(tree, tmpdir)
    finally:
        EarthLocation._site_registry = orig_sites

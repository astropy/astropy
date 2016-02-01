# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

"""Test initialization of angles not already covered by the API tests"""

import pickle
import numpy as np

from ..earth import EarthLocation, ELLIPSOIDS
from ..angles import Longitude, Latitude
from ...tests.helper import pytest, quantity_allclose
from ... import units as u

def allclose_m14(a, b, rtol=1.e-14, atol=None):
    if atol is None:
        atol = 1.e-14 * getattr(a, 'unit', 1)
    return quantity_allclose(a, b, rtol, atol)

def allclose_m8(a, b, rtol=1.e-8, atol=None):
    if atol is None:
        atol = 1.e-8 * getattr(a, 'unit', 1)
    return quantity_allclose(a, b, rtol, atol)


def isclose_m14(val, ref):
    return np.array([allclose_m14(v, r) for (v, r) in zip(val, ref)])

def isclose_m8(val, ref):
    return np.array([allclose_m8(v, r) for (v, r) in zip(val, ref)])


def vvd(val, valok, dval, func, test, status):
    """Mimic routine of erfa/src/t_erfa_c.c (to help copy & paste)"""
    assert quantity_allclose(val, valok * val.unit, atol=dval * val.unit)


def test_gc2gd():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gc2gd"""
    x, y, z = (2e6, 3e6, 5.244e6)

    status = 0  # help for copy & paste of vvd

    location = EarthLocation.from_geocentric(x, y, z, u.m)
    e, p, h = location.to_geodetic('WGS84')
    e, p, h = e.to(u.radian), p.to(u.radian), h.to(u.m)
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e2", status)
    vvd(p, 0.97160184820607853, 1e-14, "eraGc2gd", "p2", status)
    vvd(h, 331.41731754844348, 1e-8, "eraGc2gd", "h2", status)

    e, p, h = location.to_geodetic('GRS80')
    e, p, h = e.to(u.radian), p.to(u.radian), h.to(u.m)
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e2", status)
    vvd(p, 0.97160184820607853, 1e-14, "eraGc2gd", "p2", status)
    vvd(h, 331.41731754844348, 1e-8, "eraGc2gd", "h2", status)

    e, p, h = location.to_geodetic('WGS72')
    e, p, h = e.to(u.radian), p.to(u.radian), h.to(u.m)
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e3", status)
    vvd(p, 0.97160181811015119, 1e-14, "eraGc2gd", "p3", status)
    vvd(h, 333.27707261303181, 1e-8, "eraGc2gd", "h3", status)


def test_gd2gc():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gd2gc"""
    e = 3.1 * u.rad
    p = -0.5 * u.rad
    h = 2500.0 * u.m

    status = 0  # help for copy & paste of vvd

    location = EarthLocation.from_geodetic(e, p, h, ellipsoid='WGS84')
    xyz = tuple(v.to(u.m) for v in location.to_geocentric())
    vvd(xyz[0], -5599000.5577049947, 1e-7, "eraGd2gc", "0/1", status)
    vvd(xyz[1], 233011.67223479203, 1e-7, "eraGd2gc", "1/1", status)
    vvd(xyz[2], -3040909.4706983363, 1e-7, "eraGd2gc", "2/1", status)

    location = EarthLocation.from_geodetic(e, p, h, ellipsoid='GRS80')
    xyz = tuple(v.to(u.m) for v in location.to_geocentric())
    vvd(xyz[0], -5599000.5577260984, 1e-7, "eraGd2gc", "0/2", status)
    vvd(xyz[1], 233011.6722356703, 1e-7, "eraGd2gc", "1/2", status)
    vvd(xyz[2], -3040909.4706095476, 1e-7, "eraGd2gc", "2/2", status)

    location = EarthLocation.from_geodetic(e, p, h, ellipsoid='WGS72')
    xyz = tuple(v.to(u.m) for v in location.to_geocentric())
    vvd(xyz[0], -5598998.7626301490, 1e-7, "eraGd2gc", "0/3", status)
    vvd(xyz[1], 233011.5975297822, 1e-7, "eraGd2gc", "1/3", status)
    vvd(xyz[2], -3040908.6861467111, 1e-7, "eraGd2gc", "2/3", status)


class TestInput():
    def setup(self):
        self.lon = Longitude([0., 45., 90., 135., 180., -180, -90, -45], u.deg,
                             wrap_angle=180*u.deg)
        self.lat = Latitude([+0., 30., 60., +90., -90., -60., -30., 0.], u.deg)
        self.h = u.Quantity([0.1, 0.5, 1.0, -0.5, -1.0, +4.2, -11.,-.1], u.m)
        self.location = EarthLocation.from_geodetic(self.lon, self.lat, self.h)
        self.x, self.y, self.z = self.location.to_geocentric()

    def test_default_ellipsoid(self):
        assert self.location.ellipsoid == EarthLocation._ellipsoid

    def test_geo_attributes(self):
        assert all([np.all(_1 == _2)
                    for _1, _2 in zip(self.location.geodetic,
                                      self.location.to_geodetic())])
        assert all([np.all(_1 == _2)
                    for _1, _2 in zip(self.location.geocentric,
                                      self.location.to_geocentric())])

    def test_attribute_classes(self):
        """Test that attribute classes are correct (and not EarthLocation)"""
        assert type(self.location.x) is u.Quantity
        assert type(self.location.y) is u.Quantity
        assert type(self.location.z) is u.Quantity
        assert type(self.location.longitude) is Longitude
        assert type(self.location.latitude) is Latitude
        assert type(self.location.height) is u.Quantity

    def test_input(self):
        """Check input is parsed correctly"""

        # units of length should be assumed geocentric
        geocentric = EarthLocation(self.x, self.y, self.z)
        assert np.all(geocentric == self.location)
        geocentric2 = EarthLocation(self.x.value, self.y.value, self.z.value,
                                    self.x.unit)
        assert np.all(geocentric2 == self.location)
        geodetic = EarthLocation(self.lon, self.lat, self.h)
        assert np.all(geodetic == self.location)
        geodetic2 = EarthLocation(self.lon.to(u.degree).value,
                                  self.lat.to(u.degree).value,
                                  self.h.to(u.m).value)
        assert np.all(geodetic2 == self.location)
        geodetic3 = EarthLocation(self.lon, self.lat)
        assert allclose_m14(geodetic3.longitude.value,
                            self.location.longitude.value)
        assert allclose_m14(geodetic3.latitude.value,
                            self.location.latitude.value)
        assert not np.any(isclose_m14(geodetic3.height.value,
                                      self.location.height.value))
        geodetic4 = EarthLocation(self.lon, self.lat, self.h[-1])
        assert allclose_m14(geodetic4.longitude.value,
                            self.location.longitude.value)
        assert allclose_m14(geodetic4.latitude.value,
                            self.location.latitude.value)
        assert allclose_m14(geodetic4.height[-1].value,
                            self.location.height[-1].value)
        assert not np.any(isclose_m14(geodetic4.height[:-1].value,
                                      self.location.height[:-1].value))
        # check length unit preservation
        geocentric5 = EarthLocation(self.x, self.y, self.z, u.pc)
        assert geocentric5.unit is u.pc
        assert geocentric5.x.unit is u.pc
        assert geocentric5.height.unit is u.pc
        assert allclose_m14(geocentric5.x.to(self.x.unit).value, self.x.value)
        geodetic5 = EarthLocation(self.lon, self.lat, self.h.to(u.pc))
        assert geodetic5.unit is u.pc
        assert geodetic5.x.unit is u.pc
        assert geodetic5.height.unit is u.pc
        assert allclose_m14(geodetic5.x.to(self.x.unit).value, self.x.value)

    def test_invalid_input(self):
        """Check invalid input raises exception"""
        # incomprehensible by either raises TypeError
        with pytest.raises(TypeError):
            EarthLocation(self.lon, self.y, self.z)

        # wrong units
        with pytest.raises(u.UnitsError):
            EarthLocation.from_geocentric(self.lon, self.lat, self.lat)
        # inconsistent units
        with pytest.raises(u.UnitsError):
            EarthLocation.from_geocentric(self.h, self.lon, self.lat)
        # floats without a unit
        with pytest.raises(TypeError):
            EarthLocation.from_geocentric(self.x.value, self.y.value,
                                          self.z.value)
        # inconsistent shape
        with pytest.raises(ValueError):
            EarthLocation.from_geocentric(self.x, self.y, self.z[:5])

        # inconsistent units
        with pytest.raises(u.UnitsError):
            EarthLocation.from_geodetic(self.x, self.y, self.z)
        # inconsistent shape
        with pytest.raises(ValueError):
            EarthLocation.from_geodetic(self.lon, self.lat, self.h[:5])

    def test_slicing(self):
        # test on WGS72 location, so we can check the ellipsoid is passed on
        locwgs72 = EarthLocation.from_geodetic(self.lon, self.lat, self.h,
                                               ellipsoid='WGS72')
        loc_slice1 = locwgs72[4]
        assert isinstance(loc_slice1, EarthLocation)
        assert loc_slice1.unit is locwgs72.unit
        assert loc_slice1.ellipsoid == locwgs72.ellipsoid == 'WGS72'
        assert not loc_slice1.shape
        with pytest.raises(TypeError):
            loc_slice1[0]
        with pytest.raises(IndexError):
            len(loc_slice1)

        loc_slice2 = locwgs72[4:6]
        assert isinstance(loc_slice2, EarthLocation)
        assert len(loc_slice2) == 2
        assert loc_slice2.unit is locwgs72.unit
        assert loc_slice2.ellipsoid == locwgs72.ellipsoid
        assert loc_slice2.shape == (2,)
        loc_x = locwgs72['x']
        assert type(loc_x) is u.Quantity
        assert loc_x.shape == locwgs72.shape
        assert loc_x.unit is locwgs72.unit

    def test_invalid_ellipsoid(self):
        # unknown ellipsoid
        with pytest.raises(ValueError):
            EarthLocation.from_geodetic(self.lon, self.lat, self.h,
                                        ellipsoid='foo')
        with pytest.raises(TypeError):
            EarthLocation(self.lon, self.lat, self.h, ellipsoid='foo')

        with pytest.raises(ValueError):
            self.location.ellipsoid = 'foo'

        with pytest.raises(ValueError):
            self.location.to_geodetic('foo')

    @pytest.mark.parametrize('ellipsoid', ELLIPSOIDS)
    def test_ellipsoid(self, ellipsoid):
        """Test that different ellipsoids are understood, and differ"""
        # check that heights differ for different ellipsoids
        # need different tolerance, since heights are relative to ~6000 km
        lon, lat, h = self.location.to_geodetic(ellipsoid)
        if ellipsoid == self.location.ellipsoid:
            assert allclose_m8(h.value, self.h.value)
        else:
            # Some heights are very similar for some; some lon, lat identical.
            assert not np.all(isclose_m8(h.value, self.h.value))

        # given lon, lat, height, check that x,y,z differ
        location = EarthLocation.from_geodetic(self.lon, self.lat, self.h,
                                               ellipsoid=ellipsoid)
        if ellipsoid == self.location.ellipsoid:
            assert allclose_m14(location.z.value, self.z.value)
        else:
            assert not np.all(isclose_m14(location.z.value, self.z.value))


def test_pickling():
    """Regression test against #4304."""
    el = EarthLocation(0.*u.m, 6000*u.km, 6000*u.km)
    s = pickle.dumps(el)
    el2 = pickle.loads(s)
    assert el == el2

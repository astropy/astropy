# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

"""Test initalization of angles not already covered by the API tests"""

import numpy as np
from numpy.testing.utils import assert_allclose

from ..earth import EarthLocation
from ..angles import Longitude, Latitude
from ...tests.helper import pytest
from ... import units as u


def vvd(val, valok, dval, func, test, status):
    """Mimic routine of erfa/src/t_erfa_c.c (to help copy & paste)"""
    assert np.allclose(val, valok, atol=dval)


def test_gc2gd():
    """Test that we reproduce erfa/src/t_erfa_c.c t_gc2gd"""
    x, y, z = (2e6, 3e6, 5.244e6)

    status = 0  # help for copy & paste of vvd

    location = EarthLocation.from_geocentric(x, y, z, u.m)
    e, p, h = location.to_geodetic('WGS84')
    e, p, h = e.to(u.radian).value, p.to(u.radian).value, h.to(u.m).value
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e2", status)
    vvd(p, 0.97160184820607853, 1e-14, "eraGc2gd", "p2", status)
    vvd(h, 331.41731754844348, 1e-8, "eraGc2gd", "h2", status)

    e, p, h = location.to_geodetic('GRS80')
    e, p, h = e.to(u.radian).value, p.to(u.radian).value, h.to(u.m).value
    vvd(e, 0.98279372324732907, 1e-14, "eraGc2gd", "e2", status)
    vvd(p, 0.97160184820607853, 1e-14, "eraGc2gd", "p2", status)
    vvd(h, 331.41731754844348, 1e-8, "eraGc2gd", "h2", status)

    e, p, h = location.to_geodetic('WGS72')
    e, p, h = e.to(u.radian).value, p.to(u.radian).value, h.to(u.m).value
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
    xyz = location.to_geocentric()
    vvd(xyz[0], -5599000.5577049947, 1e-7, "eraGd2gc", "0/1", status)
    vvd(xyz[1], 233011.67223479203, 1e-7, "eraGd2gc", "1/1", status)
    vvd(xyz[2], -3040909.4706983363, 1e-7, "eraGd2gc", "2/1", status)

    location = EarthLocation.from_geodetic(e, p, h, ellipsoid='GRS80')
    xyz = location.to_geocentric()
    vvd(xyz[0], -5599000.5577260984, 1e-7, "eraGd2gc", "0/2", status)
    vvd(xyz[1], 233011.6722356703, 1e-7, "eraGd2gc", "1/2", status)
    vvd(xyz[2], -3040909.4706095476, 1e-7, "eraGd2gc", "2/2", status)

    location = EarthLocation.from_geodetic(e, p, h, ellipsoid='WGS72')
    xyz = location.to_geocentric()
    vvd(xyz[0], -5598998.7626301490, 1e-7, "eraGd2gc", "0/3", status)
    vvd(xyz[1], 233011.5975297822, 1e-7, "eraGd2gc", "1/3", status)
    vvd(xyz[2], -3040908.6861467111, 1e-7, "eraGd2gc", "2/3", status)


class TestInput():
    def setup(self):
        self.lon = Longitude([0., 45., 90., 135., 180., -180, -90, -45], u.deg)
        self.lat = Latitude([+0., 30., 60., +90., -90., -60., -30., 0.], u.deg)
        self.h = u.Quantity([0.1, 0.5, 1.0, -0.5, -1.0, +4.2, -11.,-.1], u.m)
        self.location = EarthLocation.from_geodetic(self.lon, self.lat, self.h)
        self.x, self.y, self.z = self.location.to_geocentric()

    def test_input_validation(self):
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
        assert_allclose(geodetic3.longitude, self.location.longitude)
        assert_allclose(geodetic3.latitude, self.location.latitude)
        assert not np.any(np.isclose(geodetic3.height.value,
                                     self.location.height.value))
        geodetic4 = EarthLocation(self.lon, self.lat, self.h[-1])
        assert_allclose(geodetic4.longitude, self.location.longitude)
        assert_allclose(geodetic4.latitude, self.location.latitude)
        assert_allclose(geodetic4.height[-1], self.location.height[-1])
        assert not np.any(np.isclose(geodetic4.height[:-1].value,
                                     self.location.height[:-1].value))

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

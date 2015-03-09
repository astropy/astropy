# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .. import models
from numpy.testing import utils
from ...tests.helper import pytest


def test_RotateNative2Celestial():
    lon, lat, lon_pole = 42, 43, 44
    model = models.RotateNative2Celestial(lon, lat, lon_pole)
    model.lon = model.lon + 1
    model.lon_pole = model.lon_pole + 1
    model.lat = model.lat + 1
    assert model.lon == lon + 1
    assert model.lat == lat + 1
    assert model.lon_pole == lon_pole + 1


def test_native_celestial_native():
    lon, lat, lon_pole = 42, 43, 44
    n2c = models.RotateNative2Celestial(lon, lat, lon_pole)
    c2n = models.RotateCelestial2Native(lon, lat, lon_pole)

    nnphi, ntheta = 33, 44
    calpha, cdelta = n2c(nnphi, ntheta)
    nnphi2, ntheta2 = c2n(calpha, cdelta)
    utils.assert_allclose(nnphi2, nnphi)
    utils.assert_allclose(ntheta2, ntheta)

    assert n2c.inverse(nnphi, ntheta) == c2n(nnphi, ntheta)
    assert c2n.inverse(nnphi, ntheta) == n2c(nnphi, ntheta)


def test_native_celestial_lat90():
    n2c = models.RotateNative2Celestial(1, 90, 0)
    alpha, delta = n2c(1, 1)
    utils.assert_allclose( delta, 1)
    utils.assert_allclose(alpha, 182)


def test_Rotation2D():
    model = models.Rotation2D(angle=90)
    x, y = model(1, 0)
    utils.assert_allclose([x, y], [0, 1], atol=1e-10)


def test_Rotation2D_inverse():
    model = models.Rotation2D(angle=234.23494)
    x, y = model.inverse(*model(1, 0))
    utils.assert_allclose([x, y], [1, 0], atol=1e-10)

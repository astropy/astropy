# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .. import models
from numpy.testing import utils
from ...tests.helper import pytest


def test_RotateNative2Celestial():
    phi, theta, psi = 42, 43, 44
    model = models.RotateNative2Celestial(phi, theta, psi)
    model.phi = model.phi + 1
    model.psi = model.psi + 1
    model.theta = model.theta + 1
    assert model.phi == phi + 1
    assert model.theta == theta + 1
    assert model.psi == psi + 1


def test_native_celestial_native():
    phi, theta, psi = 42, 43, 44
    n2c = models.RotateNative2Celestial(phi, theta, psi)
    c2n = models.RotateCelestial2Native(phi, theta, psi)

    nnphi, ntheta = 33, 44
    calpha, cdelta = n2c(nnphi, ntheta)
    nnphi2, ntheta2 = c2n(calpha, cdelta)
    utils.assert_allclose(nnphi2, nnphi)
    utils.assert_allclose(ntheta2, ntheta)

    assert n2c.inverse()(nnphi, ntheta) == c2n(nnphi, ntheta)
    assert c2n.inverse()(nnphi, ntheta) == n2c(nnphi, ntheta)


def test_native_celestial_theta90():
    n2c = models.RotateNative2Celestial(1, 90, 0)
    alpha, delta = n2c(1, 1)
    utils.assert_allclose( delta, 1)
    utils.assert_allclose(alpha, 182)


def test_MatrixRotation2D():
    model = models.MatrixRotation2D(angle=90)
    x, y = model(1, 0)
    utils.assert_allclose([x, y], [0, -1], atol=1e-10)

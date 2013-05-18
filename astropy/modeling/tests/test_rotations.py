# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .. import rotations
from numpy.testing.utils import assert_allclose
from ...tests.helper import pytest

def test_RotateNative2Celestial():
    phi, theta, psi = 42, 43, 44
    model = rotations.RotateNative2Celestial(phi, theta, psi)
    model.phi = model.phi + 1
    model.psi = model.psi + 1
    model.theta = model.theta + 1
    assert model.phi == phi + 1
    assert model.theta == theta + 1
    assert model.psi == psi + 1


def test_native_celestial_native():
    phi, theta, psi = 42, 43, 44
    n2c = rotations.RotateNative2Celestial(phi, theta, psi)
    c2n = rotations.RotateCelestial2Native(phi, theta, psi)

    nphi, ntheta = 33, 44
    calpha, cdelta = n2c(nphi, ntheta)
    nphi2, ntheta2 = c2n(calpha, cdelta)
    assert_allclose(nphi2, nphi)
    assert_allclose(ntheta2, ntheta)

    assert n2c.inverse(phi, theta, psi)(nphi, ntheta) == c2n(nphi, ntheta)
    assert c2n.inverse(phi, theta, psi)(nphi, ntheta) == n2c(nphi, ntheta)

def test_MatrixRotation2D():
    model = rotations.MatrixRotation2D(angle=90)
    x, y = model(1, 0)
    assert_allclose([x, y], [0, -1], atol=1e-10)

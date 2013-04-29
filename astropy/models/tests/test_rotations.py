# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ..rotations import (RotateNative2Celestial, RotateCelestial2Native,
                         MatrixRotation2D)
from ...tests.compat import assert_allclose

def test_native_celestial_native():
    phi, theta, psi = 42, 43, 44
    n2c = RotateNative2Celestial(phi, theta, psi)
    c2n = RotateCelestial2Native(phi, theta, psi)

    nphi, ntheta = 33, 44
    calpha, cdelta = n2c(nphi, ntheta)
    nphi2, ntheta2 = c2n(calpha, cdelta)
    assert_allclose(nphi2, nphi)
    assert_allclose(ntheta2, ntheta)


def test_MatrixRotation2D():
    angle = 42
    # TODO: this fails
    # m = MatrixRotation2D(angle=angle)
    # result = m(3, 4)
    # print result
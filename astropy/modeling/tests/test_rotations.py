# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from math import cos, sin
import numpy as np
from numpy.testing import utils
from .. import models

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


def test_euler_angle_rotations():
    x = (0, 0)
    y = (90, 0)
    z = (0, 90)
    negx = (180, 0)
    negy = (-90, 0)

    # rotate y into minus z
    model = models.EulerAngleRotation(0, 90, 0, 'zxz')
    utils.assert_allclose(model(*z), y, atol=10**-12)
    # rotate z into minus x
    model = models.EulerAngleRotation(0, 90, 0, 'zyz')
    utils.assert_allclose(model(*z), negx, atol=10**-12)
    # rotate x into minus y
    model = models.EulerAngleRotation(0, 90, 0, 'yzy')
    utils.assert_allclose(model(*x), negy, atol=10**-12)


def test_euler_against_sky():
    """
    Test sky transformations against euler angle rotations
    """
    # Write the Euler angles in terms of FITS WCS angles on the sky

    # This corresponds to FITS WCS CRVAL1 = 5.63 deg
    phi = 5.63 + 90
    # This corresponds to FITS WCS CRVAL2 = -72.63 deg
    theta = 90 - -72.63
    # This corresponds to FITS WCS LONPOLE = 180 deg
    psi = 180 - 90

    # Create rotations between celestial and Native systems
    c2n = models.RotateCelestial2Native(5.63, -72.63, 180)
    n2c = models.RotateNative2Celestial(5.63, -72.63, 180)

    # And the corresponding rotation in terms of Euler angles
    # This corresponds to c2n since it's rotating coordinate systems,
    # not vectors
    model= models.EulerAngleRotation(phi, theta, -psi, 'zxz')

    utils.assert_allclose(model(1, 23.1), c2n(1, 23.1))
    ra, dec = model.inverse(1, 23.1)
    if ra < 0:
        ra += 360
    utils.assert_allclose((ra, dec), n2c(1, 23.1))



euler_axes_order = ['zxz', 'zyz', 'yzy', 'yxy', 'xyx', 'xzx']

@pytest.mark.parametrize(('axes_order'), euler_axes_order)
def test_euler_angles(axes_order):
    """
    Tests against all Euler sequences.
    The rotation matrices definitions come from Wikipedia.
    """
    phi = np.deg2rad(23.4)
    theta = np.deg2rad(12.2)
    psi = np.deg2rad(34)
    c1 =cos(phi)
    c2 = cos(theta)
    c3 = cos(psi)
    s1 = sin(phi)
    s2 = sin(theta)
    s3 = sin(psi)

    matrices = {'zxz': np.array([[(c1*c3 - c2*s1*s3), (-c1*s3 - c2*c3*s1), (s1*s2)],
                        [(c3*s1 + c1*c2*s3), (c1*c2*c3 - s1*s3), (-c1*s2)],
                        [(s2*s3),             (c3*s2),            (c2)   ]]),
                'zyz': np.array([[(c1*c2*c3 -s1*s3), (-c3*s1 - c1*c2*s3), (c1*s2)],
                        [(c1*s3 +c2*c3*s1), (c1*c3 - c2*s1*s3),  (s1*s2)],
                        [(-c3*s2),           (s2*s3),             (c2)]]),
                'yzy': np.array([[(c1*c2*c3 - s1*s3), (-c1*s2), (c3*s1+c1*c2*s3)],
                                  [(c3*s2), (c2), (s2*s3)],
                                  [(-c1*s3 - c2*c3*s1), (s1*s2), (c1*c3-c2*s1*s3)]]),
                'yxy': np.array([[(c1*c3 - c2*s1*s3), (s1*s2), (c1*s3+c2*c3*s1)],
                                 [(s2*s3), (c2), (-c3*s2)],
                                 [(-c3*s1 - c1*c2*s3), (c1*s2), (c1*c2*c3 - s1*s3)]]),
                'xyx': np.array([[(c2), (s2*s3), (c3*s2)],
                                 [(s1*s2), (c1*c3 - c2*s1*s3), (-c1*s3 - c2*c3*s1)],
                                 [(-c1*s2), (c3*s1 + c1*c2*s3), (c1*c2*c3 - s1*s3)]]),
                'xzx': np.array([[(c2), (-c3*s2), (s2*s3)],
                                 [(c1*s2), (c1*c2*c3 - s1*s3), (-c3*s1 - c1*c2*s3)],
                                 [(s1*s2), (c1*s3 + c2*c3*s1), (c1*c3 - c2*s1*s3)]])

                }
    model = models.EulerAngleRotation(23.4, 12.2, 34, axes_order)
    mat = model._create_matrix(phi, theta, psi, axes_order)


    utils.assert_allclose(mat.T, matrices[axes_order])#get_rotation_matrix(axes_order))


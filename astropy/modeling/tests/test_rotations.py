# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from math import cos, sin

import pytest
import numpy as np
from numpy.testing import assert_allclose

from .. import models
from ...wcs import wcs


@pytest.mark.parametrize(('inp'), [(0, 0), (4000, -20.56), (-2001.5, 45.9),
                                   (0, 90), (0, -90), (np.mgrid[:4, :6])])
def test_against_wcslib(inp):
    w = wcs.WCS()
    crval = [202.4823228, 47.17511893]
    w.wcs.crval = crval
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    lonpole = 180
    tan = models.Pix2Sky_TAN()
    n2c = models.RotateNative2Celestial(crval[0], crval[1], lonpole)
    c2n = models.RotateCelestial2Native(crval[0], crval[1], lonpole)
    m = tan | n2c
    minv = c2n | tan.inverse

    radec = w.wcs_pix2world(inp[0], inp[1], 1)
    xy = w.wcs_world2pix(radec[0], radec[1], 1)

    assert_allclose(m(*inp), radec, atol=1e-12)
    assert_allclose(minv(*radec), xy, atol=1e-12)


@pytest.mark.parametrize(('inp'), [(0, 0), (40, -20.56), (21.5, 45.9)])
def test_roundtrip_sky_rotaion(inp):
    lon, lat, lon_pole = 42, 43, 44
    n2c = models.RotateNative2Celestial(lon, lat, lon_pole)
    c2n = models.RotateCelestial2Native(lon, lat, lon_pole)
    assert_allclose(n2c.inverse(*n2c(*inp)), inp, atol=1e-13)
    assert_allclose(c2n.inverse(*c2n(*inp)), inp, atol=1e-13)


def test_native_celestial_lat90():
    n2c = models.RotateNative2Celestial(1, 90, 0)
    alpha, delta = n2c(1, 1)
    assert_allclose(delta, 1)
    assert_allclose(alpha, 182)


def test_Rotation2D():
    model = models.Rotation2D(angle=90)
    x, y = model(1, 0)
    assert_allclose([x, y], [0, 1], atol=1e-10)


def test_Rotation2D_inverse():
    model = models.Rotation2D(angle=234.23494)
    x, y = model.inverse(*model(1, 0))
    assert_allclose([x, y], [1, 0], atol=1e-10)


def test_euler_angle_rotations():
    x = (0, 0)
    y = (90, 0)
    z = (0, 90)
    negx = (180, 0)
    negy = (-90, 0)

    # rotate y into minus z
    model = models.EulerAngleRotation(0, 90, 0, 'zxz')
    assert_allclose(model(*z), y, atol=10**-12)
    # rotate z into minus x
    model = models.EulerAngleRotation(0, 90, 0, 'zyz')
    assert_allclose(model(*z), negx, atol=10**-12)
    # rotate x into minus y
    model = models.EulerAngleRotation(0, 90, 0, 'yzy')
    assert_allclose(model(*x), negy, atol=10**-12)


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
    c1 = cos(phi)
    c2 = cos(theta)
    c3 = cos(psi)
    s1 = sin(phi)
    s2 = sin(theta)
    s3 = sin(psi)

    matrices = {'zxz': np.array([[(c1*c3 - c2*s1*s3), (-c1*s3 - c2*c3*s1), (s1*s2)],
                        [(c3*s1 + c1*c2*s3), (c1*c2*c3 - s1*s3), (-c1*s2)],
                        [(s2*s3), (c3*s2), (c2)]]),
                'zyz': np.array([[(c1*c2*c3 - s1*s3), (-c3*s1 - c1*c2*s3), (c1*s2)],
                        [(c1*s3 + c2*c3*s1), (c1*c3 - c2*s1*s3), (s1*s2)],
                        [(-c3*s2), (s2*s3), (c2)]]),
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

    assert_allclose(mat.T, matrices[axes_order])  # get_rotation_matrix(axes_order))

# Licensed under a 3-clause BSD style license - see LICENSE.rst

# pylint: disable=invalid-name
import unittest.mock as mk
from math import cos, sin

import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.modeling import models, rotations
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import wcs


@pytest.mark.parametrize(('inp'), [(0, 0), (4000, -20.56), (-2001.5, 45.9),
                                   (0, 90), (0, -90), (np.mgrid[:4, :6]),
                                   ([[1, 2, 3], [4,   5,  6]],
                                    [[7, 8, 9], [10, 11, 12]]),
                                   ([[[1,   2,  3,  4],
                                      [5,   6,  7,  8],
                                      [9,  10, 11, 12]],
                                     [[13, 14, 15, 16],
                                      [17, 18, 19, 20],
                                      [21, 22, 23, 24]]],
                                    [[[25, 26, 27, 28],
                                      [29, 30, 31, 32],
                                      [33, 34, 35, 36]],
                                     [[37, 38, 39, 40],
                                      [41, 42, 43, 44],
                                      [45, 46, 47, 48]]])])
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


@pytest.mark.parametrize(('inp'), [(1e-5, 1e-4), (40, -20.56), (21.5, 45.9),
                                   ([[1, 2, 3], [4,   5,  6]],
                                    [[7, 8, 9], [10, 11, 12]]),
                                   ([[[1,   2,  3,  4],
                                      [5,   6,  7,  8],
                                      [9,  10, 11, 12]],
                                     [[13, 14, 15, 16],
                                      [17, 18, 19, 20],
                                      [21, 22, 23, 24]]],
                                    [[[25, 26, 27, 28],
                                      [29, 30, 31, 32],
                                      [33, 34, 35, 36]],
                                     [[37, 38, 39, 40],
                                      [41, 42, 43, 44],
                                      [45, 46, 47, 48]]])])
def test_roundtrip_sky_rotation(inp):
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


def test_Rotation2D_quantity():
    model = models.Rotation2D(angle=90*u.deg)
    x, y = model(1*u.deg, 0*u.arcsec)
    assert_quantity_allclose([x, y], [0, 1]*u.deg, atol=1e-10*u.deg)


def test_Rotation2D_inverse():
    model = models.Rotation2D(angle=234.23494)
    x, y = model.inverse(*model(1, 0))
    assert_allclose([x, y], [1, 0], atol=1e-10)


def test_Rotation2D_errors():
    model = models.Rotation2D(angle=90*u.deg)

    # Bad evaluation input shapes
    x = np.array([1, 2])
    y = np.array([1, 2, 3])
    message = "Expected input arrays to have the same shape"

    with pytest.raises(ValueError) as err:
        model.evaluate(x, y, model.angle)
    assert str(err.value) == message
    with pytest.raises(ValueError) as err:
        model.evaluate(y, x, model.angle)
    assert str(err.value) == message

    # Bad evaluation units
    x = np.array([1, 2])
    y = np.array([1, 2])
    message = "x and y must have compatible units"
    with pytest.raises(u.UnitsError) as err:
        model.evaluate(x * u.m, y, model.angle)
    assert str(err.value) == message


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
    mat = rotations._create_matrix([phi, theta, psi], axes_order)

    assert_allclose(mat.T, matrices[axes_order])  # get_rotation_matrix(axes_order))


def test_rotation_3d():
    """
    A sanity test - when V2_REF = 0 and V3_REF = 0,
    for V2, V3 close to the origin
    ROLL_REF should be approximately PA_V3 .

    (Test taken from JWST SIAF report.)
    """
    def _roll_angle_from_matrix(matrix, v2, v3):
        X = (
            -(matrix[2, 0] * np.cos(v2) + matrix[2, 1] * np.sin(v2)) *
            np.sin(v3) + matrix[2, 2] * np.cos(v3)
        )
        Y = (
            (matrix[0, 0] * matrix[1, 2] - matrix[1, 0] * matrix[0, 2]) * np.cos(v2) +
            (matrix[0, 1] * matrix[1, 2] - matrix[1, 1] * matrix[0, 2]) * np.sin(v2)
        )
        new_roll = np.rad2deg(np.arctan2(Y, X))
        if new_roll < 0:
            new_roll += 360
        return new_roll

    # reference points on sky and in a coordinate frame associated
    # with the telescope
    ra_ref = 165  # in deg
    dec_ref = 54  # in deg
    v2_ref = 0
    v3_ref = 0
    pa_v3 = 37  # in deg

    v2 = np.deg2rad(2.7e-6)  # in deg.01 # in arcsec
    v3 = np.deg2rad(2.7e-6)  # in deg .01 # in arcsec
    angles = [v2_ref, -v3_ref, pa_v3, dec_ref, -ra_ref]
    axes = "zyxyz"
    M = rotations._create_matrix(np.deg2rad(angles) * u.deg, axes)
    roll_angle = _roll_angle_from_matrix(M, v2, v3)
    assert_allclose(roll_angle, pa_v3, atol=1e-3)


def test_spherical_rotation():
    """
    Test taken from JWST INS report - converts
    JWST telescope (V2, V3) coordinates to RA, DEC.
    """
    ra_ref = 165  # in deg
    dec_ref = 54  # in deg
    v2_ref = -503.654472 / 3600  # in deg
    v3_ref = -318.742464 / 3600  # in deg
    r0 = 37  # in deg

    v2 = 210  # in deg
    v3 = -75  # in deg
    expected_ra_dec = (107.12810484789563, -35.97940247128502)  # in deg
    angles = np.array([v2_ref, -v3_ref, r0, dec_ref, -ra_ref])
    axes = "zyxyz"
    v2s = rotations.RotationSequence3D(angles, axes_order=axes)
    x, y, z = rotations.spherical2cartesian(v2, v3)
    x1, y1, z1 = v2s(x, y, z)
    radec = rotations.cartesian2spherical(x1, y1, z1)
    assert_allclose(radec, expected_ra_dec, atol=1e-10)

    v2s = rotations.SphericalRotationSequence(angles, axes_order=axes)
    radec = v2s(v2, v3)
    assert_allclose(radec, expected_ra_dec, atol=1e-10)


def test_RotationSequence3D_errors():
    # Bad axes_order labels
    with pytest.raises(ValueError, match=r"Unrecognized axis label .* should be one of .*"):
        rotations.RotationSequence3D(mk.MagicMock(), axes_order="abc")

    # Bad number of angles
    with pytest.raises(ValueError) as err:
        rotations.RotationSequence3D([1, 2, 3, 4], axes_order="zyx")
    assert str(err.value) == "The number of angles 4 should match the number of axes 3."

    # Bad evaluation input shapes
    model = rotations.RotationSequence3D([1, 2, 3], axes_order="zyx")
    message = "Expected input arrays to have the same shape"
    with pytest.raises(ValueError) as err:
        model.evaluate(np.array([1, 2, 3]),
                       np.array([1, 2]),
                       np.array([1, 2]),
                       [1, 2, 3])
    assert str(err.value) == message
    with pytest.raises(ValueError) as err:
        model.evaluate(np.array([1, 2]),
                       np.array([1, 2, 3]),
                       np.array([1, 2]),
                       [1, 2, 3])
    assert str(err.value) == message
    with pytest.raises(ValueError) as err:
        model.evaluate(np.array([1, 2]),
                       np.array([1, 2]),
                       np.array([1, 2, 3]),
                       [1, 2, 3])
    assert str(err.value) == message


def test_RotationSequence3D_inverse():
    model = rotations.RotationSequence3D([1, 2, 3], axes_order="zyx")

    assert_allclose(model.inverse.angles.value, [-3, -2, -1])
    assert model.inverse.axes_order == "xyz"


def test_EulerAngleRotation_errors():
    # Bad length of axes_order
    with pytest.raises(TypeError) as err:
        rotations.EulerAngleRotation(mk.MagicMock(), mk.MagicMock(), mk.MagicMock(),
                                     axes_order="xyzx")
    assert str(err.value) == "Expected axes_order to be a character sequence of length 3, got xyzx"

    # Bad axes_order labels
    with pytest.raises(ValueError, match=r"Unrecognized axis label .* should be one of .*"):
        rotations.EulerAngleRotation(mk.MagicMock(), mk.MagicMock(), mk.MagicMock(),
                                     axes_order="abc")

    # Bad units
    message = "All parameters should be of the same type - float or Quantity."
    with pytest.raises(TypeError) as err:
        rotations.EulerAngleRotation(1 * u.m, 2, 3,
                                     axes_order="xyz")
    assert str(err.value) == message
    with pytest.raises(TypeError) as err:
        rotations.EulerAngleRotation(1, 2 * u.m, 3,
                                     axes_order="xyz")
    assert str(err.value) == message
    with pytest.raises(TypeError) as err:
        rotations.EulerAngleRotation(1, 2, 3 * u.m,
                                     axes_order="xyz")
    assert str(err.value) == message


def test_EulerAngleRotation_inverse():
    model = rotations.EulerAngleRotation(1, 2, 3, "xyz")

    assert_allclose(model.inverse.phi, -3)
    assert_allclose(model.inverse.theta, -2)
    assert_allclose(model.inverse.psi, -1)
    assert model.inverse.axes_order == "zyx"


def test__SkyRotation_errors():
    # Bad units
    message = "All parameters should be of the same type - float or Quantity."
    with pytest.raises(TypeError) as err:
        rotations._SkyRotation(1 * u.m, 2, 3)
    assert str(err.value) == message
    with pytest.raises(TypeError) as err:
        rotations._SkyRotation(1, 2 * u.m, 3)
    assert str(err.value) == message
    with pytest.raises(TypeError) as err:
        rotations._SkyRotation(1, 2, 3 * u.m)
    assert str(err.value) == message


def test__SkyRotation__evaluate():
    model = rotations._SkyRotation(1, 2, 3)

    phi = mk.MagicMock()
    theta = mk.MagicMock()
    lon = mk.MagicMock()
    lat = mk.MagicMock()
    lon_pole = mk.MagicMock()

    alpha = 5
    delta = mk.MagicMock()
    with mk.patch.object(rotations._EulerRotation, 'evaluate',
                         autospec=True, return_value=(alpha, delta)) as mkEval:
        assert (365, delta) == model._evaluate(phi, theta, lon, lat, lon_pole)
        assert mkEval.call_args_list == [mk.call(model, phi, theta, lon, lat, lon_pole, 'zxz')]

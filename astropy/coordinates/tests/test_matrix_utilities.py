# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.coordinates.matrix_utilities import (
    angle_axis,
    is_O3,
    is_rotation,
    rotation_matrix,
)


def test_rotation_matrix():
    assert_allclose(rotation_matrix(0 * u.deg, "x"), np.eye(3))

    assert_allclose(
        rotation_matrix(90 * u.deg, "y"), [[0, 0, -1], [0, 1, 0], [1, 0, 0]], atol=1e-12
    )

    assert_allclose(
        rotation_matrix(-90 * u.deg, "z"),
        [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
        atol=1e-12,
    )

    assert_allclose(
        rotation_matrix(45 * u.deg, "x"), rotation_matrix(45 * u.deg, [1, 0, 0])
    )
    assert_allclose(
        rotation_matrix(125 * u.deg, "y"), rotation_matrix(125 * u.deg, [0, 1, 0])
    )
    assert_allclose(
        rotation_matrix(-30 * u.deg, "z"), rotation_matrix(-30 * u.deg, [0, 0, 1])
    )

    assert_allclose(
        np.dot(rotation_matrix(180 * u.deg, [1, 1, 0]), [1, 0, 0]),
        [0, 1, 0],
        atol=1e-12,
    )

    # make sure it also works for very small angles
    assert_allclose(
        rotation_matrix(0.000001 * u.deg, "x"),
        rotation_matrix(0.000001 * u.deg, [1, 0, 0]),
    )


def test_angle_axis():
    m1 = rotation_matrix(35 * u.deg, "x")
    an1, ax1 = angle_axis(m1)

    assert an1 - 35 * u.deg < 1e-10 * u.deg
    assert_allclose(ax1, [1, 0, 0])

    m2 = rotation_matrix(-89 * u.deg, [1, 1, 0])
    an2, ax2 = angle_axis(m2)

    assert an2 - 89 * u.deg < 1e-10 * u.deg
    assert_allclose(ax2, [-(2**-0.5), -(2**-0.5), 0])


def test_is_O3():
    """Test the matrix checker ``is_O3``."""
    # Normal rotation matrix
    m1 = rotation_matrix(35 * u.deg, "x")
    assert is_O3(m1)
    # and (M, 3, 3)
    n1 = np.tile(m1, (2, 1, 1))
    assert tuple(is_O3(n1)) == (True, True)  # (show the broadcasting)
    # Test atol parameter
    nn1 = np.tile(0.5 * m1, (2, 1, 1))
    assert tuple(is_O3(nn1)) == (False, False)  # (show the broadcasting)
    assert tuple(is_O3(nn1, atol=1)) == (True, True)  # (show the broadcasting)

    # reflection
    m2 = m1.copy()
    m2[0, 0] *= -1
    assert is_O3(m2)
    # and (M, 3, 3)
    n2 = np.stack((m1, m2))
    assert tuple(is_O3(n2)) == (True, True)  # (show the broadcasting)

    # Not any sort of O(3)
    m3 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert not is_O3(m3)
    # and (M, 3, 3)
    n3 = np.stack((m1, m3))
    assert tuple(is_O3(n3)) == (True, False)  # (show the broadcasting)


def test_is_rotation():
    """Test the rotation matrix checker ``is_rotation``."""
    # Normal rotation matrix
    m1 = rotation_matrix(35 * u.deg, "x")
    assert is_rotation(m1)
    assert is_rotation(m1, allow_improper=True)  # (a less restrictive test)
    # and (M, 3, 3)
    n1 = np.tile(m1, (2, 1, 1))
    assert tuple(is_rotation(n1)) == (True, True)  # (show the broadcasting)
    # Test atol parameter
    nn1 = np.tile(0.5 * m1, (2, 1, 1))
    assert tuple(is_rotation(nn1)) == (False, False)  # (show the broadcasting)
    assert tuple(is_rotation(nn1, atol=10)) == (True, True)  # (show the broadcasting)

    # Improper rotation (unit rotation + reflection)
    m2 = np.identity(3)
    m2[0, 0] = -1
    assert not is_rotation(m2)
    assert is_rotation(m2, allow_improper=True)
    # and (M, 3, 3)
    n2 = np.stack((m1, m2))
    assert tuple(is_rotation(n2)) == (True, False)  # (show the broadcasting)

    # Not any sort of rotation
    m3 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert not is_rotation(m3)
    assert not is_rotation(m3, allow_improper=True)
    # and (M, 3, 3)
    n3 = np.stack((m1, m3))
    assert tuple(is_rotation(n3)) == (True, False)  # (show the broadcasting)

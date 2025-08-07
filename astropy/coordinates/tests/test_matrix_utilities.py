# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from astropy import units as u
from astropy.coordinates import Angle, is_rotation_or_reflection, rotation_matrix
from astropy.coordinates.matrix_utilities import (
    angle_axis,
    is_O3,
    is_rotation,
)
from astropy.utils.exceptions import AstropyDeprecationWarning

ROTATION_MATRIX = rotation_matrix(35 * u.deg, "x")
IMPROPER_ROTATION_MATRIX = ROTATION_MATRIX.copy()
IMPROPER_ROTATION_MATRIX[0, 0] *= -1
NON_ORTHOGONAL_MATRIX = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])


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


@pytest.mark.parametrize(
    "angle", [Angle(0 * u.deg), 0 * u.deg, np.array(0), 0, "0 deg"], ids=type
)
def test_rotation_angle_input_types(angle):
    assert_array_equal(rotation_matrix(angle), np.eye(3), strict=True)


def test_angle_axis():
    m1 = rotation_matrix(35 * u.deg, "x")
    with pytest.warns(AstropyDeprecationWarning):
        an1, ax1 = angle_axis(m1)

    assert an1 - 35 * u.deg < 1e-10 * u.deg
    assert_allclose(ax1, [1, 0, 0])

    m2 = rotation_matrix(-89 * u.deg, [1, 1, 0])
    with pytest.warns(AstropyDeprecationWarning):
        an2, ax2 = angle_axis(m2)

    assert an2 - 89 * u.deg < 1e-10 * u.deg
    assert_allclose(ax2, [-(2**-0.5), -(2**-0.5), 0])


@pytest.mark.parametrize(
    "matrix,expectation",
    [
        pytest.param(ROTATION_MATRIX, np.True_, id="rotation"),
        pytest.param(IMPROPER_ROTATION_MATRIX, np.True_, id="improper rotation"),
        pytest.param(
            np.stack((ROTATION_MATRIX, IMPROPER_ROTATION_MATRIX)),
            (np.True_, np.True_),
            id="proper and improper rotation",
        ),
        pytest.param(NON_ORTHOGONAL_MATRIX, np.False_, id="non-orthogonal matrix"),
        pytest.param(
            np.stack((ROTATION_MATRIX, NON_ORTHOGONAL_MATRIX)),
            (np.True_, np.False_),
            id="rotation and non-orthogonal matrix",
        ),
    ],
)
def test_is_rotation_or_reflection(matrix, expectation):
    assert_array_equal(is_rotation_or_reflection(matrix), expectation, strict=True)


@pytest.mark.parametrize(
    "atol,expectation",
    [
        pytest.param(None, (np.True_, np.False_), id="default atol"),
        pytest.param(1, (np.True_, np.True_), id="very large atol"),
    ],
)
def test_is_rotation_or_reflection_atol(atol, expectation):
    assert_array_equal(
        is_rotation_or_reflection(
            np.stack([ROTATION_MATRIX, 0.5 * ROTATION_MATRIX]), atol
        ),
        expectation,
        strict=True,
    )


def test_is_O3():
    with pytest.warns(
        AstropyDeprecationWarning, match=r"Use is_rotation_or_reflection instead\.$"
    ):
        is_O3(ROTATION_MATRIX)


def test_is_rotation():
    """Test the rotation matrix checker ``is_rotation``."""
    # Normal rotation matrix
    m1 = rotation_matrix(35 * u.deg, "x")
    with pytest.warns(AstropyDeprecationWarning):
        assert is_rotation(m1)
    with pytest.warns(AstropyDeprecationWarning):
        assert is_rotation(m1, allow_improper=True)  # (a less restrictive test)
    # and (M, 3, 3)
    n1 = np.tile(m1, (2, 1, 1))
    with pytest.warns(AstropyDeprecationWarning):
        assert tuple(is_rotation(n1)) == (True, True)  # (show the broadcasting)
    # Test atol parameter
    nn1 = np.tile(0.5 * m1, (2, 1, 1))
    with pytest.warns(AstropyDeprecationWarning):
        assert tuple(is_rotation(nn1)) == (False, False)  # (show the broadcasting)
    with pytest.warns(AstropyDeprecationWarning):
        assert tuple(is_rotation(nn1, atol=10)) == (True, True)

    # Improper rotation (unit rotation + reflection)
    m2 = np.identity(3)
    m2[0, 0] = -1
    with pytest.warns(AstropyDeprecationWarning):
        assert not is_rotation(m2)
    with pytest.warns(AstropyDeprecationWarning):
        assert is_rotation(m2, allow_improper=True)
    # and (M, 3, 3)
    n2 = np.stack((m1, m2))
    with pytest.warns(AstropyDeprecationWarning):
        assert tuple(is_rotation(n2)) == (True, False)  # (show the broadcasting)

    # Not any sort of rotation
    m3 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    with pytest.warns(AstropyDeprecationWarning):
        assert not is_rotation(m3)
    with pytest.warns(AstropyDeprecationWarning):
        assert not is_rotation(m3, allow_improper=True)
    # and (M, 3, 3)
    n3 = np.stack((m1, m3))
    with pytest.warns(AstropyDeprecationWarning):
        assert tuple(is_rotation(n3)) == (True, False)  # (show the broadcasting)

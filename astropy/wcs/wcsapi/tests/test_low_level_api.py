import numpy as np
import pytest
from numpy.testing import assert_equal

from astropy.wcs.wcsapi.low_level_api import BaseLowLevelWCS, validate_physical_types


def test_validate_physical_types():
    # Check valid cases
    validate_physical_types(["pos.eq.ra", "pos.eq.ra"])
    validate_physical_types(["spect.dopplerVeloc.radio", "custom:spam"])
    validate_physical_types(["time", None])

    # Make sure validation is case sensitive
    with pytest.raises(
        ValueError, match=r"'Pos\.eq\.dec' is not a valid IOVA UCD1\+ physical type"
    ):
        validate_physical_types(["pos.eq.ra", "Pos.eq.dec"])

    # Make sure nonsense types are picked up
    with pytest.raises(
        ValueError, match=r"'spam' is not a valid IOVA UCD1\+ physical type"
    ):
        validate_physical_types(["spam"])


class MatrixLowLevelWCS(BaseLowLevelWCS):
    """
    Minimal low-level WCS whose only meaningful property is the axis
    correlation matrix, used to test the default reverse matrix.
    """

    def __init__(self, matrix):
        self._matrix = np.asarray(matrix, dtype=bool)

    @property
    def pixel_n_dim(self):
        return self._matrix.shape[1]

    @property
    def world_n_dim(self):
        return self._matrix.shape[0]

    @property
    def axis_correlation_matrix(self):
        return self._matrix

    @property
    def world_axis_physical_types(self):
        return [None] * self.world_n_dim

    @property
    def world_axis_units(self):
        return [""] * self.world_n_dim

    # Abstract members that are not needed for these tests
    pixel_to_world_values = world_to_pixel_values = None
    world_axis_object_components = world_axis_object_classes = None


# Each tuple gives the forward axis correlation matrix and the expected default
# reverse matrix (as 0/1, converted to bool in the test).
@pytest.mark.parametrize(
    ("forward", "expected"),
    [
        # Independent axes: the transpose
        ([[1, 0], [0, 1]], [[1, 0], [0, 1]]),
        # Fully coupled celestial pair
        ([[1, 1], [1, 1]], [[1, 1], [1, 1]]),
        # Spectral cube: two independent blocks
        ([[1, 1, 0], [1, 1, 0], [0, 0, 1]], [[1, 1, 0], [1, 1, 0], [0, 0, 1]]),
        # Same cube with the world axes in a different order
        ([[0, 0, 1], [1, 1, 0], [1, 1, 0]], [[0, 1, 1], [0, 1, 1], [1, 0, 0]]),
        # Blocks scattered across non-contiguous rows and columns
        ([[0, 1, 1], [1, 0, 0], [0, 1, 1]], [[0, 1, 0], [1, 0, 1], [1, 0, 1]]),
        # Triangular: w0 = f(p0), w1 = g(p0, p1), so p1 needs w0 as well as w1
        ([[1, 0], [1, 1]], [[1, 1], [1, 1]]),
        # Two pixel, three world axes with a triangular structure
        ([[1, 0], [1, 1], [1, 1]], [[1, 1, 1], [1, 1, 1]]),
        # Two pixel, three world axes, all coupled
        ([[1, 1], [1, 1], [1, 1]], [[1, 1, 1], [1, 1, 1]]),
        # Triangular block interleaved with an independent axis
        ([[1, 1, 0], [0, 0, 1], [0, 1, 0]], [[1, 0, 1], [1, 0, 1], [0, 1, 0]]),
        # World axis that depends on no pixel axis is never required
        ([[1, 0], [0, 1], [0, 0]], [[1, 0, 0], [0, 1, 0]]),
        # Pixel axis that no world axis depends on requires nothing
        ([[1, 0, 0], [0, 1, 0]], [[1, 0], [0, 1], [0, 0]]),
    ],
)
def test_default_reverse_axis_correlation_matrix(forward, expected):
    wcs = MatrixLowLevelWCS(forward)
    reverse = wcs.reverse_axis_correlation_matrix
    assert reverse.dtype == bool
    assert reverse.shape == (wcs.pixel_n_dim, wcs.world_n_dim)
    assert_equal(reverse, np.array(expected, dtype=bool))


def test_default_reverse_axis_correlation_matrix_all_true():
    # With no information about the forward matrix, everything is required
    class AllTrueWCS(MatrixLowLevelWCS):
        axis_correlation_matrix = BaseLowLevelWCS.axis_correlation_matrix

    wcs = AllTrueWCS(np.ones((3, 2)))
    assert_equal(wcs.reverse_axis_correlation_matrix, np.ones((2, 3), dtype=bool))

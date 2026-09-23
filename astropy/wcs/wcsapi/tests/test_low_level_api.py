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
    correlation matrix, used to test the default inverse matrix.
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

    def pixel_to_world_values(self, *pixel_arrays):
        raise NotImplementedError

    def world_to_pixel_values(self, *world_arrays):
        raise NotImplementedError

    @property
    def world_axis_object_components(self):
        raise NotImplementedError

    @property
    def world_axis_object_classes(self):
        raise NotImplementedError


T, F = True, False


@pytest.mark.parametrize(
    ("forward", "expected"),
    [
        # Independent axes: the transpose
        ([[T, F], [F, T]], [[T, F], [F, T]]),
        # Fully coupled celestial pair
        ([[T, T], [T, T]], [[T, T], [T, T]]),
        # Spectral cube: two independent blocks
        ([[T, T, F], [T, T, F], [F, F, T]], [[T, T, F], [T, T, F], [F, F, T]]),
        # Same cube with the world axes in a different order
        ([[F, F, T], [T, T, F], [T, T, F]], [[F, T, T], [F, T, T], [T, F, F]]),
        # Blocks scattered across non-contiguous rows and columns
        ([[F, T, T], [T, F, F], [F, T, T]], [[F, T, F], [T, F, T], [T, F, T]]),
        # Triangular: w0 = f(p0), w1 = g(p0, p1), so p1 needs w0 as well as w1
        ([[T, F], [T, T]], [[T, T], [T, T]]),
        # Two pixel, three world axes with a triangular structure
        ([[T, F], [T, T], [T, T]], [[T, T, T], [T, T, T]]),
        # Two pixel, three world axes, all coupled
        ([[T, T], [T, T], [T, T]], [[T, T, T], [T, T, T]]),
        # Triangular block interleaved with an independent axis
        ([[T, T, F], [F, F, T], [F, T, F]], [[T, F, T], [T, F, T], [F, T, F]]),
        # World axis that depends on no pixel axis is never required
        ([[T, F], [F, T], [F, F]], [[T, F, F], [F, T, F]]),
        # Pixel axis that no world axis depends on requires nothing
        ([[T, F, F], [F, T, F]], [[T, F], [F, T], [F, F]]),
    ],
)
def test_default_inverse_axis_correlation_matrix(forward, expected):
    wcs = MatrixLowLevelWCS(forward)
    inverse = wcs.inverse_axis_correlation_matrix
    assert inverse.dtype == bool
    assert inverse.shape == (wcs.pixel_n_dim, wcs.world_n_dim)
    assert_equal(inverse, expected)


def test_default_inverse_axis_correlation_matrix_all_true():
    # With no information about the forward matrix, everything is required
    class AllTrueWCS(MatrixLowLevelWCS):
        axis_correlation_matrix = BaseLowLevelWCS.axis_correlation_matrix

    wcs = AllTrueWCS(np.ones((3, 2)))
    assert_equal(wcs.inverse_axis_correlation_matrix, np.ones((2, 3), dtype=bool))

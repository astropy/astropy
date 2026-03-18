"""
Direct tests for the _convolveNd_c Cython extension in astropy/convolution/_convolve.pyx.

These tests call _convolveNd_c directly, bypassing the public convolve()
API, to verify the compiled extension's behavior at the interface level.
"""

import numpy as np
import pytest

from astropy.convolution._convolve import _convolveNd_c


def test_1d_simple_convolution():
    """Basic 1D convolution with a uniform kernel."""
    array = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=float, order="C")
    kernel = np.array([1.0, 1.0, 1.0], dtype=float, order="C") / 3.0
    result = np.zeros(array.shape, dtype=float, order="C")

    _convolveNd_c(result, array, kernel, False, True, 1)

    # Interior points only - boundary=None leaves edges as zero
    np.testing.assert_allclose(result[1], 2.0, rtol=1e-6)
    np.testing.assert_allclose(result[2], 3.0, rtol=1e-6)
    np.testing.assert_allclose(result[3], 4.0, rtol=1e-6)


def test_1d_nan_interpolate():
    """NaN interpolation: NaN replaced by kernel-weighted neighbors."""
    array = np.array([1.0, np.nan, 3.0, 4.0, 5.0], dtype=float, order="C")
    kernel = np.array([1.0, 1.0, 1.0], dtype=float, order="C") / 3.0
    result = np.zeros(array.shape, dtype=float, order="C")

    _convolveNd_c(result, array, kernel, True, True, 1)

    # Index 1: NaN interpolated from neighbors 1.0 and 3.0 -> 2.0
    np.testing.assert_allclose(result[1], 2.0, rtol=1e-6)
    # Index 2: neighbors are nan and 3.0, 4.0 -> 3.5
    np.testing.assert_allclose(result[2], 3.5, rtol=1e-6)


def test_1d_no_nan_interpolate():
    """With nan_interpolate=False, NaN propagates into result."""
    array = np.array([1.0, np.nan, 3.0, 4.0, 5.0], dtype=float, order="C")
    kernel = np.array([1.0, 1.0, 1.0], dtype=float, order="C") / 3.0
    result = np.zeros(array.shape, dtype=float, order="C")

    _convolveNd_c(result, array, kernel, False, True, 1)

    assert np.isnan(result[1])
    assert np.isnan(result[2])


def test_2d_simple_convolution():
    """Basic 2D convolution with identity kernel."""
    array = np.array(
        [[1.0, 2.0, 3.0],
         [4.0, 5.0, 6.0],
         [7.0, 8.0, 9.0]], dtype=float, order="C"
    )
    kernel = np.array(
        [[0.0, 0.0, 0.0],
         [0.0, 1.0, 0.0],
         [0.0, 0.0, 0.0]], dtype=float, order="C"
    )
    result = np.zeros(array.shape, dtype=float, order="C")

    _convolveNd_c(result, array, kernel, False, True, 1)

    # Identity kernel - center element only
    np.testing.assert_allclose(result[1, 1], 5.0, rtol=1e-6)


def test_result_shape_matches_input():
    """Output shape must match input array shape."""
    array = np.ones((7,), dtype=float, order="C")
    kernel = np.array([1.0, 1.0, 1.0], dtype=float, order="C") / 3.0
    result = np.zeros(array.shape, dtype=float, order="C")

    _convolveNd_c(result, array, kernel, False, True, 1)

    assert result.shape == array.shape


def test_embed_result_within_padded_region_false():
    """embed_result_within_padded_region=False with padded input."""
    # Simulate boundary='fill' - array is pre-padded, result is full size
    array = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0], dtype=float, order="C")
    kernel = np.array([1.0, 1.0, 1.0], dtype=float, order="C") / 3.0
    result = np.zeros((5,), dtype=float, order="C")

    _convolveNd_c(result, array, kernel, False, False, 1)

    np.testing.assert_allclose(result[1], 2.0, rtol=1e-6)

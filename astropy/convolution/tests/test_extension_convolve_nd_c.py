# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np
import numpy.testing as npt
import pytest

from astropy.convolution._convolve import _convolveNd_c


# these helprer functions preprocess the inputs before passing to the extension
def _as_float64_c_contiguous(array):
    """Return a C-contiguous ``float64`` array for the extension call."""
    return np.ascontiguousarray(array, dtype=np.float64)


def _pad_for_direct_call(image, kernel):
    """Pad an image the same way with half-kernel width."""
    pad_width = tuple((size // 2, size // 2) for size in kernel.shape)
    padded = np.pad(image, pad_width, mode="constant", constant_values=0.0)
    return _as_float64_c_contiguous(padded)


def _allocate_result(image, padded_image, embed):
    """Allocate the result buffer with the shape expected by the C code."""

    shape = padded_image.shape if embed else image.shape
    return np.zeros(shape, dtype=np.float64, order="C")


def _call_direct(image, kernel, *, nan_interpolate=False, embed=False, n_threads=1):
    """Prepare inputs and call ``_convolveNd_c`` directly."""

    image = _as_float64_c_contiguous(image)
    kernel = _as_float64_c_contiguous(kernel)
    padded_image = _pad_for_direct_call(image, kernel)
    result = _allocate_result(image, padded_image, embed)

    _convolveNd_c(
        result,
        padded_image,
        kernel,
        nan_interpolate,
        embed,
        n_threads,
    )

    return result


def _trim_padding(embedded_result, kernel):
    """Remove the outer padding region from an embedded result array after convolution."""

    slices = []
    for size in kernel.shape:
        width = size // 2
        slices.append(slice(width, -width if width else None))
    return embedded_result[tuple(slices)]


# this function acts as the oracle, a reference against which the outputs of the _convolveNd_c will be validated
def _reference_direct_convolution(image, kernel, *, nan_interpolate=False, embed=False):
    """Pure NumPy reference for the direct extension contract."""
    image = _as_float64_c_contiguous(image)
    kernel = _as_float64_c_contiguous(kernel)
    padded_image = _pad_for_direct_call(image, kernel)
    flipped_kernel = kernel[tuple(slice(None, None, -1) for _ in range(kernel.ndim))]
    pad_width = tuple(size // 2 for size in kernel.shape)
    result = _allocate_result(image, padded_image, embed)

    for image_index in np.ndindex(image.shape):
        padded_index = tuple(index + pad for index, pad in zip(image_index, pad_width))

        window_slices = tuple(
            slice(
                padded_index[axis] - pad_width[axis],
                padded_index[axis] - pad_width[axis] + kernel.shape[axis],
            )
            for axis in range(image.ndim)
        )
        window = padded_image[window_slices]

        if nan_interpolate:
            valid = ~np.isnan(window)
            if np.any(valid):
                top = np.sum(window[valid] * flipped_kernel[valid])
                bot = np.sum(flipped_kernel[valid])
            else:
                top = 0.0
                bot = 0.0

            # The C code copies the center input value when the effective
            # kernel weight sum is exactly zero.
            if bot == 0.0:
                value = padded_image[padded_index]
            else:
                value = top / bot
        else:
            value = np.sum(window * flipped_kernel)

        result_index = padded_index if embed else image_index
        result[result_index] = value

    return result


@pytest.mark.parametrize(
    ("image", "kernel"),
    [
        pytest.param(
            np.array([10.0, 20.0, 30.0, 40.0, 50.0]),
            np.array([1.0, 2.0, 1.0]),
            id="1d",
        ),
        pytest.param(
            np.arange(1.0, 17.0).reshape(4, 4),
            np.array([[0.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 0.0]]),
            id="2d",
        ),
        pytest.param(
            np.arange(1.0, 65.0).reshape(4, 4, 4),
            np.ones((3, 3, 3), dtype=np.float64),
            id="3d",
        ),
    ],
)
def test_dispatch_paths_match_reference_without_nan_interpolation(image, kernel):
    """Each dispatch path(1D, 2D, 3D) should agree with the NumPy reference result."""
    result = _call_direct(image, kernel, nan_interpolate=False, embed=False)
    expected = _reference_direct_convolution(
        image,
        kernel,
        nan_interpolate=False,
        embed=False,
    )

    npt.assert_allclose(result, expected)


def test_asymmetric_kernel_is_flipped_before_multiplication():
    """The direct extension should perform convolution, not correlation."""
    image = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    kernel = np.array([1.0, 2.0, 3.0])

    result = _call_direct(image, kernel, nan_interpolate=False, embed=False)

    # ``convolve.c`` flips the kernel, so [1, 2, 3] is applied as [3, 2, 1].
    expected = np.array([4.0, 10.0, 16.0, 22.0, 22.0])
    npt.assert_allclose(result, expected)


@pytest.mark.parametrize(
    ("image", "kernel"),
    [
        pytest.param(
            np.array([3.0, 7.0, 2.0, 9.0]), np.array([1.0, 2.0, 1.0]), id="1d"
        ),
        pytest.param(
            np.arange(1.0, 10.0).reshape(3, 3),
            np.ones((3, 3), dtype=np.float64),
            id="2d",
        ),
        pytest.param(
            np.arange(1.0, 28.0).reshape(3, 3, 3),
            np.ones((3, 3, 3), dtype=np.float64),
            id="3d",
        ),
    ],
)
def test_embed_true_matches_embed_false_in_inner_region_and_leaves_border_zero(
    image, kernel
):
    """Embedded results should land in the padded coordinates and nowhere else."""
    embedded = _call_direct(image, kernel, nan_interpolate=False, embed=True)
    non_embedded = _call_direct(image, kernel, nan_interpolate=False, embed=False)

    # The inner region should be exactly the same data we get from embed=False.
    npt.assert_allclose(_trim_padding(embedded, kernel), non_embedded)

    # The C loop only writes valid image coordinates, so the outer border should
    # still contain the zeros from the initial result allocation.
    inner_slices = tuple(
        slice(size // 2, -(size // 2) if size // 2 else None) for size in kernel.shape
    )
    padding_mask = np.ones(embedded.shape, dtype=bool)
    padding_mask[inner_slices] = False
    npt.assert_array_equal(embedded[padding_mask], 0.0)


@pytest.mark.parametrize(
    ("image", "kernel", "probe_index"),
    [
        pytest.param(np.array([1.0, np.nan, 3.0, 4.0, 5.0]), np.ones(3), (0,), id="1d"),
        pytest.param(
            np.array(
                [
                    [1.0, 2.0, 3.0, 4.0, 5.0],
                    [6.0, 7.0, 8.0, 9.0, 10.0],
                    [11.0, 12.0, np.nan, 14.0, 15.0],
                    [16.0, 17.0, 18.0, 19.0, 20.0],
                    [21.0, 22.0, 23.0, 24.0, 25.0],
                ]
            ),
            np.ones((3, 3)),
            (2, 2),
            id="2d",
        ),
        pytest.param(
            np.ones((5, 5, 5), dtype=np.float64),
            np.ones((3, 3, 3), dtype=np.float64),
            (2, 2, 2),
            id="3d",
        ),
    ],
)
def test_nan_propagates_when_nan_interpolation_is_disabled(image, kernel, probe_index):
    """Without NaN interpolation, any NaN in the window should poison the sum."""
    if image.ndim == 3:
        image = image.copy()
        image[2, 2, 2] = np.nan

    result = _call_direct(image, kernel, nan_interpolate=False, embed=False)
    expected = _reference_direct_convolution(
        image,
        kernel,
        nan_interpolate=False,
        embed=False,
    )

    assert np.isnan(result[probe_index])
    npt.assert_allclose(result, expected, equal_nan=True)


@pytest.mark.parametrize(
    ("image", "kernel", "probe_index"),
    [
        pytest.param(np.array([1.0, 1.0, np.nan, 1.0, 1.0]), np.ones(3), (2,), id="1d"),
        pytest.param(
            np.ones((5, 5), dtype=np.float64),
            np.ones((3, 3), dtype=np.float64),
            (2, 2),
            id="2d",
        ),
        pytest.param(
            np.ones((5, 5, 5), dtype=np.float64),
            np.ones((3, 3, 3), dtype=np.float64),
            (2, 2, 2),
            id="3d",
        ),
    ],
)
def test_nan_interpolation_skips_nan_values_and_renormalizes(
    image, kernel, probe_index
):
    """With NaN interpolation enabled, valid neighbors should be reweighted."""
    image = image.copy()
    image[probe_index] = np.nan

    result = _call_direct(image, kernel, nan_interpolate=True, embed=False)
    expected = _reference_direct_convolution(
        image,
        kernel,
        nan_interpolate=True,
        embed=False,
    )

    # A window of ones with one missing value should still evaluate to 1.0.
    assert result[probe_index] == pytest.approx(1.0)
    npt.assert_allclose(result, expected, equal_nan=True)


def test_all_nan_window_copies_the_center_nan_value():
    """An all-NaN window should fall back to the center input value."""
    image = np.array([1.0, np.nan, np.nan, np.nan, 5.0])
    kernel = np.ones(3, dtype=np.float64)

    result = _call_direct(image, kernel, nan_interpolate=True, embed=False)
    expected = _reference_direct_convolution(
        image,
        kernel,
        nan_interpolate=True,
        embed=False,
    )

    # The center window is [nan, nan, nan], so ``bot`` becomes zero and the
    # C code copies the center value, which is also NaN in this case.
    assert np.isnan(result[2])
    npt.assert_allclose(result, expected, equal_nan=True)


def test_zero_effective_weight_sum_copies_the_center_input_value():
    """The ``bot == 0`` fallback should also work when valid pixels remain."""
    image = np.array([5.0, np.nan, 20.0, 30.0, 40.0])
    kernel = np.array([1.0, -1.0, 0.0])

    result = _call_direct(image, kernel, nan_interpolate=True, embed=False)
    expected = _reference_direct_convolution(
        image,
        kernel,
        nan_interpolate=True,
        embed=False,
    )

    # At index 2 the valid weights sum to zero after the NaN is skipped, so the
    # extension should copy the center input value instead of dividing by zero.
    assert result[2] == pytest.approx(20.0)
    npt.assert_allclose(result, expected, equal_nan=True)


@pytest.mark.parametrize(
    ("image", "kernel"),
    [
        pytest.param(np.array([1.0, 2.0, 3.0]), np.array([1.0, 0.0, -1.0]), id="1d"),
        pytest.param(
            np.arange(1.0, 10.0).reshape(3, 3),
            np.ones((3, 3), dtype=np.float64),
            id="2d",
        ),
        pytest.param(
            np.arange(1.0, 28.0).reshape(3, 3, 3),
            np.ones((3, 3, 3), dtype=np.float64),
            id="3d",
        ),
    ],
)
def test_smallest_valid_input_shape_still_computes(image, kernel):
    """Images that are just large enough for the kernel should still work."""
    result = _call_direct(image, kernel, nan_interpolate=False, embed=False)
    expected = _reference_direct_convolution(
        image,
        kernel,
        nan_interpolate=False,
        embed=False,
    )

    npt.assert_allclose(result, expected)


@pytest.mark.parametrize(
    ("image", "kernel"),
    [
        pytest.param(np.linspace(1.0, 9.0, 9), np.array([1.0, 2.0, 1.0]), id="1d"),
        pytest.param(
            np.arange(1.0, 101.0).reshape(10, 10),
            np.ones((3, 3), dtype=np.float64),
            id="2d",
        ),
        pytest.param(
            np.arange(1.0, 513.0).reshape(8, 8, 8),
            np.ones((3, 3, 3), dtype=np.float64),
            id="3d",
        ),
    ],
)
def test_multi_threaded_results_match_single_threaded(image, kernel):
    """Changing ``n_threads`` should not change the numerical result."""
    single_thread = _call_direct(
        image,
        kernel,
        nan_interpolate=False,
        embed=False,
        n_threads=1,
    )
    multi_thread = _call_direct(
        image,
        kernel,
        nan_interpolate=False,
        embed=False,
        n_threads=4,
    )

    npt.assert_allclose(single_thread, multi_thread)

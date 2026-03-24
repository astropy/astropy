# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np
import numpy.testing as npt
import pytest

from astropy.convolution._convolve import _convolveNd_c


def _run_extension_with_zero_padding(
    image, kernel, *, nan_interpolate=False, embed=False, n_threads=1
):
    """Call ``_convolveNd_c`` with a zero-padded input."""
    image = np.ascontiguousarray(image, dtype=np.float64)
    kernel = np.ascontiguousarray(kernel, dtype=np.float64)
    pad_width = tuple((size // 2, size // 2) for size in kernel.shape)
    padded_image = np.pad(
        image,
        pad_width=pad_width,
    )
    result_shape = padded_image.shape if embed else image.shape
    result = np.zeros(result_shape, dtype=np.float64, order="C")

    _convolveNd_c(
        result,
        padded_image,
        kernel,
        nan_interpolate,
        embed,
        n_threads,
    )

    return result


@pytest.mark.parametrize(
    ("image", "kernel", "expected"),
    [
        pytest.param(
            np.array([1.0, 2.0, 3.0], dtype=np.float64),
            np.array([1.0, 0.0, -1.0], dtype=np.float64),
            np.array([2.0, 2.0, -2.0], dtype=np.float64),
            id="1d",
        ),
        pytest.param(
            np.arange(1.0, 10.0, dtype=np.float64).reshape(3, 3),
            np.ones((3, 3), dtype=np.float64),
            np.array(
                [
                    [12.0, 21.0, 16.0],
                    [27.0, 45.0, 33.0],
                    [24.0, 39.0, 28.0],
                ],
                dtype=np.float64,
            ),
            id="2d",
        ),
        pytest.param(
            np.arange(1.0, 28.0, dtype=np.float64).reshape(3, 3, 3),
            np.ones((3, 3, 3), dtype=np.float64),
            np.array(
                [
                    [[60.0, 96.0, 68.0], [108.0, 171.0, 120.0], [84.0, 132.0, 92.0]],
                    [
                        [144.0, 225.0, 156.0],
                        [243.0, 378.0, 261.0],
                        [180.0, 279.0, 192.0],
                    ],
                    [
                        [132.0, 204.0, 140.0],
                        [216.0, 333.0, 228.0],
                        [156.0, 240.0, 164.0],
                    ],
                ],
                dtype=np.float64,
            ),
            id="3d",
        ),
    ],
)
def test_zero_padded_convolution_matches_expected_values(image, kernel, expected):
    """Small representative inputs should produce the expected output values."""
    result = _run_extension_with_zero_padding(image, kernel)

    npt.assert_allclose(result, expected)


@pytest.mark.parametrize(
    ("image", "kernel", "expected"),
    [
        pytest.param(
            np.arange(1.0, 6.0, dtype=np.float64),
            np.arange(1.0, 4.0, dtype=np.float64),
            np.array([4.0, 10.0, 16.0, 22.0, 22.0], dtype=np.float64),
            id="1d",
        ),
        pytest.param(
            np.arange(1.0, 26.0, dtype=np.float64).reshape(5, 5),
            np.arange(1.0, 10.0, dtype=np.float64).reshape(3, 3),
            np.array(
                [
                    [32.0, 68.0, 89.0, 110.0, 96.0],
                    [114.0, 219.0, 264.0, 309.0, 252.0],
                    [249.0, 444.0, 489.0, 534.0, 417.0],
                    [384.0, 669.0, 714.0, 759.0, 582.0],
                    [440.0, 734.0, 773.0, 812.0, 600.0],
                ],
                dtype=np.float64,
            ),
            id="2d",
        ),
        pytest.param(
            np.arange(1.0, 28.0, dtype=np.float64).reshape(3, 3, 3),
            np.arange(1.0, 28.0, dtype=np.float64).reshape(3, 3, 3),
            np.array(
                [
                    [
                        [268.0, 490.0, 396.0],
                        [654.0, 1140.0, 882.0],
                        [700.0, 1174.0, 876.0],
                    ],
                    [
                        [1050.0, 1788.0, 1350.0],
                        [2196.0, 3654.0, 2700.0],
                        [2022.0, 3300.0, 2394.0],
                    ],
                    [
                        [1996.0, 3190.0, 2268.0],
                        [3570.0, 5676.0, 4014.0],
                        [2860.0, 4522.0, 3180.0],
                    ],
                ],
                dtype=np.float64,
            ),
            id="3d",
        ),
    ],
)
def test_asymmetric_kernel_is_flipped_correctly(image, kernel, expected):
    """
    Kernel flip is fundamental to convolution. The extension should flip
    the kernel correctly to compute the correct output. For a symmetric
    kernel like ``[1, 1, 1]`` the flip does not matter since it reads the
    same in both directions, but for an asymmetric kernel like ``[1, 2, 3]``
    the flip changes the output .
    """
    result = _run_extension_with_zero_padding(image, kernel)
    npt.assert_allclose(result, expected)


@pytest.mark.parametrize(
    ("image", "kernel", "expected"),
    [
        pytest.param(
            np.array([3.0, 7.0, 2.0, 9.0], dtype=np.float64),
            np.array([1.0, 2.0, 1.0], dtype=np.float64),
            np.array([0.0, 13.0, 19.0, 20.0, 20.0, 0.0], dtype=np.float64),
            id="1d",
        ),
        pytest.param(
            np.arange(1.0, 10.0, dtype=np.float64).reshape(3, 3),
            np.ones((3, 3), dtype=np.float64),
            np.array(
                [
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 12.0, 21.0, 16.0, 0.0],
                    [0.0, 27.0, 45.0, 33.0, 0.0],
                    [0.0, 24.0, 39.0, 28.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                ],
                dtype=np.float64,
            ),
            id="2d",
        ),
    ],
)
def test_embed_true_returns_expected_padded_output(image, kernel, expected):
    """With ``embed=True``, the result is written into the interior of the padded output."""
    result = _run_extension_with_zero_padding(image, kernel, embed=True)

    npt.assert_allclose(result, expected)


# fmt: off
@pytest.mark.parametrize(
    ("image", "kernel", "expected"),
    [
        pytest.param(
            np.array([1.0, np.nan, 3.0, 4.0, 5.0], dtype=np.float64),
            np.ones(3, dtype=np.float64),
            np.array([np.nan, np.nan, np.nan, 12.0, 9.0], dtype=np.float64),
            id="1d",
        ),
        pytest.param(
            np.array(
                [
                    [1.0, 2.0, 3.0, 4.0],
                    [5.0, np.nan, 7.0, 8.0],
                    [9.0, 10.0, 11.0, 12.0],
                    [13.0, 14.0, 15.0, 16.0],
                ],
                dtype=np.float64,
            ),
            np.ones((3, 3), dtype=np.float64),
            np.array(
                [
                    [np.nan, np.nan, np.nan, 22.0],
                    [np.nan, np.nan, np.nan, 45.0],
                    [np.nan, np.nan, np.nan, 69.0],
                    [46.0, 72.0, 78.0, 54.0],
                ],
                dtype=np.float64,
            ),
            id="2d",
        ),
        pytest.param(
            np.array(
                [
                    [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
                    [[1.0, 1.0, 1.0], [1.0, np.nan, 1.0], [1.0, 1.0, 1.0]],
                    [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
                ],
                dtype=np.float64,
            ),
           np.ones((3, 3, 3), dtype=np.float64)
,
           np.full((3, 3, 3), np.nan, dtype=np.float64),
            id="3d",
        ),
    ],
)
def test_nan_without_interpolation_matches_expected_values(image, kernel, expected):
    """Without NaN interpolation, any NaN in the window propagates.

    Example:

    for padded ``image = [0, 1, nan, 1, 0]`` and ``kernel = [1, 1, 1]``,
    ``_convolveNd_c`` yields ``result = [nan, nan, nan]`` because each
    output window contains the ``nan`` term and the extension
    intentionally does not skip NaN values unless
    ``nan_interpolate=True``.
    """

    result = _run_extension_with_zero_padding(image, kernel)

    npt.assert_allclose(result, expected, equal_nan=True)
# fmt: on


@pytest.mark.parametrize(
    ("image", "kernel", "expected"),
    [
        pytest.param(
            np.array([1.0, 1.0, np.nan, 1.0, 1.0], dtype=np.float64),
            np.ones(3, dtype=np.float64),
            np.array([2.0 / 3.0, 1.0, 1.0, 1.0, 2.0 / 3.0], dtype=np.float64),
            id="1d",
        ),
        pytest.param(
            np.array(
                [
                    [1.0, 1.0, 1.0],
                    [1.0, np.nan, 1.0],
                    [1.0, 1.0, 1.0],
                ],
                dtype=np.float64,
            ),
            np.ones((3, 3), dtype=np.float64),
            np.array(
                [
                    [0.375, 0.625, 0.375],
                    [0.625, 1.0, 0.625],
                    [0.375, 0.625, 0.375],
                ],
                dtype=np.float64,
            ),
            id="2d",
        ),
    ],
)
def test_nan_interpolation_matches_expected_values(image, kernel, expected):
    """With NaN interpolation enabled, valid neighbors are reweighted."""
    result = _run_extension_with_zero_padding(
        image,
        kernel,
        nan_interpolate=True,
    )

    npt.assert_allclose(result, expected)


def test_all_nan_window_copies_the_center_nan_value():
    """An all-NaN window should copy the center value back into the result when nan_interpolate is enabled."""
    image = np.array([1.0, np.nan, np.nan, np.nan, 5.0], dtype=np.float64)
    kernel = np.ones(3, dtype=np.float64)

    result = _run_extension_with_zero_padding(
        image,
        kernel,
        nan_interpolate=True,
    )
    expected = np.array([0.5, 1.0, np.nan, 5.0, 2.5], dtype=np.float64)

    npt.assert_allclose(result, expected, equal_nan=True)


def test_zero_effective_weight_sum_copies_the_center_input_value():
    """A zero effective kernel sum should copy the center input value when nan_interpolation is enabled."""

    image = np.array([5.0, np.nan, 20.0, 30.0, 40.0], dtype=np.float64)
    kernel = np.array([1.0, -1.0, 0.0], dtype=np.float64)

    result = _run_extension_with_zero_padding(
        image,
        kernel,
        nan_interpolate=True,
    )
    expected = np.array([5.0, 20.0, 20.0, 30.0, 40.0], dtype=np.float64)

    npt.assert_allclose(result, expected)


@pytest.mark.parametrize(
    ("image", "kernel"),
    [
        pytest.param(
            np.linspace(1.0, 9.0, 9, dtype=np.float64),
            np.array([1.0, 2.0, 1.0], dtype=np.float64),
            id="1d",
        ),
        pytest.param(
            np.arange(1.0, 101.0, dtype=np.float64).reshape(10, 10),
            np.ones((3, 3), dtype=np.float64),
            id="2d",
        ),
        pytest.param(
            np.arange(1.0, 513.0, dtype=np.float64).reshape(8, 8, 8),
            np.ones((3, 3, 3), dtype=np.float64),
            id="3d",
        ),
    ],
)
def test_multi_threaded_results_match_single_threaded(image, kernel):
    """Changing ``n_threads`` should not change the numerical result."""
    single_thread = _run_extension_with_zero_padding(
        image,
        kernel,
        n_threads=1,
    )
    multi_thread = _run_extension_with_zero_padding(
        image,
        kernel,
        n_threads=4,
    )

    npt.assert_allclose(single_thread, multi_thread)

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools
from contextlib import nullcontext

import numpy as np
import numpy.ma as ma
import pytest
from numpy.testing import (
    assert_allclose,
    assert_array_almost_equal,
    assert_array_almost_equal_nulp,
)

from astropy import units as u
from astropy.convolution.convolve import convolve, convolve_fft
from astropy.convolution.kernels import Gaussian2DKernel
from astropy.utils.compat.optional_deps import HAS_PANDAS, HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning

VALID_DTYPES = (">f4", "<f4", ">f8", "<f8")
VALID_DTYPE_MATRIX = list(itertools.product(VALID_DTYPES, VALID_DTYPES))

BOUNDARY_OPTIONS = [None, "fill", "wrap", "extend"]
NANHANDLING_OPTIONS = ["interpolate", "fill"]
NORMALIZE_OPTIONS = [True, False]
PRESERVE_NAN_OPTIONS = [True, False]

BOUNDARIES_AND_CONVOLUTIONS = list(
    zip(itertools.cycle((convolve,)), BOUNDARY_OPTIONS)
) + [(convolve_fft, "wrap"), (convolve_fft, "fill")]


class TestConvolve1D:
    def test_list(self):
        """
        Test that convolve works correctly when inputs are lists
        """

        x = [1, 4, 5, 6, 5, 7, 8]
        y = [0.2, 0.6, 0.2]
        z = convolve(x, y, boundary=None)
        assert_array_almost_equal_nulp(
            z, np.array([0.0, 3.6, 5.0, 5.6, 5.6, 6.8, 0.0]), 10
        )

    def test_tuple(self):
        """
        Test that convolve works correctly when inputs are tuples
        """

        x = (1, 4, 5, 6, 5, 7, 8)
        y = (0.2, 0.6, 0.2)
        z = convolve(x, y, boundary=None)
        assert_array_almost_equal_nulp(
            z, np.array([0.0, 3.6, 5.0, 5.6, 5.6, 6.8, 0.0]), 10
        )

    @pytest.mark.parametrize(
        ("boundary", "nan_treatment", "normalize_kernel", "preserve_nan", "dtype"),
        itertools.product(
            BOUNDARY_OPTIONS,
            NANHANDLING_OPTIONS,
            NORMALIZE_OPTIONS,
            PRESERVE_NAN_OPTIONS,
            VALID_DTYPES,
        ),
    )
    def test_quantity(
        self, boundary, nan_treatment, normalize_kernel, preserve_nan, dtype
    ):
        """
        Test that convolve works correctly when input array is a Quantity
        """

        x = np.array([1, 4, 5, 6, 5, 7, 8], dtype=dtype) * u.ph
        y = np.array([0.2, 0.6, 0.2], dtype=dtype)
        z = convolve(
            x,
            y,
            boundary=boundary,
            nan_treatment=nan_treatment,
            normalize_kernel=normalize_kernel,
            preserve_nan=preserve_nan,
        )

        assert x.unit == z.unit

    @pytest.mark.parametrize(
        ("boundary", "nan_treatment", "normalize_kernel", "preserve_nan", "dtype"),
        itertools.product(
            BOUNDARY_OPTIONS,
            NANHANDLING_OPTIONS,
            NORMALIZE_OPTIONS,
            PRESERVE_NAN_OPTIONS,
            VALID_DTYPES,
        ),
    )
    def test_input_unmodified(
        self, boundary, nan_treatment, normalize_kernel, preserve_nan, dtype
    ):
        """
        Test that convolve works correctly when inputs are lists
        """

        array = [1.0, 4.0, 5.0, 6.0, 5.0, 7.0, 8.0]
        kernel = [0.2, 0.6, 0.2]
        x = np.array(array, dtype=dtype)
        y = np.array(kernel, dtype=dtype)

        # Make pseudoimmutable
        x.flags.writeable = False
        y.flags.writeable = False

        convolve(
            x,
            y,
            boundary=boundary,
            nan_treatment=nan_treatment,
            normalize_kernel=normalize_kernel,
            preserve_nan=preserve_nan,
        )

        assert np.all(np.array(array, dtype=dtype) == x)
        assert np.all(np.array(kernel, dtype=dtype) == y)

    @pytest.mark.parametrize(
        ("boundary", "nan_treatment", "normalize_kernel", "preserve_nan", "dtype"),
        itertools.product(
            BOUNDARY_OPTIONS,
            NANHANDLING_OPTIONS,
            NORMALIZE_OPTIONS,
            PRESERVE_NAN_OPTIONS,
            VALID_DTYPES,
        ),
    )
    def test_input_unmodified_with_nan(
        self, boundary, nan_treatment, normalize_kernel, preserve_nan, dtype
    ):
        """
        Test that convolve doesn't modify the input data
        """

        array = [1.0, 4.0, 5.0, np.nan, 5.0, 7.0, 8.0]
        kernel = [0.2, 0.6, 0.2]
        x = np.array(array, dtype=dtype)
        y = np.array(kernel, dtype=dtype)

        # Make pseudoimmutable
        x.flags.writeable = False
        y.flags.writeable = False

        # make copies for post call comparison
        x_copy = x.copy()
        y_copy = y.copy()

        convolve(
            x,
            y,
            boundary=boundary,
            nan_treatment=nan_treatment,
            normalize_kernel=normalize_kernel,
            preserve_nan=preserve_nan,
        )

        # ( NaN == NaN ) = False
        # Only compare non NaN values for canonical equivalence
        # and then check NaN explicitly with np.isnan()
        array_is_nan = np.isnan(array)
        kernel_is_nan = np.isnan(kernel)
        array_not_nan = ~array_is_nan
        kernel_not_nan = ~kernel_is_nan
        assert np.all(x_copy[array_not_nan] == x[array_not_nan])
        assert np.all(y_copy[kernel_not_nan] == y[kernel_not_nan])
        assert np.all(np.isnan(x[array_is_nan]))
        assert np.all(np.isnan(y[kernel_is_nan]))

    @pytest.mark.parametrize(("dtype_array", "dtype_kernel"), VALID_DTYPE_MATRIX)
    def test_dtype(self, dtype_array, dtype_kernel):
        """
        Test that 32- and 64-bit floats are correctly handled
        """

        x = np.array([1.0, 2.0, 3.0], dtype=dtype_array)

        y = np.array([0.0, 1.0, 0.0], dtype=dtype_kernel)

        z = convolve(x, y)

        assert x.dtype == z.dtype

    @pytest.mark.parametrize(
        (
            "convfunc",
            "boundary",
        ),
        BOUNDARIES_AND_CONVOLUTIONS,
    )
    def test_unity_1_none(self, boundary, convfunc):
        """
        Test that a unit kernel with a single element returns the same array
        """

        x = np.array([1.0, 2.0, 3.0], dtype=">f8")

        y = np.array([1.0], dtype=">f8")

        z = convfunc(x, y, boundary=boundary)

        np.testing.assert_allclose(z, x)

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_unity_3(self, boundary):
        """
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None).
        """

        x = np.array([1.0, 2.0, 3.0], dtype=">f8")

        y = np.array([0.0, 1.0, 0.0], dtype=">f8")

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(z == np.array([0.0, 2.0, 0.0], dtype=">f8"))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_uniform_3(self, boundary):
        """
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements
        """

        x = np.array([1.0, 0.0, 3.0], dtype=">f8")

        y = np.array([1.0, 1.0, 1.0], dtype=">f8")

        z = convolve(x, y, boundary=boundary, normalize_kernel=False)

        if boundary is None:
            assert np.all(z == np.array([0.0, 4.0, 0.0], dtype=">f8"))
        elif boundary == "fill":
            assert np.all(z == np.array([1.0, 4.0, 3.0], dtype=">f8"))
        elif boundary == "wrap":
            assert np.all(z == np.array([4.0, 4.0, 4.0], dtype=">f8"))
        else:
            assert np.all(z == np.array([2.0, 4.0, 6.0], dtype=">f8"))

    @pytest.mark.parametrize(
        ("boundary", "nan_treatment", "normalize_kernel", "preserve_nan"),
        itertools.product(
            BOUNDARY_OPTIONS,
            NANHANDLING_OPTIONS,
            NORMALIZE_OPTIONS,
            PRESERVE_NAN_OPTIONS,
        ),
    )
    def test_unity_3_withnan(
        self, boundary, nan_treatment, normalize_kernel, preserve_nan
    ):
        """
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        """

        x = np.array([1.0, np.nan, 3.0], dtype=">f8")

        y = np.array([0.0, 1.0, 0.0], dtype=">f8")

        # Since the kernel is actually only one pixel wide (because of the
        # zeros) the NaN value doesn't get interpolated over so a warning is
        # expected.
        if nan_treatment == "interpolate" and not preserve_nan:
            ctx = pytest.warns(
                AstropyUserWarning,
                match="nan_treatment='interpolate', however, NaN values detected",
            )
        else:
            ctx = nullcontext()

        with ctx:
            z = convolve(
                x,
                y,
                boundary=boundary,
                nan_treatment=nan_treatment,
                normalize_kernel=normalize_kernel,
                preserve_nan=preserve_nan,
            )

        if preserve_nan:
            assert np.isnan(z[1])

        x = np.nan_to_num(z)
        z = np.nan_to_num(z)

        if boundary is None:
            assert np.all(z == np.array([0.0, 0.0, 0.0], dtype=">f8"))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize(
        ("boundary", "nan_treatment", "normalize_kernel", "preserve_nan"),
        itertools.product(
            BOUNDARY_OPTIONS,
            NANHANDLING_OPTIONS,
            NORMALIZE_OPTIONS,
            PRESERVE_NAN_OPTIONS,
        ),
    )
    def test_uniform_3_withnan(
        self, boundary, nan_treatment, normalize_kernel, preserve_nan
    ):
        """
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements. This version includes a NaN
        value in the original array.
        """

        x = np.array([1.0, np.nan, 3.0], dtype=">f8")

        y = np.array([1.0, 1.0, 1.0], dtype=">f8")

        z = convolve(
            x,
            y,
            boundary=boundary,
            nan_treatment=nan_treatment,
            normalize_kernel=normalize_kernel,
            preserve_nan=preserve_nan,
        )

        if preserve_nan:
            assert np.isnan(z[1])

        z = np.nan_to_num(z)

        # boundary, nan_treatment, normalize_kernel
        rslt = {
            (None, "interpolate", True): [0, 2, 0],
            (None, "interpolate", False): [0, 6, 0],
            (None, "fill", True): [0, 4 / 3.0, 0],
            (None, "fill", False): [0, 4, 0],
            ("fill", "interpolate", True): [1 / 2.0, 2, 3 / 2.0],
            ("fill", "interpolate", False): [3 / 2.0, 6, 9 / 2.0],
            ("fill", "fill", True): [1 / 3.0, 4 / 3.0, 3 / 3.0],
            ("fill", "fill", False): [1, 4, 3],
            ("wrap", "interpolate", True): [2, 2, 2],
            ("wrap", "interpolate", False): [6, 6, 6],
            ("wrap", "fill", True): [4 / 3.0, 4 / 3.0, 4 / 3.0],
            ("wrap", "fill", False): [4, 4, 4],
            ("extend", "interpolate", True): [1, 2, 3],
            ("extend", "interpolate", False): [3, 6, 9],
            ("extend", "fill", True): [2 / 3.0, 4 / 3.0, 6 / 3.0],
            ("extend", "fill", False): [2, 4, 6],
        }[boundary, nan_treatment, normalize_kernel]
        if preserve_nan:
            rslt[1] = 0

        assert_array_almost_equal_nulp(z, np.array(rslt, dtype=">f8"), 10)

    @pytest.mark.parametrize(
        ("boundary", "normalize_kernel"),
        itertools.product(BOUNDARY_OPTIONS, NORMALIZE_OPTIONS),
    )
    def test_zero_sum_kernel(self, boundary, normalize_kernel):
        """
        Test that convolve works correctly with zero sum kernels.
        """

        if normalize_kernel:
            pytest.xfail("You can't normalize by a zero sum kernel")

        x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        y = [-1, -1, -1, -1, 8, -1, -1, -1, -1]
        assert np.isclose(sum(y), 0, atol=1e-8)

        z = convolve(x, y, boundary=boundary, normalize_kernel=normalize_kernel)

        # boundary, normalize_kernel == False
        rslt = {
            (None): [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "fill": [-6.0, -3.0, -1.0, 0.0, 0.0, 10.0, 21.0, 33.0, 46.0],
            "wrap": [-36.0, -27.0, -18.0, -9.0, 0.0, 9.0, 18.0, 27.0, 36.0],
            "extend": [-10.0, -6.0, -3.0, -1.0, 0.0, 1.0, 3.0, 6.0, 10.0],
        }[boundary]

        assert_array_almost_equal_nulp(z, np.array(rslt, dtype=">f8"), 10)

    @pytest.mark.parametrize(
        ("boundary", "normalize_kernel"),
        itertools.product(BOUNDARY_OPTIONS, NORMALIZE_OPTIONS),
    )
    def test_int_masked_kernel(self, boundary, normalize_kernel):
        """
        Test that convolve works correctly with integer masked kernels.
        """

        if normalize_kernel:
            pytest.xfail("You can't normalize by a zero sum kernel")

        x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        y = ma.array(
            [-1, -1, -1, -1, 8, -1, -1, -1, -1],
            mask=[1, 0, 0, 0, 0, 0, 0, 0, 0],
            fill_value=0.0,
        )

        z = convolve(x, y, boundary=boundary, normalize_kernel=normalize_kernel)

        # boundary, normalize_kernel == False
        rslt = {
            (None): [0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0],
            "fill": [-1.0, 3.0, 6.0, 8.0, 9.0, 10.0, 21.0, 33.0, 46.0],
            "wrap": [-31.0, -21.0, -11.0, -1.0, 9.0, 10.0, 20.0, 30.0, 40.0],
            "extend": [-5.0, 0.0, 4.0, 7.0, 9.0, 10.0, 12.0, 15.0, 19.0],
        }[boundary]

        assert_array_almost_equal_nulp(z, np.array(rslt, dtype=">f8"), 10)

    @pytest.mark.parametrize("preserve_nan", PRESERVE_NAN_OPTIONS)
    def test_int_masked_array(self, preserve_nan):
        """
        Test that convolve works correctly with integer masked arrays.
        """

        x = ma.array([3, 5, 7, 11, 13], mask=[0, 0, 1, 0, 0], fill_value=0.0)
        y = np.array([1.0, 1.0, 1.0], dtype=">f8")

        z = convolve(x, y, preserve_nan=preserve_nan)

        if preserve_nan:
            assert np.isnan(z[2])
            z[2] = 8

        assert_array_almost_equal_nulp(z, (8 / 3.0, 4, 8, 12, 8), 10)


class TestConvolve2D:
    def test_list(self):
        """
        Test that convolve works correctly when inputs are lists
        """
        x = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]

        z = convolve(x, x, boundary="fill", fill_value=1, normalize_kernel=True)
        assert_array_almost_equal_nulp(z, x, 10)
        z = convolve(x, x, boundary="fill", fill_value=1, normalize_kernel=False)
        assert_array_almost_equal_nulp(z, np.array(x, float) * 9, 10)

    @pytest.mark.parametrize(("dtype_array", "dtype_kernel"), VALID_DTYPE_MATRIX)
    def test_dtype(self, dtype_array, dtype_kernel):
        """
        Test that 32- and 64-bit floats are correctly handled
        """

        x = np.array(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], dtype=dtype_array
        )

        y = np.array(
            [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]], dtype=dtype_kernel
        )

        z = convolve(x, y)

        assert x.dtype == z.dtype

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_unity_1x1_none(self, boundary):
        """
        Test that a 1x1 unit kernel returns the same array
        """

        x = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], dtype=">f8")

        y = np.array([[1.0]], dtype=">f8")

        z = convolve(x, y, boundary=boundary)

        assert np.all(z == x)

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_unity_3x3(self, boundary):
        """
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None).
        """

        x = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], dtype=">f8")

        y = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]], dtype=">f8")

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(
                z
                == np.array(
                    [[0.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 0.0]], dtype=">f8"
                )
            )
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_uniform_3x3(self, boundary):
        """
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel.
        """

        x = np.array([[0.0, 0.0, 3.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]], dtype=">f8")

        y = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], dtype=">f8")

        z = convolve(x, y, boundary=boundary, normalize_kernel=False)

        if boundary is None:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[0.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 0.0]], dtype=">f8"
                ),
                10,
            )
        elif boundary == "fill":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[1.0, 4.0, 3.0], [3.0, 6.0, 5.0], [3.0, 3.0, 2.0]], dtype=">f8"
                ),
                10,
            )
        elif boundary == "wrap":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[6.0, 6.0, 6.0], [6.0, 6.0, 6.0], [6.0, 6.0, 6.0]], dtype=">f8"
                ),
                10,
            )
        else:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[2.0, 7.0, 12.0], [4.0, 6.0, 8.0], [6.0, 5.0, 4.0]], dtype=">f8"
                ),
                10,
            )

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_unity_3x3_withnan(self, boundary):
        """
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None). This version includes a NaN value in the original
        array.
        """

        x = np.array(
            [[1.0, 2.0, 3.0], [4.0, np.nan, 6.0], [7.0, 8.0, 9.0]], dtype=">f8"
        )

        y = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]], dtype=">f8")

        z = convolve(x, y, boundary=boundary, nan_treatment="fill", preserve_nan=True)

        assert np.isnan(z[1, 1])
        x = np.nan_to_num(z)
        z = np.nan_to_num(z)

        if boundary is None:
            assert np.all(
                z
                == np.array(
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=">f8"
                )
            )
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_uniform_3x3_withnanfilled(self, boundary):
        """
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        """

        x = np.array(
            [[0.0, 0.0, 4.0], [1.0, np.nan, 0.0], [0.0, 3.0, 0.0]], dtype=">f8"
        )

        y = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], dtype=">f8")

        z = convolve(
            x, y, boundary=boundary, nan_treatment="fill", normalize_kernel=False
        )

        if boundary is None:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[0.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 0.0]], dtype=">f8"
                ),
                10,
            )
        elif boundary == "fill":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[1.0, 5.0, 4.0], [4.0, 8.0, 7.0], [4.0, 4.0, 3.0]], dtype=">f8"
                ),
                10,
            )
        elif boundary == "wrap":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[8.0, 8.0, 8.0], [8.0, 8.0, 8.0], [8.0, 8.0, 8.0]], dtype=">f8"
                ),
                10,
            )
        elif boundary == "extend":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[2.0, 9.0, 16.0], [5.0, 8.0, 11.0], [8.0, 7.0, 6.0]], dtype=">f8"
                ),
                10,
            )
        else:
            raise ValueError("Invalid boundary specification")

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_uniform_3x3_withnaninterped(self, boundary):
        """
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        """

        x = np.array(
            [[0.0, 0.0, 4.0], [1.0, np.nan, 0.0], [0.0, 3.0, 0.0]], dtype=">f8"
        )

        y = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], dtype=">f8")

        z = convolve(
            x, y, boundary=boundary, nan_treatment="interpolate", normalize_kernel=True
        )

        if boundary is None:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]], dtype=">f8"
                ),
                10,
            )
        elif boundary == "fill":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [1.0 / 8, 5.0 / 8, 4.0 / 8],
                        [4.0 / 8, 8.0 / 8, 7.0 / 8],
                        [4.0 / 8, 4.0 / 8, 3.0 / 8],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        elif boundary == "wrap":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], dtype=">f8"
                ),
                10,
            )
        elif boundary == "extend":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [2.0 / 8, 9.0 / 8, 16.0 / 8],
                        [5.0 / 8, 8.0 / 8, 11.0 / 8],
                        [8.0 / 8, 7.0 / 8, 6.0 / 8],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        else:
            raise ValueError("Invalid boundary specification")

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_non_normalized_kernel_2D(self, boundary):
        x = np.array([[0.0, 0.0, 4.0], [1.0, 2.0, 0.0], [0.0, 3.0, 0.0]], dtype="float")

        y = np.array(
            [[1.0, -1.0, 1.0], [-1.0, 0.0, -1.0], [1.0, -1.0, 1.0]], dtype="float"
        )

        z = convolve(
            x, y, boundary=boundary, nan_treatment="fill", normalize_kernel=False
        )

        if boundary is None:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype="float"
                ),
                10,
            )
        elif boundary == "fill":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[1.0, -5.0, 2.0], [1.0, 0.0, -3.0], [-2.0, -1.0, -1.0]],
                    dtype="float",
                ),
                10,
            )
        elif boundary == "wrap":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[0.0, -8.0, 6.0], [5.0, 0.0, -4.0], [2.0, 3.0, -4.0]],
                    dtype="float",
                ),
                10,
            )
        elif boundary == "extend":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [[2.0, -1.0, -2.0], [0.0, 0.0, 1.0], [2.0, -4.0, 2.0]],
                    dtype="float",
                ),
                10,
            )
        else:
            raise ValueError("Invalid boundary specification")


class TestConvolve3D:
    def test_list(self):
        """
        Test that convolve works correctly when inputs are lists
        """
        x = [
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
        ]

        z = convolve(x, x, boundary="fill", fill_value=1, normalize_kernel=False)
        assert_array_almost_equal_nulp(z / 27, x, 10)

    @pytest.mark.parametrize(("dtype_array", "dtype_kernel"), VALID_DTYPE_MATRIX)
    def test_dtype(self, dtype_array, dtype_kernel):
        """
        Test that 32- and 64-bit floats are correctly handled
        """

        x = np.array(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], dtype=dtype_array
        )

        y = np.array(
            [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]], dtype=dtype_kernel
        )

        z = convolve(x, y)

        assert x.dtype == z.dtype

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_unity_1x1x1_none(self, boundary):
        """
        Test that a 1x1x1 unit kernel returns the same array
        """

        x = np.array(
            [
                [[1.0, 2.0, 1.0], [2.0, 3.0, 1.0], [3.0, 2.0, 5.0]],
                [[4.0, 3.0, 1.0], [5.0, 0.0, 2.0], [6.0, 1.0, 1.0]],
                [[7.0, 0.0, 2.0], [8.0, 2.0, 3.0], [9.0, 2.0, 2.0]],
            ],
            dtype=">f8",
        )

        y = np.array([[[1.0]]], dtype=">f8")

        z = convolve(x, y, boundary=boundary)

        assert np.all(z == x)

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_unity_3x3x3(self, boundary):
        """
        Test that a 3x3x3 unit kernel returns the same array (except when
        boundary is None).
        """

        x = np.array(
            [
                [[1.0, 2.0, 1.0], [2.0, 3.0, 1.0], [3.0, 2.0, 5.0]],
                [[4.0, 3.0, 1.0], [5.0, 3.0, 2.0], [6.0, 1.0, 1.0]],
                [[7.0, 0.0, 2.0], [8.0, 2.0, 3.0], [9.0, 2.0, 2.0]],
            ],
            dtype=">f8",
        )

        y = np.zeros((3, 3, 3), dtype=">f8")
        y[1, 1, 1] = 1.0

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(
                z
                == np.array(
                    [
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    ],
                    dtype=">f8",
                )
            )
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_uniform_3x3x3(self, boundary):
        """
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel.
        """

        x = np.array(
            [
                [[1.0, 2.0, 1.0], [2.0, 3.0, 1.0], [3.0, 2.0, 5.0]],
                [[4.0, 3.0, 1.0], [5.0, 3.0, 2.0], [6.0, 1.0, 1.0]],
                [[7.0, 0.0, 2.0], [8.0, 2.0, 3.0], [9.0, 2.0, 2.0]],
            ],
            dtype=">f8",
        )

        y = np.ones((3, 3, 3), dtype=">f8")

        z = convolve(x, y, boundary=boundary, normalize_kernel=False)

        if boundary is None:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 81.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        elif boundary == "fill":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[23.0, 28.0, 16.0], [35.0, 46.0, 25.0], [25.0, 34.0, 18.0]],
                        [[40.0, 50.0, 23.0], [63.0, 81.0, 36.0], [46.0, 60.0, 27.0]],
                        [[32.0, 40.0, 16.0], [50.0, 61.0, 22.0], [36.0, 44.0, 16.0]],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        elif boundary == "wrap":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[81.0, 81.0, 81.0], [81.0, 81.0, 81.0], [81.0, 81.0, 81.0]],
                        [[81.0, 81.0, 81.0], [81.0, 81.0, 81.0], [81.0, 81.0, 81.0]],
                        [[81.0, 81.0, 81.0], [81.0, 81.0, 81.0], [81.0, 81.0, 81.0]],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        else:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[65.0, 54.0, 43.0], [75.0, 66.0, 57.0], [85.0, 78.0, 71.0]],
                        [[96.0, 71.0, 46.0], [108.0, 81.0, 54.0], [120.0, 91.0, 62.0]],
                        [
                            [127.0, 88.0, 49.0],
                            [141.0, 96.0, 51.0],
                            [155.0, 104.0, 53.0],
                        ],
                    ],
                    dtype=">f8",
                ),
                10,
            )

    @pytest.mark.parametrize(
        ("boundary", "nan_treatment"),
        itertools.product(BOUNDARY_OPTIONS, NANHANDLING_OPTIONS),
    )
    def test_unity_3x3x3_withnan(self, boundary, nan_treatment):
        """
        Test that a 3x3x3 unit kernel returns the same array (except when
        boundary is None). This version includes a NaN value in the original
        array.
        """

        x = np.array(
            [
                [[1.0, 2.0, 1.0], [2.0, 3.0, 1.0], [3.0, 2.0, 5.0]],
                [[4.0, 3.0, 1.0], [5.0, np.nan, 2.0], [6.0, 1.0, 1.0]],
                [[7.0, 0.0, 2.0], [8.0, 2.0, 3.0], [9.0, 2.0, 2.0]],
            ],
            dtype=">f8",
        )

        y = np.zeros((3, 3, 3), dtype=">f8")
        y[1, 1, 1] = 1.0

        z = convolve(
            x, y, boundary=boundary, nan_treatment=nan_treatment, preserve_nan=True
        )

        assert np.isnan(z[1, 1, 1])
        x = np.nan_to_num(z)
        z = np.nan_to_num(z)

        if boundary is None:
            assert np.all(
                z
                == np.array(
                    [
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    ],
                    dtype=">f8",
                )
            )
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_uniform_3x3x3_withnan_filled(self, boundary):
        """
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        """

        x = np.array(
            [
                [[1.0, 2.0, 1.0], [2.0, 3.0, 1.0], [3.0, 2.0, 5.0]],
                [[4.0, 3.0, 1.0], [5.0, np.nan, 2.0], [6.0, 1.0, 1.0]],
                [[7.0, 0.0, 2.0], [8.0, 2.0, 3.0], [9.0, 2.0, 2.0]],
            ],
            dtype=">f8",
        )

        y = np.ones((3, 3, 3), dtype=">f8")

        z = convolve(
            x, y, boundary=boundary, nan_treatment="fill", normalize_kernel=False
        )

        if boundary is None:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 78.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        elif boundary == "fill":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[20.0, 25.0, 13.0], [32.0, 43.0, 22.0], [22.0, 31.0, 15.0]],
                        [[37.0, 47.0, 20.0], [60.0, 78.0, 33.0], [43.0, 57.0, 24.0]],
                        [[29.0, 37.0, 13.0], [47.0, 58.0, 19.0], [33.0, 41.0, 13.0]],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        elif boundary == "wrap":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[78.0, 78.0, 78.0], [78.0, 78.0, 78.0], [78.0, 78.0, 78.0]],
                        [[78.0, 78.0, 78.0], [78.0, 78.0, 78.0], [78.0, 78.0, 78.0]],
                        [[78.0, 78.0, 78.0], [78.0, 78.0, 78.0], [78.0, 78.0, 78.0]],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        elif boundary == "extend":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[62.0, 51.0, 40.0], [72.0, 63.0, 54.0], [82.0, 75.0, 68.0]],
                        [[93.0, 68.0, 43.0], [105.0, 78.0, 51.0], [117.0, 88.0, 59.0]],
                        [
                            [124.0, 85.0, 46.0],
                            [138.0, 93.0, 48.0],
                            [152.0, 101.0, 50.0],
                        ],
                    ],
                    dtype=">f8",
                ),
                10,
            )
        else:
            raise ValueError("Invalid Boundary Option")

    @pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
    def test_uniform_3x3x3_withnan_interped(self, boundary):
        """
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        """

        x = np.array(
            [
                [[1.0, 2.0, 1.0], [2.0, 3.0, 1.0], [3.0, 2.0, 5.0]],
                [[4.0, 3.0, 1.0], [5.0, np.nan, 2.0], [6.0, 1.0, 1.0]],
                [[7.0, 0.0, 2.0], [8.0, 2.0, 3.0], [9.0, 2.0, 2.0]],
            ],
            dtype=">f8",
        )

        y = np.ones((3, 3, 3), dtype=">f8")

        z = convolve(
            x, y, boundary=boundary, nan_treatment="interpolate", normalize_kernel=True
        )

        kernsum = y.sum() - 1  # one nan is missing
        mid = x[np.isfinite(x)].sum() / kernsum

        if boundary is None:
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 78.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    ],
                    dtype=">f8",
                )
                / kernsum,
                10,
            )
        elif boundary == "fill":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[20.0, 25.0, 13.0], [32.0, 43.0, 22.0], [22.0, 31.0, 15.0]],
                        [[37.0, 47.0, 20.0], [60.0, 78.0, 33.0], [43.0, 57.0, 24.0]],
                        [[29.0, 37.0, 13.0], [47.0, 58.0, 19.0], [33.0, 41.0, 13.0]],
                    ],
                    dtype=">f8",
                )
                / kernsum,
                10,
            )
        elif boundary == "wrap":
            assert_array_almost_equal_nulp(z, np.tile(mid.astype(">f8"), [3, 3, 3]), 10)
        elif boundary == "extend":
            assert_array_almost_equal_nulp(
                z,
                np.array(
                    [
                        [[62.0, 51.0, 40.0], [72.0, 63.0, 54.0], [82.0, 75.0, 68.0]],
                        [[93.0, 68.0, 43.0], [105.0, 78.0, 51.0], [117.0, 88.0, 59.0]],
                        [
                            [124.0, 85.0, 46.0],
                            [138.0, 93.0, 48.0],
                            [152.0, 101.0, 50.0],
                        ],
                    ],
                    dtype=">f8",
                )
                / kernsum,
                10,
            )
        else:
            raise ValueError("Invalid Boundary Option")


@pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
def test_asymmetric_kernel(boundary):
    """
    Regression test for #6264: make sure that asymmetric convolution
    functions go the right direction
    """

    x = np.array([3.0, 0.0, 1.0], dtype=">f8")

    y = np.array([1, 2, 3], dtype=">f8")

    z = convolve(x, y, boundary=boundary, normalize_kernel=False)

    if boundary == "fill":
        assert_array_almost_equal_nulp(z, np.array([6.0, 10.0, 2.0], dtype="float"), 10)
    elif boundary is None:
        assert_array_almost_equal_nulp(z, np.array([0.0, 10.0, 0.0], dtype="float"), 10)
    elif boundary == "extend":
        assert_array_almost_equal_nulp(
            z, np.array([15.0, 10.0, 3.0], dtype="float"), 10
        )
    elif boundary == "wrap":
        assert_array_almost_equal_nulp(z, np.array([9.0, 10.0, 5.0], dtype="float"), 10)


@pytest.mark.parametrize("ndims", (1, 2, 3))
def test_convolution_consistency(ndims):
    np.random.seed(0)
    array = np.random.randn(*([3] * ndims))
    np.random.seed(0)
    kernel = np.random.rand(*([3] * ndims))

    conv_f = convolve_fft(array, kernel, boundary="fill")
    conv_d = convolve(array, kernel, boundary="fill")

    assert_array_almost_equal_nulp(conv_f, conv_d, 30)


def test_astropy_convolution_against_numpy():
    x = np.array([1, 2, 3])
    y = np.array([5, 4, 3, 2, 1])

    assert_array_almost_equal(
        np.convolve(y, x, "same"), convolve(y, x, normalize_kernel=False)
    )
    assert_array_almost_equal(
        np.convolve(y, x, "same"), convolve_fft(y, x, normalize_kernel=False)
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
def test_astropy_convolution_against_scipy():
    from scipy.signal import fftconvolve

    x = np.array([1, 2, 3])
    y = np.array([5, 4, 3, 2, 1])

    assert_array_almost_equal(
        fftconvolve(y, x, "same"), convolve(y, x, normalize_kernel=False)
    )
    assert_array_almost_equal(
        fftconvolve(y, x, "same"), convolve_fft(y, x, normalize_kernel=False)
    )


@pytest.mark.skipif(not HAS_PANDAS, reason="Requires pandas")
def test_regression_6099():
    import pandas

    wave = np.array(np.linspace(5000, 5100, 10))
    boxcar = 3
    nonseries_result = convolve(wave, np.ones((boxcar,)) / boxcar)

    wave_series = pandas.Series(wave)
    series_result = convolve(wave_series, np.ones((boxcar,)) / boxcar)

    assert_array_almost_equal(nonseries_result, series_result)


def test_invalid_array_convolve():
    kernel = np.ones(3) / 3.0

    with pytest.raises(TypeError):
        convolve("glork", kernel)


@pytest.mark.parametrize("boundary", BOUNDARY_OPTIONS)
def test_non_square_kernel_asymmetric(boundary):
    # Regression test for a bug that occurred when using non-square kernels in
    # 2D when using boundary=None
    kernel = np.array([[1, 2, 3, 2, 1], [0, 1, 2, 1, 0], [0, 0, 0, 0, 0]])
    image = np.zeros((13, 13))
    image[6, 6] = 1
    result = convolve(image, kernel, normalize_kernel=False, boundary=boundary)
    assert_allclose(result[5:8, 4:9], kernel)


@pytest.mark.parametrize(
    ("boundary", "normalize_kernel"),
    itertools.product(BOUNDARY_OPTIONS, NORMALIZE_OPTIONS),
)
def test_uninterpolated_nan_regions(boundary, normalize_kernel):
    # Issue #8086
    # Test NaN interpolation of contiguous NaN regions with kernels of size
    # identical and greater than that of the region of NaN values.

    # Test case: kernel.shape == NaN_region.shape
    kernel = Gaussian2DKernel(1, 5, 5)
    nan_centroid = np.full(kernel.shape, np.nan)
    image = np.pad(
        nan_centroid, pad_width=kernel.shape[0] * 2, mode="constant", constant_values=1
    )
    with pytest.warns(
        AstropyUserWarning,
        match=r"nan_treatment='interpolate', however, NaN values detected "
        r"post convolution. A contiguous region of NaN values, larger "
        r"than the kernel size, are present in the input array. "
        r"Increase the kernel size to avoid this.",
    ):
        result = convolve(
            image,
            kernel,
            boundary=boundary,
            nan_treatment="interpolate",
            normalize_kernel=normalize_kernel,
        )
        assert np.any(np.isnan(result))

    # Test case: kernel.shape > NaN_region.shape
    nan_centroid = np.full(
        (kernel.shape[0] - 1, kernel.shape[1] - 1), np.nan
    )  # 1 smaller than kerenel
    image = np.pad(
        nan_centroid, pad_width=kernel.shape[0] * 2, mode="constant", constant_values=1
    )
    result = convolve(
        image,
        kernel,
        boundary=boundary,
        nan_treatment="interpolate",
        normalize_kernel=normalize_kernel,
    )
    assert ~np.any(np.isnan(result))  # Note: negation


def test_regressiontest_issue9168():
    """
    Issue #9168 pointed out that kernels can be (unitless) quantities, which
    leads to crashes when inplace modifications are made to arrays in
    convolve/convolve_fft, so we now strip the quantity aspects off of kernels.
    """

    x = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
    )

    kernel_fwhm = 1 * u.arcsec
    pixel_size = 1 * u.arcsec

    kernel = Gaussian2DKernel(x_stddev=kernel_fwhm / pixel_size)

    convolve_fft(x, kernel, boundary="fill", fill_value=np.nan, preserve_nan=True)
    convolve(x, kernel, boundary="fill", fill_value=np.nan, preserve_nan=True)


def test_convolve_nan_zero_sum_kernel():
    with pytest.raises(
        ValueError,
        match=(
            "Setting nan_treatment='interpolate' "
            "requires the kernel to be normalized, but the "
            "input kernel has a sum close to zero. For a "
            "zero-sum kernel and data with NaNs, set "
            "nan_treatment='fill'."
        ),
    ):
        convolve([1, np.nan, 3], [-1, 2, -1], normalize_kernel=False)

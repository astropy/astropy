# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools

import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal

from astropy.convolution.convolve import convolve, convolve_fft
from astropy.convolution.kernels import (AiryDisk2DKernel, Box1DKernel, Box2DKernel, CustomKernel,
                                         Gaussian1DKernel, Gaussian2DKernel, Kernel1D, Kernel2D,
                                         Model1DKernel, Model2DKernel, RickerWavelet1DKernel,
                                         RickerWavelet2DKernel, Ring2DKernel, Tophat2DKernel,
                                         Trapezoid1DKernel, TrapezoidDisk2DKernel)
from astropy.convolution.utils import KernelSizeError
from astropy.modeling.models import Box2D, Gaussian1D, Gaussian2D
from astropy.utils.compat.optional_deps import HAS_SCIPY  # noqa
from astropy.utils.exceptions import AstropyUserWarning

WIDTHS_ODD = [3, 5, 7, 9]
WIDTHS_EVEN = [2, 4, 8, 16]
MODES = ['center', 'linear_interp', 'oversample', 'integrate']
KERNEL_TYPES = [Gaussian1DKernel, Gaussian2DKernel,
                Box1DKernel, Box2DKernel,
                Trapezoid1DKernel, TrapezoidDisk2DKernel,
                RickerWavelet1DKernel, Tophat2DKernel, AiryDisk2DKernel,
                Ring2DKernel]


NUMS = [1, 1., np.float32(1.), np.float64(1.)]


# Test data
delta_pulse_1D = np.zeros(81)
delta_pulse_1D[40] = 1

delta_pulse_2D = np.zeros((81, 81))
delta_pulse_2D[40, 40] = 1

random_data_1D = np.random.rand(61)
random_data_2D = np.random.rand(61, 61)


class TestKernels:
    """
    Test class for the built-in convolution kernels.
    """

    @pytest.mark.skipif('not HAS_SCIPY')
    @pytest.mark.parametrize(('width'), WIDTHS_ODD)
    def test_scipy_filter_gaussian(self, width):
        """
        Test GaussianKernel against SciPy ndimage gaussian filter.
        """
        from scipy.ndimage import gaussian_filter

        gauss_kernel_1D = Gaussian1DKernel(width)
        gauss_kernel_1D.normalize()
        gauss_kernel_2D = Gaussian2DKernel(width)
        gauss_kernel_2D.normalize()

        astropy_1D = convolve(delta_pulse_1D, gauss_kernel_1D, boundary='fill')
        astropy_2D = convolve(delta_pulse_2D, gauss_kernel_2D, boundary='fill')

        scipy_1D = gaussian_filter(delta_pulse_1D, width)
        scipy_2D = gaussian_filter(delta_pulse_2D, width)

        assert_almost_equal(astropy_1D, scipy_1D, decimal=12)
        assert_almost_equal(astropy_2D, scipy_2D, decimal=12)

    @pytest.mark.skipif('not HAS_SCIPY')
    @pytest.mark.parametrize(('width'), WIDTHS_ODD)
    def test_scipy_filter_gaussian_laplace(self, width):
        """
        Test RickerWavelet kernels against SciPy ndimage gaussian laplace filters.
        """
        from scipy.ndimage import gaussian_laplace

        ricker_kernel_1D = RickerWavelet1DKernel(width)
        ricker_kernel_2D = RickerWavelet2DKernel(width)

        astropy_1D = convolve(delta_pulse_1D, ricker_kernel_1D, boundary='fill', normalize_kernel=False)
        astropy_2D = convolve(delta_pulse_2D, ricker_kernel_2D, boundary='fill', normalize_kernel=False)

        with pytest.raises(Exception) as exc:
            astropy_1D = convolve(delta_pulse_1D, ricker_kernel_1D, boundary='fill', normalize_kernel=True)
        assert 'sum is close to zero' in exc.value.args[0]

        with pytest.raises(Exception) as exc:
            astropy_2D = convolve(delta_pulse_2D, ricker_kernel_2D, boundary='fill', normalize_kernel=True)
        assert 'sum is close to zero' in exc.value.args[0]

        # The Laplace of Gaussian filter is an inverted Ricker Wavelet filter.
        scipy_1D = -gaussian_laplace(delta_pulse_1D, width)
        scipy_2D = -gaussian_laplace(delta_pulse_2D, width)

        # There is a slight deviation in the normalization. They differ by a
        # factor of ~1.0000284132604045. The reason is not known.
        assert_almost_equal(astropy_1D, scipy_1D, decimal=5)
        assert_almost_equal(astropy_2D, scipy_2D, decimal=5)

    @pytest.mark.parametrize(('kernel_type', 'width'), list(itertools.product(KERNEL_TYPES, WIDTHS_ODD)))
    def test_delta_data(self, kernel_type, width):
        """
        Test smoothing of an image with a single positive pixel
        """
        if kernel_type == AiryDisk2DKernel and not HAS_SCIPY:
            pytest.skip("Omitting AiryDisk2DKernel, which requires SciPy")
        if not kernel_type == Ring2DKernel:
            kernel = kernel_type(width)
        else:
            kernel = kernel_type(width, width * 0.2)

        if kernel.dimension == 1:
            c1 = convolve_fft(delta_pulse_1D, kernel, boundary='fill', normalize_kernel=False)
            c2 = convolve(delta_pulse_1D, kernel, boundary='fill', normalize_kernel=False)
            assert_almost_equal(c1, c2, decimal=12)
        else:
            c1 = convolve_fft(delta_pulse_2D, kernel, boundary='fill', normalize_kernel=False)
            c2 = convolve(delta_pulse_2D, kernel, boundary='fill', normalize_kernel=False)
            assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('kernel_type', 'width'), list(itertools.product(KERNEL_TYPES, WIDTHS_ODD)))
    def test_random_data(self, kernel_type, width):
        """
        Test smoothing of an image made of random noise
        """
        if kernel_type == AiryDisk2DKernel and not HAS_SCIPY:
            pytest.skip("Omitting AiryDisk2DKernel, which requires SciPy")
        if not kernel_type == Ring2DKernel:
            kernel = kernel_type(width)
        else:
            kernel = kernel_type(width, width * 0.2)

        if kernel.dimension == 1:
            c1 = convolve_fft(random_data_1D, kernel, boundary='fill', normalize_kernel=False)
            c2 = convolve(random_data_1D, kernel, boundary='fill', normalize_kernel=False)
            assert_almost_equal(c1, c2, decimal=12)
        else:
            c1 = convolve_fft(random_data_2D, kernel, boundary='fill', normalize_kernel=False)
            c2 = convolve(random_data_2D, kernel, boundary='fill', normalize_kernel=False)
            assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('width'), WIDTHS_ODD)
    def test_uniform_smallkernel(self, width):
        """
        Test smoothing of an image with a single positive pixel

        Instead of using kernel class, uses a simple, small kernel
        """
        kernel = np.ones([width, width])

        c2 = convolve_fft(delta_pulse_2D, kernel, boundary='fill')
        c1 = convolve(delta_pulse_2D, kernel, boundary='fill')
        assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('width'), WIDTHS_ODD)
    def test_smallkernel_vs_Box2DKernel(self, width):
        """
        Test smoothing of an image with a single positive pixel
        """
        kernel1 = np.ones([width, width]) / width ** 2
        kernel2 = Box2DKernel(width)

        c2 = convolve_fft(delta_pulse_2D, kernel2, boundary='fill')
        c1 = convolve_fft(delta_pulse_2D, kernel1, boundary='fill')

        assert_almost_equal(c1, c2, decimal=12)

    def test_convolve_1D_kernels(self):
        """
        Check if convolving two kernels with each other works correctly.
        """
        gauss_1 = Gaussian1DKernel(3)
        gauss_2 = Gaussian1DKernel(4)
        test_gauss_3 = Gaussian1DKernel(5)

        with pytest.warns(AstropyUserWarning, match=r'Both array and kernel '
                          r'are Kernel instances'):
            gauss_3 = convolve(gauss_1, gauss_2)

        assert np.all(np.abs((gauss_3 - test_gauss_3).array) < 0.01)

    def test_convolve_2D_kernels(self):
        """
        Check if convolving two kernels with each other works correctly.
        """
        gauss_1 = Gaussian2DKernel(3)
        gauss_2 = Gaussian2DKernel(4)
        test_gauss_3 = Gaussian2DKernel(5)

        with pytest.warns(AstropyUserWarning, match=r'Both array and kernel '
                          r'are Kernel instances'):
            gauss_3 = convolve(gauss_1, gauss_2)

        assert np.all(np.abs((gauss_3 - test_gauss_3).array) < 0.01)

    @pytest.mark.parametrize(('number'), NUMS)
    def test_multiply_scalar(self, number):
        """
        Check if multiplying a kernel with a scalar works correctly.
        """
        gauss = Gaussian1DKernel(3)
        gauss_new = number * gauss
        assert_almost_equal(gauss_new.array, gauss.array * number, decimal=12)

    @pytest.mark.parametrize(('number'), NUMS)
    def test_multiply_scalar_type(self, number):
        """
        Check if multiplying a kernel with a scalar works correctly.
        """
        gauss = Gaussian1DKernel(3)
        gauss_new = number * gauss
        assert type(gauss_new) is Gaussian1DKernel

    @pytest.mark.parametrize(('number'), NUMS)
    def test_rmultiply_scalar_type(self, number):
        """
        Check if multiplying a kernel with a scalar works correctly.
        """
        gauss = Gaussian1DKernel(3)
        gauss_new = gauss * number
        assert type(gauss_new) is Gaussian1DKernel

    def test_multiply_kernel1d(self):
        """Test that multiplying two 1D kernels raises an exception."""
        gauss = Gaussian1DKernel(3)
        with pytest.raises(Exception):
            gauss * gauss

    def test_multiply_kernel2d(self):
        """Test that multiplying two 2D kernels raises an exception."""
        gauss = Gaussian2DKernel(3)
        with pytest.raises(Exception):
            gauss * gauss

    def test_multiply_kernel1d_kernel2d(self):
        """
        Test that multiplying a 1D kernel with a 2D kernel raises an
        exception.
        """
        with pytest.raises(Exception):
            Gaussian1DKernel(3) * Gaussian2DKernel(3)

    def test_add_kernel_scalar(self):
        """Test that adding a scalar to a kernel raises an exception."""
        with pytest.raises(Exception):
            Gaussian1DKernel(3) + 1

    def test_model_1D_kernel(self):
        """
        Check Model1DKernel against Gaussian1Dkernel
        """
        stddev = 5.
        gauss = Gaussian1D(1. / np.sqrt(2 * np.pi * stddev**2), 0, stddev)
        model_gauss_kernel = Model1DKernel(gauss, x_size=21)
        gauss_kernel = Gaussian1DKernel(stddev, x_size=21)
        assert_almost_equal(model_gauss_kernel.array, gauss_kernel.array,
                            decimal=12)

    def test_model_2D_kernel(self):
        """
        Check Model2DKernel against Gaussian2Dkernel
        """
        stddev = 5.
        gauss = Gaussian2D(1. / (2 * np.pi * stddev**2), 0, 0, stddev, stddev)
        model_gauss_kernel = Model2DKernel(gauss, x_size=21)
        gauss_kernel = Gaussian2DKernel(stddev, x_size=21)
        assert_almost_equal(model_gauss_kernel.array, gauss_kernel.array,
                            decimal=12)

    def test_custom_1D_kernel(self):
        """
        Check CustomKernel against Box1DKernel.
        """
        # Define one dimensional array:
        array = np.ones(5)
        custom = CustomKernel(array)
        custom.normalize()
        box = Box1DKernel(5)

        c2 = convolve(delta_pulse_1D, custom, boundary='fill')
        c1 = convolve(delta_pulse_1D, box, boundary='fill')
        assert_almost_equal(c1, c2, decimal=12)

    def test_custom_2D_kernel(self):
        """
        Check CustomKernel against Box2DKernel.
        """
        # Define one dimensional array:
        array = np.ones((5, 5))
        custom = CustomKernel(array)
        custom.normalize()
        box = Box2DKernel(5)

        c2 = convolve(delta_pulse_2D, custom, boundary='fill')
        c1 = convolve(delta_pulse_2D, box, boundary='fill')
        assert_almost_equal(c1, c2, decimal=12)

    def test_custom_1D_kernel_list(self):
        """
        Check if CustomKernel works with lists.
        """
        custom = CustomKernel([1, 1, 1, 1, 1])
        assert custom.is_bool is True

    def test_custom_2D_kernel_list(self):
        """
        Check if CustomKernel works with lists.
        """
        custom = CustomKernel([[1, 1, 1],
                               [1, 1, 1],
                               [1, 1, 1]])
        assert custom.is_bool is True

    def test_custom_1D_kernel_zerosum(self):
        """
        Check if CustomKernel works when the input array/list
        sums to zero.
        """
        array = [-2, -1, 0, 1, 2]

        custom = CustomKernel(array)

        with pytest.warns(AstropyUserWarning, match=r'kernel cannot be '
                          r'normalized because it sums to zero'):
            custom.normalize()

        assert custom.truncation == 0.
        assert custom._kernel_sum == 0.

    def test_custom_2D_kernel_zerosum(self):
        """
        Check if CustomKernel works when the input array/list
        sums to zero.
        """
        array = [[0, -1, 0], [-1, 4, -1], [0, -1, 0]]

        custom = CustomKernel(array)

        with pytest.warns(AstropyUserWarning, match=r'kernel cannot be '
                          r'normalized because it sums to zero'):
            custom.normalize()

        assert custom.truncation == 0.
        assert custom._kernel_sum == 0.

    def test_custom_kernel_odd_error(self):
        """
        Check if CustomKernel raises if the array size is odd.
        """
        with pytest.raises(KernelSizeError):
            CustomKernel([1, 1, 1, 1])

    def test_add_1D_kernels(self):
        """
        Check if adding of two 1D kernels works.
        """
        box_1 = Box1DKernel(5)
        box_2 = Box1DKernel(3)
        box_3 = Box1DKernel(1)
        box_sum_1 = box_1 + box_2 + box_3
        box_sum_2 = box_2 + box_3 + box_1
        box_sum_3 = box_3 + box_1 + box_2
        ref = [1/5., 1/5. + 1/3., 1 + 1/3. + 1/5., 1/5. + 1/3., 1/5.]
        assert_almost_equal(box_sum_1.array, ref, decimal=12)
        assert_almost_equal(box_sum_2.array, ref, decimal=12)
        assert_almost_equal(box_sum_3.array, ref, decimal=12)

        # Assert that the kernels haven't changed
        assert_almost_equal(box_1.array, [0.2, 0.2, 0.2, 0.2, 0.2], decimal=12)
        assert_almost_equal(box_2.array, [1/3., 1/3., 1/3.], decimal=12)
        assert_almost_equal(box_3.array, [1], decimal=12)

    def test_add_2D_kernels(self):
        """
        Check if adding of two 1D kernels works.
        """
        box_1 = Box2DKernel(3)
        box_2 = Box2DKernel(1)
        box_sum_1 = box_1 + box_2
        box_sum_2 = box_2 + box_1
        ref = [[1 / 9., 1 / 9., 1 / 9.],
               [1 / 9., 1 + 1 / 9., 1 / 9.],
               [1 / 9., 1 / 9., 1 / 9.]]
        ref_1 = [[1 / 9., 1 / 9., 1 / 9.],
                 [1 / 9., 1 / 9., 1 / 9.],
                 [1 / 9., 1 / 9., 1 / 9.]]
        assert_almost_equal(box_2.array, [[1]], decimal=12)
        assert_almost_equal(box_1.array, ref_1, decimal=12)
        assert_almost_equal(box_sum_1.array, ref, decimal=12)
        assert_almost_equal(box_sum_2.array, ref, decimal=12)

    def test_Gaussian1DKernel_even_size(self):
        """
        Check if even size for GaussianKernel works.
        """
        gauss = Gaussian1DKernel(3, x_size=10)
        assert gauss.array.size == 10

    def test_Gaussian2DKernel_even_size(self):
        """
        Check if even size for GaussianKernel works.
        """
        gauss = Gaussian2DKernel(3, x_size=10, y_size=10)
        assert gauss.array.shape == (10, 10)

    # https://github.com/astropy/astropy/issues/3605
    def test_Gaussian2DKernel_rotated(self):
        gauss = Gaussian2DKernel(
            x_stddev=3, y_stddev=1.5, theta=0.7853981633974483,
            x_size=5, y_size=5)  # rotated 45 deg ccw
        ans = [[0.02267712, 0.02464785, 0.02029238, 0.01265463, 0.00597762],
               [0.02464785, 0.03164847, 0.03078144, 0.02267712, 0.01265463],
               [0.02029238, 0.03078144, 0.03536777, 0.03078144, 0.02029238],
               [0.01265463, 0.02267712, 0.03078144, 0.03164847, 0.02464785],
               [0.00597762, 0.01265463, 0.02029238, 0.02464785, 0.02267712]]
        assert_allclose(gauss, ans, rtol=0.001)  # Rough comparison at 0.1 %

    def test_normalize_peak(self):
        """
        Check if normalize works with peak mode.
        """
        custom = CustomKernel([1, 2, 3, 2, 1])
        custom.normalize(mode='peak')
        assert custom.array.max() == 1

    def test_check_kernel_attributes(self):
        """
        Check if kernel attributes are correct.
        """
        box = Box2DKernel(5)

        # Check truncation
        assert box.truncation == 0

        # Check model
        assert isinstance(box.model, Box2D)

        # Check center
        assert box.center == [2, 2]

        # Check normalization
        box.normalize()
        assert_almost_equal(box._kernel_sum, 1., decimal=12)

        # Check separability
        assert box.separable

    @pytest.mark.parametrize(('kernel_type', 'mode'), list(itertools.product(KERNEL_TYPES, MODES)))
    def test_discretize_modes(self, kernel_type, mode):
        """
        Check if the different modes result in kernels that work with convolve.
        Use only small kernel width, to make the test pass quickly.
        """
        if kernel_type == AiryDisk2DKernel and not HAS_SCIPY:
            pytest.skip("Omitting AiryDisk2DKernel, which requires SciPy")
        if not kernel_type == Ring2DKernel:
            kernel = kernel_type(3)
        else:
            kernel = kernel_type(3, 3 * 0.2)

        if kernel.dimension == 1:
            c1 = convolve_fft(delta_pulse_1D, kernel, boundary='fill', normalize_kernel=False)
            c2 = convolve(delta_pulse_1D, kernel, boundary='fill', normalize_kernel=False)
            assert_almost_equal(c1, c2, decimal=12)
        else:
            c1 = convolve_fft(delta_pulse_2D, kernel, boundary='fill', normalize_kernel=False)
            c2 = convolve(delta_pulse_2D, kernel, boundary='fill', normalize_kernel=False)
            assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('width'), WIDTHS_EVEN)
    def test_box_kernels_even_size(self, width):
        """
        Check if BoxKernel work properly with even sizes.
        """
        kernel_1D = Box1DKernel(width)
        assert kernel_1D.shape[0] % 2 != 0
        assert kernel_1D.array.sum() == 1.

        kernel_2D = Box2DKernel(width)
        assert np.all([_ % 2 != 0 for _ in kernel_2D.shape])
        assert kernel_2D.array.sum() == 1.

    def test_kernel_normalization(self):
        """
        Test that repeated normalizations do not change the kernel [#3747].
        """

        kernel = CustomKernel(np.ones(5))
        kernel.normalize()
        data = np.copy(kernel.array)

        kernel.normalize()
        assert_allclose(data, kernel.array)

        kernel.normalize()
        assert_allclose(data, kernel.array)

    def test_kernel_normalization_mode(self):
        """
        Test that an error is raised if mode is invalid.
        """
        with pytest.raises(ValueError):
            kernel = CustomKernel(np.ones(3))
            kernel.normalize(mode='invalid')

    def test_kernel1d_int_size(self):
        """
        Test that an error is raised if ``Kernel1D`` ``x_size`` is not
        an integer.
        """
        with pytest.raises(TypeError):
            Gaussian1DKernel(3, x_size=1.2)

    def test_kernel2d_int_xsize(self):
        """
        Test that an error is raised if ``Kernel2D`` ``x_size`` is not
        an integer.
        """
        with pytest.raises(TypeError):
            Gaussian2DKernel(3, x_size=1.2)

    def test_kernel2d_int_ysize(self):
        """
        Test that an error is raised if ``Kernel2D`` ``y_size`` is not
        an integer.
        """
        with pytest.raises(TypeError):
            Gaussian2DKernel(3, x_size=5, y_size=1.2)

    def test_kernel1d_initialization(self):
        """
        Test that an error is raised if an array or model is not
        specified for ``Kernel1D``.
        """
        with pytest.raises(TypeError):
            Kernel1D()

    def test_kernel2d_initialization(self):
        """
        Test that an error is raised if an array or model is not
        specified for ``Kernel2D``.
        """
        with pytest.raises(TypeError):
            Kernel2D()

    def test_array_keyword_not_allowed(self):
        """
        Regression test for issue #10439
        """
        x = np.ones([10, 10])
        with pytest.raises(TypeError, match=r".* allowed .*"):
            AiryDisk2DKernel(2, array=x)
            Box1DKernel(2, array=x)
            Box2DKernel(2, array=x)
            Gaussian1DKernel(2, array=x)
            Gaussian2DKernel(2, array=x)
            RickerWavelet1DKernel(2, array=x)
            RickerWavelet2DKernel(2, array=x)
            Model1DKernel(Gaussian1D(1, 0, 2), array=x)
            Model2DKernel(Gaussian2D(1, 0, 0, 2, 2), array=x)
            Ring2DKernel(9, 8, array=x)
            Tophat2DKernel(2, array=x)
            Trapezoid1DKernel(2, array=x)
            Trapezoid1DKernel(2, array=x)

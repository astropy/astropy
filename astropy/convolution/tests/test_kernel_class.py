# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from ...tests.helper import pytest
from ..convolve import convolve, convolve_fft
from ..kernels import *

from numpy.testing import assert_almost_equal

import itertools

try:
    from scipy.ndimage import filters
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

widths = [3, 5, 7, 9]
kernel_types = [Gaussian1DKernel, Gaussian2DKernel,
                Box1DKernel, Box2DKernel,
                Trapezoid1DKernel, TrapezoidDisk2DKernel,
                MexicanHat1DKernel, Tophat2DKernel]

# Test data
delta_pulse_1D = np.zeros(81)
delta_pulse_1D[40] = 1

delta_pulse_2D = np.zeros((81, 81))
delta_pulse_2D[40, 40] = 1

random_data_1D = np.random.rand(61)
random_data_2D = np.random.rand(61, 61)


class TestKernels(object):
    """
    Test class for the built-in convolution kernels.
    """

    @pytest.mark.skipif('not HAS_SCIPY')
    @pytest.mark.parametrize(('width'), widths)
    def test_scipy_filter_gaussian(self, width):
        """
        Test GaussianKernel against SciPy ndimage gaussian filter.
        """
        gauss_kernel_1D = Gaussian1DKernel(width)
        gauss_kernel_1D.normalize()
        gauss_kernel_2D = Gaussian2DKernel(width)
        gauss_kernel_2D.normalize()

        astropy_1D = convolve(delta_pulse_1D, gauss_kernel_1D, boundary='fill')
        astropy_2D = convolve(delta_pulse_2D, gauss_kernel_2D, boundary='fill')

        scipy_1D = filters.gaussian_filter(delta_pulse_1D, width)
        scipy_2D = filters.gaussian_filter(delta_pulse_2D, width)

        assert_almost_equal(astropy_1D, scipy_1D, decimal=12)
        assert_almost_equal(astropy_2D, scipy_2D, decimal=12)

    @pytest.mark.skipif('not HAS_SCIPY')
    @pytest.mark.parametrize(('width'), widths)
    def test_scipy_filter_gaussian_laplace(self, width):
        """
        Test MexicanHat kernels against SciPy ndimage gaussian laplace filters.
        """
        mexican_kernel_1D = MexicanHat1DKernel(width)
        mexican_kernel_2D = MexicanHat2DKernel(width)

        astropy_1D = convolve(delta_pulse_1D, mexican_kernel_1D, boundary='fill')
        astropy_2D = convolve(delta_pulse_2D, mexican_kernel_2D, boundary='fill')

        scipy_1D = filters.gaussian_laplace(delta_pulse_1D, width)
        scipy_2D = filters.gaussian_laplace(delta_pulse_2D, width)

        # There is a slight deviation in the normalization. They differ by a
        # factor of ~1.0000284132604045. The reason is not known.
        assert_almost_equal(astropy_1D, scipy_1D, decimal=5)
        assert_almost_equal(astropy_2D, scipy_2D, decimal=5)

    @pytest.mark.parametrize(('kernel_type', 'width'), list(itertools.product(kernel_types, widths)))
    def test_delta_data(self, kernel_type,  width):
        """
        Test smoothing of an image with a single positive pixel
        """
        kernel = kernel_type(width)

        if kernel.dimension == 1:
            c1 = convolve_fft(delta_pulse_1D, kernel, boundary='fill')
            c2 = convolve(delta_pulse_1D, kernel, boundary='fill')
            assert_almost_equal(c1, c2, decimal=12)
        else:
            c1 = convolve_fft(delta_pulse_2D, kernel, boundary='fill')
            c2 = convolve(delta_pulse_2D, kernel, boundary='fill')
            assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('kernel_type', 'width'), list(itertools.product(kernel_types, widths)))
    def test_random_data(self, kernel_type, width):
        """
        Test smoothing of an image made of random noise
        """
        kernel = kernel_type(width)

        if kernel.dimension == 1:
            c1 = convolve_fft(random_data_1D, kernel, boundary='fill')
            c2 = convolve(random_data_1D, kernel, boundary='fill')
            assert_almost_equal(c1, c2, decimal=12)
        else:
            c1 = convolve_fft(random_data_2D, kernel, boundary='fill')
            c2 = convolve(random_data_2D, kernel, boundary='fill')
            assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('width'), widths)
    def test_uniform_smallkernel(self, width):
        """
        Test smoothing of an image with a single positive pixel

        Instead of using kernel class, uses a simple, small kernel
        """
        kernel = np.ones([width, width])

        c2 = convolve_fft(delta_pulse_2D, kernel, boundary='fill')
        c1 = convolve(delta_pulse_2D, kernel, boundary='fill')
        assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('width'), widths)
    def test_smallkernel_vs_Box2DKernel(self, width):
        """
        Test smoothing of an image with a single positive pixel

        Compares a small kernel to something produced by makekernel
        """
        kernel1 = np.ones([width, width]) / width ** 2
        kernel2 = Box2DKernel(width)

        c2 = convolve_fft(delta_pulse_2D, kernel2, boundary='fill')
        c1 = convolve_fft(delta_pulse_2D, kernel1, boundary='fill')

        assert_almost_equal(c1, c2, decimal=12)

    def test_convolve_1D_kernels(self):
        """
        Check if convolving two kernels with eachother works correctly.
        """
        gauss_1 = Gaussian1DKernel(3)
        gauss_2 = Gaussian1DKernel(4)
        test_gauss_3 = Gaussian1DKernel(5)

        gauss_3 = convolve(gauss_1, gauss_2)
        assert np.all(np.abs((gauss_3 - test_gauss_3).array) < 0.01)

    def test_convolve_2D_kernels(self):
        """
        Check if convolving two kernels with eachother works correctly.
        """
        gauss_1 = Gaussian2DKernel(3)
        gauss_2 = Gaussian2DKernel(4)
        test_gauss_3 = Gaussian2DKernel(5)

        gauss_3 = convolve(gauss_1, gauss_2)
        assert np.all(np.abs((gauss_3 - test_gauss_3).array) < 0.01)

    def test_multiply_scalar(self):
        """
        Check if multiplying two kernels with eachother works correctly.
        """
        gauss = Gaussian1DKernel(3)
        assert np.all(np.abs(3 * gauss.array - 3 * gauss) < 0.000001)

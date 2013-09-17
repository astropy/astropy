# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from ....tests.helper import pytest

from ..make_kernel import make_kernel
from ..convolve import convolve, convolve_fft

from numpy.testing import assert_array_almost_equal_nulp, assert_almost_equal

import itertools

shapes = [[8, 8], [15, 15], [16, 16], [31, 31], [32, 32]]
widths = [2, 3, 4, 5]
kerneltype = ['gaussian', 'tophat', 'boxcar']


class Test2DConvolutions(object):

    @pytest.mark.parametrize(('shape', 'width', 'kerneltype'), list(itertools.product(shapes, widths, kerneltype)))
    def test_centered_makekernel(self, shape, width, kerneltype):
        """
        Test smoothing of an image with a single positive pixel
        """

        if width % 2 == 0 and kerneltype == 'boxcar':
            # this is a shifting kernel.  I don't understand how these are treated.
            return

        kernel = make_kernel(shape, width, force_odd=True, kerneltype=kerneltype)

        x = np.zeros(shape)
        xslice = [slice(sh // 2, sh // 2 + 1) for sh in shape]
        x[xslice] = 1.0

        c2 = convolve_fft(x, kernel, boundary='fill')
        c1 = convolve(x, kernel, boundary='fill')

        print shape, width, kerneltype
        assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('shape', 'width', 'kerneltype'), list(itertools.product(shapes, widths, kerneltype)))
    def test_random_makekernel(self, shape, width, kerneltype):
        """
        Test smoothing of an image made of random noise
        """

        if width % 2 == 0 and kerneltype == 'boxcar':
            # this is a shifting kernel.  I don't understand how these are treated.
            return

        kernel = make_kernel(shape, width, force_odd=True, kerneltype=kerneltype)

        x = np.random.randn(*shape)

        c2 = convolve_fft(x, kernel, boundary='fill')
        c1 = convolve(x, kernel, boundary='fill')

        print shape, width, kerneltype
        # not clear why, but these differ by a couple ulps...
        assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('shape', 'width'), list(itertools.product(shapes, widths)))
    def test_uniform_smallkernel(self, shape, width):
        """
        Test smoothing of an image with a single positive pixel

        Instead of using make_kernel, uses a simple, small kernel
        """

        if width % 2 == 0:
            # convolve does not accept odd-shape kernels
            return

        kernel = np.ones([width, width])

        x = np.zeros(shape)
        xslice = [slice(sh // 2, sh // 2 + 1) for sh in shape]
        x[xslice] = 1.0

        c2 = convolve_fft(x, kernel, boundary='fill')
        c1 = convolve(x, kernel, boundary='fill')

        print shape, width
        assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('shape', 'width'), list(itertools.product(shapes, widths)))
    def test_smallkernel_vs_makekernel(self, shape, width):
        """
        Test smoothing of an image with a single positive pixel

        Compares a small kernel to something produced by makekernel
        """

        kernel1 = np.ones([width, width]) / np.float(width) ** 2
        kernel2 = make_kernel(shape, width, kerneltype='boxcar')

        x = np.zeros(shape)
        xslice = [slice(sh // 2, sh // 2 + 1) for sh in shape]
        x[xslice] = 1.0

        c2 = convolve_fft(x, kernel2, boundary='fill')
        c1 = convolve_fft(x, kernel1, boundary='fill')

        print shape, width
        assert_almost_equal(c1, c2, decimal=12)

        if width % 2 == 1:
            kernel2 = make_kernel(shape, width, kerneltype='boxcar', force_odd=True)

            c2 = convolve(x, kernel2, boundary='fill')
            c1 = convolve(x, kernel1, boundary='fill')

            print shape, width

            assert_almost_equal(c1, c2, decimal=12)

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools

import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal

from astropy import units as u
from astropy.convolution.convolve import convolve, convolve_fft
from astropy.convolution.kernels import (Box2DKernel, Gaussian2DKernel,
                                         Moffat2DKernel, Tophat2DKernel)
from astropy.utils.exceptions import AstropyDeprecationWarning

SHAPES_ODD = [[15, 15], [31, 31]]
SHAPES_EVEN = [[8, 8], [16, 16], [32, 32]]  # FIXME: not used ?!
NOSHAPE = [[None, None]]
WIDTHS = [2, 3, 4, 5]

KERNELS = []

for shape in SHAPES_ODD + NOSHAPE:
    for width in WIDTHS:

        KERNELS.append(Gaussian2DKernel(width,
                                        x_size=shape[0],
                                        y_size=shape[1],
                                        mode='oversample',
                                        factor=10))

        KERNELS.append(Box2DKernel(width,
                                   x_size=shape[0],
                                   y_size=shape[1],
                                   mode='oversample',
                                   factor=10))

        KERNELS.append(Tophat2DKernel(width,
                                      x_size=shape[0],
                                      y_size=shape[1],
                                      mode='oversample',
                                      factor=10))
        KERNELS.append(Moffat2DKernel(width, 2,
                                      x_size=shape[0],
                                      y_size=shape[1],
                                      mode='oversample',
                                      factor=10))


class Test2DConvolutions:

    @pytest.mark.parametrize('kernel', KERNELS)
    def test_centered_makekernel(self, kernel):
        """
        Test smoothing of an image with a single positive pixel
        """

        shape = kernel.array.shape

        x = np.zeros(shape)
        xslice = tuple([slice(sh // 2, sh // 2 + 1) for sh in shape])
        x[xslice] = 1.0

        c2 = convolve_fft(x, kernel, boundary='fill')
        c1 = convolve(x, kernel, boundary='fill')

        assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize('kernel', KERNELS)
    def test_random_makekernel(self, kernel):
        """
        Test smoothing of an image made of random noise
        """

        shape = kernel.array.shape

        x = np.random.randn(*shape)

        c2 = convolve_fft(x, kernel, boundary='fill')
        c1 = convolve(x, kernel, boundary='fill')

        # not clear why, but these differ by a couple ulps...
        assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('shape', 'width'), list(itertools.product(SHAPES_ODD, WIDTHS)))
    def test_uniform_smallkernel(self, shape, width):
        """
        Test smoothing of an image with a single positive pixel

        Uses a simple, small kernel
        """

        if width % 2 == 0:
            # convolve does not accept odd-shape kernels
            return

        kernel = np.ones([width, width])

        x = np.zeros(shape)
        xslice = tuple([slice(sh // 2, sh // 2 + 1) for sh in shape])
        x[xslice] = 1.0

        c2 = convolve_fft(x, kernel, boundary='fill')
        c1 = convolve(x, kernel, boundary='fill')

        assert_almost_equal(c1, c2, decimal=12)

    @pytest.mark.parametrize(('shape', 'width'), list(itertools.product(SHAPES_ODD, [1, 3, 5])))
    def test_smallkernel_Box2DKernel(self, shape, width):
        """
        Test smoothing of an image with a single positive pixel

        Compares a small uniform kernel to the Box2DKernel
        """

        kernel1 = np.ones([width, width]) / float(width) ** 2
        kernel2 = Box2DKernel(width, mode='oversample', factor=10)

        x = np.zeros(shape)
        xslice = tuple([slice(sh // 2, sh // 2 + 1) for sh in shape])
        x[xslice] = 1.0

        c2 = convolve_fft(x, kernel2, boundary='fill')
        c1 = convolve_fft(x, kernel1, boundary='fill')

        assert_almost_equal(c1, c2, decimal=12)

        c2 = convolve(x, kernel2, boundary='fill')
        c1 = convolve(x, kernel1, boundary='fill')

        assert_almost_equal(c1, c2, decimal=12)


def test_gaussian_2d_kernel_quantity():
    # Make sure that the angle can be a quantity
    kernel1 = Gaussian2DKernel(x_stddev=2, y_stddev=4, theta=45 * u.deg)
    kernel2 = Gaussian2DKernel(x_stddev=2, y_stddev=4, theta=np.pi / 4)
    assert_allclose(kernel1.array, kernel2.array)


def test_deprecated_hat():

    # 'MexicanHat' was deprecated as a name for the kernels which are now
    # 'RickerWavelet'. This test ensures that the kernels are correctly
    # deprecated, and can be imported from the top-level package.

    from astropy.convolution import MexicanHat1DKernel, MexicanHat2DKernel

    with pytest.warns(AstropyDeprecationWarning):
        MexicanHat1DKernel(2)

    with pytest.warns(AstropyDeprecationWarning):
        MexicanHat2DKernel(2)

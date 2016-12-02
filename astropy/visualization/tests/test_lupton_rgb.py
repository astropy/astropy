# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for RGB Images
"""

from __future__ import division, print_function

import sys
import os
import tempfile

import numpy as np
from numpy.testing import assert_equal

from ...convolution import convolve, Gaussian2DKernel
from ...tests.helper import pytest
from .. import lupton_rgb

import matplotlib
matplotlib.use('TkAgg')

try:
    import scipy.misc
    scipy.misc.imresize  # just checking if it exists.
    HAVE_SCIPY_MISC = True
except (ImportError, AttributeError):
    HAVE_SCIPY_MISC = False

# set to True to have the finished RGBs displayed on your screen.
display = False


def my_name():
    """Return the name of the current method."""
    return sys._getframe().f_code.co_name


def saturate(image, satValue):
    """
    Return image with all points above satValue set to NaN.

    Simulates saturation on an image, so we can test 'replaceSaturatedPixels'
    """
    # TBD: FROMAFW
    # image = afwImage.makeMaskedImage(image)
    # afwDetect.FootprintSet(image, afwDetect.Threshold(satValue), "SAT")
    # arr = image.getImage().getArray()
    # arr[np.where(arr >= satValue)] = np.nan

    result = image.copy()
    saturated = image > satValue
    result[saturated] = np.nan
    return result


def random_array(dtype, N=100):
    return np.array(np.random.random(10)*100, dtype=dtype)


def test_compute_intensity_1_float():
    imageR = random_array(np.float64)
    intensity = lupton_rgb.compute_intensity(imageR)
    assert imageR.dtype == intensity.dtype
    assert_equal(imageR, intensity)


def test_compute_intensity_1_uint():
    imageR = random_array(np.uint8)
    intensity = lupton_rgb.compute_intensity(imageR)
    assert imageR.dtype == intensity.dtype
    assert_equal(imageR, intensity)


def test_compute_intensity_3_float():
    imageR = random_array(np.float64)
    imageG = random_array(np.float64)
    imageB = random_array(np.float64)
    intensity = lupton_rgb.compute_intensity(imageR, imageG, imageB)
    assert imageR.dtype == intensity.dtype
    assert_equal(intensity, (imageR+imageG+imageB)/3.0)


def test_compute_intensity_3_uint():
    imageR = random_array(np.uint8)
    imageG = random_array(np.uint8)
    imageB = random_array(np.uint8)
    intensity = lupton_rgb.compute_intensity(imageR, imageG, imageB)
    assert imageR.dtype == intensity.dtype
    assert_equal(intensity, (imageR+imageG+imageB)//3)


class TestLuptonRgb(object):
    """A test case for Rgb"""

        # TBD: Old version. converted to the above.
        # TBD: remove when it's clear the above is what was meant/works.

        # self.images.append(afwImage.ImageF(afwGeom.ExtentI(width, height)))
        # self.images.append(afwImage.ImageF(afwGeom.ExtentI(width, height)))
        # self.images.append(afwImage.ImageF(afwGeom.ExtentI(width, height)))

        # for (x, y, A, g_r, r_i) in [(15, 15, 1000,  1.0,  2.0),
        #                             (50, 45, 5500, -1.0, -0.5),
        #                             (30, 30,  600,  1.0,  2.5),
        #                             (45, 15, 20000,  1.0,  1.0),
        #                             ]:
        #     for i in range(len(self.images)):
        #         if i == B:
        #             amp = A
        #         elif i == G:
        #             amp = A*math.pow(10, 0.4*g_r)
        #         elif i == R:
        #             amp = A*math.pow(10, 0.4*r_i)

        #         self.images[i].set(x, y, amp)

        # psf = afwMath.AnalyticKernel(15, 15, afwMath.GaussianFunction2D(2.5, 1.5, 0.5))
        # convolvedImage = type(self.images[0])(self.images[0].getDimensions())
        # randomImage = type(self.images[0])(self.images[0].getDimensions())
        # rand = afwMath.Random("MT19937", 666)
        # for i in range(len(self.images)):
        #     afwMath.convolve(convolvedImage, self.images[i], psf, True, True)
        #     afwMath.randomGaussianImage(randomImage, rand)
        #     randomImage *= 2
        #     convolvedImage += randomImage
        #     self.images[i][:] = convolvedImage
        # del convolvedImage
        # del randomImage
    np.random.seed(1000)  # so we always get the same images.

    min_, range_, Q = 0, 5, 20  # asinh

    width, height = 85, 75
    images = []

    shape = (width, height)
    imageR = np.zeros(shape)
    imageG = np.zeros(shape)
    imageB = np.zeros(shape)

    # pixel locations, values and colors
    points = [[15, 15], [50, 45], [30, 30], [45, 15]]
    values = [1000, 5500, 600, 20000]
    g_r = [1.0, -1.0, 1.0, 1.0]
    r_i = [2.0, -0.5, 2.5, 1.0]

    # Put pixels in the images.
    for p, v, gr, ri in zip(points, values, g_r, r_i):
        imageR[p[0], p[1]] = v*pow(10, 0.4*ri)
        imageG[p[0], p[1]] = v*pow(10, 0.4*gr)
        imageB[p[0], p[1]] = v

    # convolve the image with a reasonable PSF, and add Gaussian background noise
    def convolve_with_noise(image, psf):
        convolvedImage = convolve(image, psf, boundary='extend', normalize_kernel=True)
        randomImage = np.random.normal(0, 2, image.shape)
        return randomImage + convolvedImage

    psf = Gaussian2DKernel(2.5)
    imageR = convolve_with_noise(imageR, psf)
    imageG = convolve_with_noise(imageG, psf)
    imageB = convolve_with_noise(imageB, psf)

    def tearDown(self):
        for im in self.images:
            del im
        del self.images

    def testStarsAsinh(self):
        """Test creating an RGB image using an asinh stretch"""
        asinhMap = lupton_rgb.AsinhMapping(self.min_, self.range_, self.Q)
        rgbImage = asinhMap.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def testStarsAsinhZscale(self):
        """Test creating an RGB image using an asinh stretch estimated using zscale"""

        map = lupton_rgb.AsinhZScaleMapping(self.imageR, self.imageG, self.imageB)
        rgbImage = map.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def testStarsAsinhZscaleIntensity(self):
        """Test creating an RGB image using an asinh stretch estimated using zscale on the intensity"""

        map = lupton_rgb.AsinhZScaleMapping(self.imageR, self.imageG, self.imageB)
        rgbImage = map.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def testStarsAsinhZscaleIntensityPedestal(self):
        """Test creating an RGB image using an asinh stretch estimated using zscale on the intensity
        where the images each have a pedestal added"""

        pedestal = [100, 400, -400]
        self.imageR += pedestal[0]
        self.imageG += pedestal[1]
        self.imageB += pedestal[2]

        map = lupton_rgb.AsinhZScaleMapping(self.imageR, self.imageG, self.imageB, pedestal=pedestal)
        rgbImage = map.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def testStarsAsinhZscaleIntensityBW(self):
        """Test creating a black-and-white image using an asinh stretch estimated
        using zscale on the intensity"""

        rgbImage = lupton_rgb.AsinhZScaleMapping(self.imageR).makeRgbImage()

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def testMakeRGB(self):
        """Test the function that does it all"""
        satValue = 1000.0
        with tempfile.NamedTemporaryFile(suffix=".png") as temp:
            red = saturate(self.imageR, satValue)
            green = saturate(self.imageG, satValue)
            blue = saturate(self.imageB, satValue)
            lupton_rgb.makeRGB(red, green, blue, self.min_, self.range_, self.Q, fileName=temp)
            assert os.path.exists(temp.name)

    def testMakeRGB_saturated_fix(self):
        satValue = 1000.0
        # TODO: Cannot test with these options yet, as that part of the code is not implemented.
        with pytest.raises(NotImplementedError):
            red = saturate(self.imageR, satValue)
            green = saturate(self.imageG, satValue)
            blue = saturate(self.imageB, satValue)
            lupton_rgb.makeRGB(red, green, blue, self.min_, self.range_, self.Q,
                               saturatedBorderWidth=1, saturatedPixelValue=2000)

    def testLinear(self):
        """Test using a specified linear stretch"""

        rgbImage = lupton_rgb.LinearMapping(-8.45, 13.44).makeRgbImage(self.imageR)

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def testLinearMinMax(self):
        """Test using a min/max linear stretch

        N.b. also checks that an image passed to the ctor is used as the default in makeRgbImage()
        """

        rgbImage = lupton_rgb.LinearMapping(image=self.imageR).makeRgbImage()

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def testZScale(self):
        """Test using a zscale stretch"""

        rgbImage = lupton_rgb.ZScaleMapping(self.imageR).makeRgbImage()

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def testWriteStars(self):
        """Test writing RGB files to disk"""
        asinhMap = lupton_rgb.AsinhMapping(self.min_, self.range_, self.Q)
        rgbImage = asinhMap.makeRgbImage(self.imageR, self.imageG, self.imageB)
        with tempfile.NamedTemporaryFile(suffix=".png") as temp:
            lupton_rgb.writeRGB(temp, rgbImage)
            assert os.path.exists(temp.name)

    def testSaturated(self):
        """Test interpolating saturated pixels"""
        pytest.skip('replaceSaturatedPixels is not implemented in astropy yet')

        satValue = 1000.0
        self.imageR = saturate(self.imageR, satValue)
        self.imageG = saturate(self.imageG, satValue)
        self.imageB = saturate(self.imageB, satValue)

        lupton_rgb.replaceSaturatedPixels(self.imageR, self.imageG, self.imageB, 1, 2000)
        # Check that we replaced those NaNs with some reasonable value
        assert np.isfinite(self.imageR.getImage().getArray()).all()
        assert np.isfinite(self.imageG.getImage().getArray()).all()
        assert np.isfinite(self.imageB.getImage().getArray()).all()

        # Prepare for generating an output file
        # TBD: FROMAFW
        self.imagesR = self.imagesR.getImage()
        self.imagesR = self.imagesG.getImage()
        self.imagesR = self.imagesB.getImage()

        asinhMap = lupton_rgb.AsinhMapping(self.min_, self.range_, self.Q)
        rgbImage = asinhMap.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testStarsResizeToSize(self):
        """Test creating an RGB image of a specified size"""

        xSize = self.imageR.shape[0]/2
        ySize = self.imageR.shape[1]/2
        asinhZ = lupton_rgb.AsinhZScaleMapping(self.imageR)
        rgbImage = asinhZ.makeRgbImage(self.imageR, self.imageG, self.imageB, xSize=xSize, ySize=ySize)

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testStarsResizeToSize_uint(self):
        """Test creating an RGB image of a specified size"""

        xSize = self.imageR.shape[0]/2
        ySize = self.imageR.shape[1]/2
        imageR = np.array(self.imageR, dtype=np.uint16)
        imageG = np.array(self.imageG, dtype=np.uint16)
        imageB = np.array(self.imageB, dtype=np.uint16)
        asinhZ = lupton_rgb.AsinhZScaleMapping(imageR)
        rgbImage = asinhZ.makeRgbImage(imageR, imageG, imageB, xSize=xSize, ySize=ySize)

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    def _testStarsResizeSpecifications(self, xSize=None, ySize=None, frac=None):
        """Test creating an RGB image changing the output """

        map = lupton_rgb.AsinhZScaleMapping(self.imageR)
        rgbImage = map.makeRgbImage(self.imageR, self.imageG, self.imageB,
                                    xSize=xSize, ySize=ySize, rescaleFactor=frac)

        # TODO: I'm not positive that I got the width/height values correct:
        # afw and numpy have different conventions!
        h, w, _ = rgbImage.shape
        if xSize is not None:
            assert int(xSize) == w
        if ySize is not None:
            assert int(ySize) == h
        if frac is not None:
            assert int(frac*self.imageR.shape[1]) == w
            assert int(frac*self.imageR.shape[0]) == h

        if display:
            lupton_rgb.displayRGB(rgbImage, title=my_name())

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testStarsResizeSpecificaions_xSize_ySize(self):
        self._testStarsResizeSpecifications(xSize=self.imageR.shape[0]/2, ySize=self.imageR.shape[1]/2)

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testStarsResizeSpecifications_twice_xSize(self):
        self._testStarsResizeSpecifications(xSize=2*self.imageR.shape[0])

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testStarsResizeSpecifications_half_xSize(self):
        self._testStarsResizeSpecifications(xSize=self.imageR.shape[0]/2)

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testStarsResizeSpecifications_half_ySize(self):
        self._testStarsResizeSpecifications(ySize=self.imageR.shape[0]/2)

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testStarsResizeSpecifications_frac_half(self):
        self._testStarsResizeSpecifications(frac=0.5)

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testStarsResizeSpecifications_frac_twice(self):
        self._testStarsResizeSpecifications(frac=2)

    @pytest.mark.skipif('not HAVE_SCIPY_MISC', reason="Resizing images requires scipy.misc")
    def testMakeRGBResize(self):
        """Test the function that does it all, including rescaling"""
        lupton_rgb.makeRGB(self.imageR, self.imageG, self.imageB, xSize=40, ySize=60)

        with tempfile.NamedTemporaryFile(suffix=".png") as temp:
            lupton_rgb.makeRGB(self.imageR, self.imageG, self.imageB, fileName=temp, rescaleFactor=0.5)
            assert os.path.exists(temp.name)

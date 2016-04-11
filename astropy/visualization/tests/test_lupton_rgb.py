# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for RGB Images
"""

from __future__ import division, print_function

import os
import unittest
import tempfile

import numpy as np
from numpy.testing import assert_equal

from ...convolution import convolve, Gaussian2DKernel
from .. import lupton_rgb

ver1, ver2, ver3 = 1, 3, 1
NO_MATPLOTLIB_STRING = "Requires matplotlib >= %d.%d.%d" % (ver1, ver2, ver3)
try:
    import matplotlib
    versionInfo = tuple(int(s.strip("rc")) for s in matplotlib.__version__.split("."))
    HAVE_MATPLOTLIB = versionInfo >= (ver1, ver2, ver3)
except ImportError:
    HAVE_MATPLOTLIB = False

try:
    import scipy.misc
    scipy.misc.imresize
    HAVE_SCIPY_MISC = True
except (ImportError, AttributeError):
    HAVE_SCIPY_MISC = False

try:
    type(display)
except NameError:
    display = False


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


class TestLuptonRgb(unittest.TestCase):
    """A test case for Rgb"""
    def setUp(self):
        self.min, self.range, self.Q = 0, 5, 20  # asinh

        width, height = 85, 75
        self.images = []

        shape = (width, height)
        self.imageR = np.zeros(shape)
        self.imageG = np.zeros(shape)
        self.imageB = np.zeros(shape)

        # pixel locations, values and colors
        points = [[15, 15], [50, 45], [30, 30], [45, 15]]
        values = [1000, 5500, 600, 20000]
        g_r = [1.0, -1.0, 1.0, 1.0]
        r_i = [2.0, -0.5, 2.5, 1.0]

        # Put pixels in the images.
        for p, v, gr, ri in zip(points, values, g_r, r_i):
            self.imageR[p[0], p[1]] = v*pow(10, 0.4*ri)
            self.imageG[p[0], p[1]] = v*pow(10, 0.4*gr)
            self.imageB[p[0], p[1]] = v

        # convolve the image with a reasonable PSF, and add Gaussian background noise
        def convolve_with_noise(image, psf):
            convolvedImage = convolve(image, psf, boundary='extend', normalize_kernel=True)
            randomImage = np.random.normal(0, 2, image.shape)
            return randomImage + convolvedImage

        psf = Gaussian2DKernel(2.5)
        self.imageR = convolve_with_noise(self.imageR, psf)
        self.imageG = convolve_with_noise(self.imageG, psf)
        self.imageB = convolve_with_noise(self.imageB, psf)

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

    def tearDown(self):
        for im in self.images:
            del im
        del self.images

    def testStarsAsinh(self):
        """Test creating an RGB image using an asinh stretch"""
        asinhMap = lupton_rgb.AsinhMapping(self.min, self.range, self.Q)
        rgbImage = asinhMap.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage)

    def testStarsAsinhZscale(self):
        """Test creating an RGB image using an asinh stretch estimated using zscale"""

        map = lupton_rgb.AsinhZScaleMapping(self.imageR, self.imageG, self.imageB)
        rgbImage = map.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage)

    def testStarsAsinhZscaleIntensity(self):
        """Test creating an RGB image using an asinh stretch estimated using zscale on the intensity"""

        map = lupton_rgb.AsinhZScaleMapping(self.imageR, self.imageG, self.imageB)
        rgbImage = map.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage)

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
            lupton_rgb.displayRGB(rgbImage)

    def testStarsAsinhZscaleIntensityBW(self):
        """Test creating a black-and-white image using an asinh stretch estimated
        using zscale on the intensity"""

        rgbImage = lupton_rgb.AsinhZScaleMapping(self.imageR).makeRgbImage()

        if display:
            lupton_rgb.displayRGB(rgbImage)

    @unittest.skipUnless(HAVE_MATPLOTLIB, NO_MATPLOTLIB_STRING)
    def testMakeRGB(self):
        """Test the function that does it all"""
        satValue = 1000.0
        with tempfile.NamedTemporaryFile(suffix=".png") as temp:
            red = saturate(self.imageR, satValue)
            green = saturate(self.imageG, satValue)
            blue = saturate(self.imageB, satValue)
            lupton_rgb.makeRGB(red, green, blue, self.min, self.range, self.Q, fileName=temp,
                               saturatedBorderWidth=1, saturatedPixelValue=2000)
            self.assertTrue(os.path.exists(temp.name))

    def testLinear(self):
        """Test using a specified linear stretch"""

        rgbImage = lupton_rgb.LinearMapping(-8.45, 13.44).makeRgbImage(self.imageR)

        if display:
            lupton_rgb.displayRGB(rgbImage)

    def testLinearMinMax(self):
        """Test using a min/max linear stretch

        N.b. also checks that an image passed to the ctor is used as the default in makeRgbImage()
        """

        rgbImage = lupton_rgb.LinearMapping(image=self.imageR).makeRgbImage()

        if display:
            lupton_rgb.displayRGB(rgbImage)

    def testZScale(self):
        """Test using a zscale stretch"""

        rgbImage = lupton_rgb.ZScaleMapping(self.imageR).makeRgbImage()

        if display:
            plt = lupton_rgb.displayRGB(rgbImage, False)
            plt.title("zscale")
            plt.show()

    @unittest.skipUnless(HAVE_MATPLOTLIB, NO_MATPLOTLIB_STRING)
    def testWriteStars(self):
        """Test writing RGB files to disk"""
        asinhMap = lupton_rgb.AsinhMapping(self.min, self.range, self.Q)
        rgbImage = asinhMap.makeRgbImage(self.imageR, self.imageG, self.imageB)
        with tempfile.NamedTemporaryFile(suffix=".png") as temp:
            lupton_rgb.writeRGB(temp, rgbImage)
            self.assertTrue(os.path.exists(temp.name))

    def testSaturated(self):
        """Test interpolating saturated pixels"""

        satValue = 1000.0
        self.imageR = saturate(self.imageR, satValue)
        self.imageG = saturate(self.imageG, satValue)
        self.imageB = saturate(self.imageB, satValue)

        lupton_rgb.replaceSaturatedPixels(self.imageR, self.imageG, self.imageB, 1, 2000)
        # Check that we replaced those NaNs with some reasonable value
        self.assertTrue(np.isfinite(self.imageR.getImage().getArray()).all())
        self.assertTrue(np.isfinite(self.imageG.getImage().getArray()).all())
        self.assertTrue(np.isfinite(self.imageB.getImage().getArray()).all())

        # Prepare for generating an output file
        # TBD: FROMAFW
        self.imagesR = self.imagesR.getImage()
        self.imagesR = self.imagesG.getImage()
        self.imagesR = self.imagesB.getImage()

        asinhMap = lupton_rgb.AsinhMapping(self.min, self.range, self.Q)
        rgbImage = asinhMap.makeRgbImage(self.imageR, self.imageG, self.imageB)

        if display:
            lupton_rgb.displayRGB(rgbImage)

    @unittest.skipUnless(HAVE_SCIPY_MISC, "Resizing images requires scipy.misc")
    def testStarsResizeToSize(self):
        """Test creating an RGB image of a specified size"""

        xSize = self.imageR.shape[0]/2
        ySize = self.imageR.shape[1]/2
        asinhZ = lupton_rgb.AsinhZScaleMapping(self.imageR)
        rgbImage = asinhZ.makeRgbImage(self.imageR, self.imageG, self.imageB, xSize=xSize, ySize=ySize)

        if display:
            lupton_rgb.displayRGB(rgbImage)

    @unittest.skipUnless(HAVE_SCIPY_MISC, "Resizing images requires scipy.misc")
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
            lupton_rgb.displayRGB(rgbImage)

    @unittest.skipUnless(HAVE_SCIPY_MISC, "Resizing images requires scipy.misc")
    def testStarsResizeSpecifications(self):
        """Test creating an RGB image changing the output """

        map = lupton_rgb.AsinhZScaleMapping(self.imageR)

        for xSize, ySize, frac in [(self.imageR.shape[0]/2, self.imageR.shape[1]/2, None),
                                   (2*self.imageR.shape[0], None,                         None),
                                   (self.imageR.shape[0]/2, None,                         None),
                                   (None,                        self.imageR.shape[1]/2, None),
                                   (None,                        None,                         0.5),
                                   (None,                        None,                         2),
                                   ]:
            rgbImage = map.makeRgbImage(self.imageR, self.imageG, self.imageB,
                                        xSize=xSize, ySize=ySize, rescaleFactor=frac)

            h, w = rgbImage.shape[0:2]
            self.assertTrue(xSize is None or xSize == w)
            self.assertTrue(ySize is None or ySize == h)
            self.assertTrue(frac is None or w == int(frac*self.imageR.shape[0]),
                            "%g == %g" % (w, int((frac if frac else 1)*self.imageR.shape[0])))

            if display:
                lupton_rgb.displayRGB(rgbImage)

    @unittest.skipUnless(HAVE_SCIPY_MISC, "Resizing images requires scipy.misc")
    @unittest.skipUnless(HAVE_MATPLOTLIB, NO_MATPLOTLIB_STRING)
    def testMakeRGBResize(self):
        """Test the function that does it all, including rescaling"""
        lupton_rgb.makeRGB(self.imageR, self.imageG, self.imageB, xSize=40, ySize=60)

        with tempfile.NamedTemporaryFile(suffix=".png") as temp:
            lupton_rgb.makeRGB(self.imageR, self.imageG, self.imageB, fileName=temp, rescaleFactor=0.5)
            self.assertTrue(os.path.exists(temp.name))

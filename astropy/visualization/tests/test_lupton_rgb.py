# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for RGB Images
"""

import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB
from astropy.visualization import lupton_rgb
from astropy.visualization.interval import ManualInterval
from astropy.visualization.stretch import LinearStretch

# Set display=True to get matplotlib imshow windows to help with debugging.
display = False


def display_rgb(rgb, title=None):
    """Display an rgb image using matplotlib (useful for debugging)"""
    import matplotlib.pyplot as plt

    plt.imshow(rgb, interpolation="nearest", origin="lower")
    if title:
        plt.title(title)
    plt.show()
    return plt


def saturate(image, satValue):
    """
    Return image with all points above satValue set to NaN.

    Simulates saturation on an image, so we can test 'replace_saturated_pixels'
    """
    result = image.copy()
    saturated = image > satValue
    result[saturated] = np.nan
    return result


def random_array(dtype, N=100):
    return np.array(np.random.random(10) * 100, dtype=dtype)


def test_compute_intensity_1_float():
    image_r = random_array(np.float64)
    intensity = lupton_rgb.compute_intensity(image_r)
    assert image_r.dtype == intensity.dtype
    assert_equal(image_r, intensity)


def test_compute_intensity_1_uint():
    image_r = random_array(np.uint8)
    intensity = lupton_rgb.compute_intensity(image_r)
    assert image_r.dtype == intensity.dtype
    assert_equal(image_r, intensity)


def test_compute_intensity_3_float():
    image_r = random_array(np.float64)
    image_g = random_array(np.float64)
    image_b = random_array(np.float64)
    intensity = lupton_rgb.compute_intensity(image_r, image_g, image_b)
    assert image_r.dtype == intensity.dtype
    assert_equal(intensity, (image_r + image_g + image_b) / 3.0)


def test_compute_intensity_3_uint():
    image_r = random_array(np.uint8)
    image_g = random_array(np.uint8)
    image_b = random_array(np.uint8)
    intensity = lupton_rgb.compute_intensity(image_r, image_g, image_b)
    assert image_r.dtype == intensity.dtype
    assert_equal(intensity, (image_r + image_g + image_b) // 3)


class TestLuptonRgb:
    """A test case for Rgb"""

    def setup_method(self, method):
        np.random.seed(1000)  # so we always get the same images.

        self.min_, self.stretch_, self.Q = 0, 5, 20  # asinh

        width, height = 85, 75
        self.width = width
        self.height = height

        shape = (width, height)
        image_r = np.zeros(shape)
        image_g = np.zeros(shape)
        image_b = np.zeros(shape)

        # pixel locations, values and colors
        points = [[15, 15], [50, 45], [30, 30], [45, 15]]
        values = [1000, 5500, 600, 20000]
        g_r = [1.0, -1.0, 1.0, 1.0]
        r_i = [2.0, -0.5, 2.5, 1.0]

        # Put pixels in the images.
        for p, v, gr, ri in zip(points, values, g_r, r_i):
            image_r[p[0], p[1]] = v * pow(10, 0.4 * ri)
            image_g[p[0], p[1]] = v * pow(10, 0.4 * gr)
            image_b[p[0], p[1]] = v

        # convolve the image with a reasonable PSF,
        # and add Gaussian background noise
        def convolve_with_noise(image, psf):
            convolvedImage = convolve(
                image, psf, boundary="extend", normalize_kernel=True
            )
            randomImage = np.random.normal(0, 2, image.shape)
            return randomImage + convolvedImage

        psf = Gaussian2DKernel(2.5)
        self.image_r = convolve_with_noise(image_r, psf)
        self.image_g = convolve_with_noise(image_g, psf)
        self.image_b = convolve_with_noise(image_b, psf)

    def test_Asinh(self):
        """Test creating an RGB image using an asinh stretch"""

        asinh_map = lupton_rgb.RGBImageMappingLupton(
            interval=ManualInterval(vmin=self.min_, vmax=None),
            stretch=lupton_rgb.LuptonAsinhStretch(self.stretch_, self.Q),
        )
        rgb_image = asinh_map.make_rgb_image(self.image_r, self.image_g, self.image_b)

        if display:
            display_rgb(rgb_image, title=sys._getframe().f_code.co_name)

    def test_Asinh_incorrect_stretch_asserts(self):
        with pytest.raises(ValueError, match=r"Stretch must be non-negative"):
            _ = lupton_rgb.LuptonAsinhStretch(-1.0, self.Q)

    def test_Asinh_incorrect_Q_asserts(self):
        with pytest.raises(ValueError, match=r"Q must be non-negative"):
            _ = lupton_rgb.LuptonAsinhStretch(self.stretch_, -1.0)

    def test_Asinh_Q_machine_floor(self):
        asinh_map = lupton_rgb.LuptonAsinhStretch(self.stretch_, 1.0e-24)
        assert_allclose(asinh_map.Q, 0.1)

    def test_Asinh_Q_ceil(self):
        asinh_map = lupton_rgb.LuptonAsinhStretch(self.stretch_, 1e11)
        assert_allclose(asinh_map.Q, 1e10)

    def test_AsinhZscale(self):
        """
        Test creating an RGB image using an asinh stretch estimated
        using zscale
        """
        map_ = lupton_rgb.RGBImageMappingLupton(
            interval=ManualInterval(vmin=self.min_, vmax=None),
            stretch=lupton_rgb.LuptonAsinhZscaleStretch(
                [self.image_r, self.image_g, self.image_b], self.Q
            ),
        )
        rgb_image = map_.make_rgb_image(self.image_r, self.image_g, self.image_b)

        if display:
            display_rgb(rgb_image, title=sys._getframe().f_code.co_name)

    def test_AsinhZscaleIntensity(self):
        """
        Test creating an RGB image using an asinh stretch estimated
        using zscale on the intensity
        """
        map_ = lupton_rgb.RGBImageMappingLupton(
            interval=ManualInterval(vmin=self.min_, vmax=None),
            stretch=lupton_rgb.LuptonAsinhZscaleStretch(
                lupton_rgb.compute_intensity(self.image_r, self.image_g, self.image_b),
                self.Q,
            ),
        )
        rgb_image = map_.make_rgb_image(self.image_r, self.image_g, self.image_b)

        if display:
            display_rgb(rgb_image, title=sys._getframe().f_code.co_name)

    def test_AsinhZscaleIntensityBW(self):
        """Test creating a black-and-white image using an asinh stretch
        estimated using zscale on the intensity"""
        map_ = lupton_rgb.RGBImageMappingLupton(
            interval=ManualInterval(vmin=self.min_, vmax=None),
            stretch=lupton_rgb.LuptonAsinhZscaleStretch(
                self.image_r,
                self.Q,
            ),
        )
        rgb_image = map_.make_rgb_image(self.image_r, self.image_r, self.image_r)
        if display:
            display_rgb(rgb_image, title=sys._getframe().f_code.co_name)

    def test_AsinhZscale_pedestal_array(self):
        """
        Test creating an RGB image using an asinh stretch estimated
        using zscale
        """
        map_ = lupton_rgb.RGBImageMappingLupton(
            interval=ManualInterval(vmin=self.min_, vmax=None),
            stretch=lupton_rgb.LuptonAsinhZscaleStretch(
                [self.image_r, self.image_g, self.image_b],
                self.Q,
                pedestal=[1.0, 1.0, 2.0],
            ),
        )
        rgb_image = map_.make_rgb_image(self.image_r, self.image_g, self.image_b)

        if display:
            display_rgb(rgb_image, title=sys._getframe().f_code.co_name)

    def test_AsinhZscale_pedestal_float(self):
        """
        Test creating an RGB image using an asinh stretch estimated
        using zscale
        """
        map_ = lupton_rgb.RGBImageMappingLupton(
            interval=ManualInterval(vmin=self.min_, vmax=None),
            stretch=lupton_rgb.LuptonAsinhZscaleStretch(
                [self.image_r, self.image_g, self.image_b],
                self.Q,
                pedestal=1.0,
            ),
        )
        rgb_image = map_.make_rgb_image(self.image_r, self.image_g, self.image_b)

        if display:
            display_rgb(rgb_image, title=sys._getframe().f_code.co_name)

    def test_AsinhZscale_pedestal_incorrect_assert(self):
        """
        Test creating an RGB image using an asinh stretch estimated
        using zscale
        """
        with pytest.raises(ValueError, match=r"pedestal must be 1 or 3 values"):
            _ = lupton_rgb.LuptonAsinhZscaleStretch(
                [self.image_r, self.image_g, self.image_b],
                self.Q,
                pedestal=[1.0, 2.0],
            )

    def test_AsinhZscale_incorrect_input_asserts(self):
        with pytest.raises(ValueError, match=r"Input 'image' must be a single"):
            _ = lupton_rgb.LuptonAsinhZscaleStretch(
                [self.image_r, self.image_g], self.Q
            )

    def test_AsinhZscale_incorrect_input_nonimage_asserts(self):
        with pytest.raises(ValueError, match=r"Input 'image' must be a single"):
            _ = lupton_rgb.LuptonAsinhZscaleStretch([1], self.Q)

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason="requires matplotlib")
    def test_make_rgb(self, tmp_path):
        """Test the function that does it all"""
        temp = tmp_path.with_suffix(".png")
        lupton_rgb.make_lupton_rgb(
            self.image_r,
            self.image_g,
            self.image_b,
            stretch=self.stretch_,
            Q=self.Q,
            minimum=self.min_,
            filename=temp,
        )
        assert temp.exists()

    def test_make_rgb_incorrect_min_input(self):
        with pytest.raises(ValueError, match=r"3 values for minimum."):
            lupton_rgb.make_lupton_rgb(
                self.image_r,
                self.image_g,
                self.image_b,
                stretch=self.stretch_,
                Q=self.Q,
                minimum=[self.min_, self.min_],
            )

    def test_make_rgb_saturated_fix(self, tmp_path):
        pytest.skip("saturation correction is not implemented")
        satValue = 1000.0
        # TODO: Cannot test with these options yet, as that part of the code
        # is not implemented.
        temp = tmp_path.with_suffix(".png")
        red = saturate(self.image_r, satValue)
        green = saturate(self.image_g, satValue)
        blue = saturate(self.image_b, satValue)
        lupton_rgb.make_lupton_rgb(
            red,
            green,
            blue,
            minimum=self.min_,
            stretch=self.stretch_,
            Q=self.Q,
            saturated_border_width=1,
            saturated_pixel_value=2000,
            filename=temp,
        )

    def test_linear(self):
        """Test using a specified linear stretch"""

        map_ = lupton_rgb.RGBImageMappingLupton(
            interval=ManualInterval(vmin=self.min_, vmax=None),
            stretch=LinearStretch(-8.45, 13.44),
        )
        rgb_image = map_.make_rgb_image(self.image_r, self.image_g, self.image_b)
        if display:
            display_rgb(rgb_image, title=sys._getframe().f_code.co_name)

    def test_saturated(self):
        """Test interpolationolating saturated pixels"""
        pytest.skip("replaceSaturatedPixels is not implemented in astropy yet")

        satValue = 1000.0
        self.image_r = saturate(self.image_r, satValue)
        self.image_g = saturate(self.image_g, satValue)
        self.image_b = saturate(self.image_b, satValue)

        lupton_rgb.replaceSaturatedPixels(
            self.image_r, self.image_g, self.image_b, 1, 2000
        )
        # Check that we replaced those NaNs with some reasonable value
        assert np.isfinite(self.image_r.getImage().getArray()).all()
        assert np.isfinite(self.image_g.getImage().getArray()).all()
        assert np.isfinite(self.image_b.getImage().getArray()).all()

        # Prepare for generating an output file
        self.imagesR = self.imagesR.getImage()
        self.imagesR = self.imagesG.getImage()
        self.imagesR = self.imagesB.getImage()

        asinhMap = lupton_rgb.AsinhMapping(self.min_, self.stretch_, self.Q)
        rgb_image = asinhMap.make_rgb_image(self.image_r, self.image_g, self.image_b)
        if display:
            display_rgb(rgb_image, title=sys._getframe().f_code.co_name)

    def test_different_shapes_asserts(self):
        with pytest.raises(ValueError, match=r"shapes must match"):
            # just swap the dimensions to get a differently-shaped 'r'
            image_r = self.image_r.reshape(self.height, self.width)
            lupton_rgb.make_lupton_rgb(image_r, self.image_g, self.image_b)

    def test_incorrect_input_compute_intensity_asserts(self):
        with pytest.raises(ValueError, match=r"specify either a single image"):
            lupton_rgb.compute_intensity(self.image_r, self.image_g)

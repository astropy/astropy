# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import tempfile

import numpy as np
from numpy.testing import assert_allclose
import pytest

from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB, HAS_PLT
from astropy.visualization import log_linear_rgb

# Set DISPLAY=True to get matplotlib imshow windows to help with debugging.
DISPLAY = False

# Override any debugging if not HAS_PLT
if HAS_PLT & DISPLAY:
    import matplotlib.pyplot as plt
elif not HAS_PLT:
    DISPLAY = False

MINSC = 0.0
MAXSC = 5.0e4
MIN = [0.0, 0.0, 0.0]
MAX = [5.0e4, 5.0e4, 4.0e4]
IX = 16
IY = 16

SCALEA = 1500.0
SHAPE = (85, 75)

# Gaussian pixel centers, peak values, and colors
_points = [[15, 15], [50, 45], [30, 30], [45, 15]]
_values = [1000, 5500, 600, 20000]
_stddevs = [1.5, 2.0, 1.0, 2.5]
_g_r = [1.0, -1.0, 1.0, 1.0]
_r_i = [2.0, -0.5, 2.5, 1.0]

_grid = np.indices(SHAPE)
IMAGER = np.zeros(SHAPE)
IMAGEG = np.zeros(SHAPE)
IMAGEB = np.zeros(SHAPE)
for p, v, std, gr, ri in zip(_points, _values, _stddevs, _g_r, _r_i):
    _gaus = np.exp(
        -(
            (_grid[0] - p[0]) ** 2 / (2.0 * std**2)
            + (_grid[1] - p[1]) ** 2 / (2.0 * std**2)
        )
    )
    IMAGER += v * np.power(10, 0.4 * ri) * _gaus
    IMAGEG += v * np.power(10, 0.4 * gr) * _gaus
    IMAGEB += v * _gaus

RSEED = 0
rng = np.random.default_rng(RSEED)
IMAGER = IMAGER + rng.normal(0, 2, SHAPE)
IMAGEG = IMAGEG + rng.normal(0, 2, SHAPE)
IMAGEB = IMAGEB + rng.normal(0, 2, SHAPE)

INCORRECT_OUTPUT_TYPES = [bool, str, np.int64, np.cdouble, "str"]


def _display_rgb(rgb, title=None):
    """Display an rgb image using matplotlib (useful for debugging)"""
    plt.imshow(rgb, interpolation="nearest", origin="lower")
    if title:
        plt.title(title)
    plt.show()
    return plt


def test_image_mapping():
    """Test creating an RGB image using a linear stretch,
    using RGBImageMapping()"""
    stretch = log_linear_rgb.LinearStretch()
    map_ = log_linear_rgb.RGBImageMapping(MIN, MAX, stretch)
    rgb_image = map_.make_rgb_image(
        IMAGER, IMAGEG, IMAGEB, output_image_format=np.float64
    )
    for i, (min_, max_, iref_) in enumerate(
        zip(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.5000598388671327],
            [0.08093024185629245, 0.032216094791227695, 0.016040737174622725],
        )
    ):
        assert_allclose(rgb_image[:, :, i].min(), min_)
        assert_allclose(rgb_image[:, :, i].max(), max_)
        assert_allclose(rgb_image[IX, IY, i], iref_)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_linear():
    """Test creating an RGB image using a linear stretch,
    using individual routines"""
    rgb_image = log_linear_rgb.make_linear_rgb(
        IMAGER, IMAGEG, IMAGEB, MIN, MAX, output_image_format=np.float64
    )
    for i, (min_, max_, iref_) in enumerate(
        zip(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.5000598388671327],
            [0.08093024185629245, 0.032216094791227695, 0.016040737174622725],
        )
    ):
        assert_allclose(rgb_image[:, :, i].min(), min_)
        assert_allclose(rgb_image[:, :, i].max(), max_)
        assert_allclose(rgb_image[IX, IY, i], iref_)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_log():
    """Test creating an RGB image using an log stretch"""
    rgb_image = log_linear_rgb.make_log_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        MIN,
        MAX,
        SCALEA,
        output_image_format=np.float64,
    )
    for i, (min_, max_, iref_) in enumerate(
        zip(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.9053360156408082],
            [0.6572779418489928, 0.5330153105260111, 0.4404384627801792],
        )
    ):
        assert_allclose(rgb_image[:, :, i].min(), min_)
        assert_allclose(rgb_image[:, :, i].max(), max_)
        assert_allclose(rgb_image[IX, IY, i], iref_)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_int8():
    """Test creating an RGB image with 8-bit output format"""
    rgb_image = log_linear_rgb.make_log_rgb(
        IMAGER, IMAGEG, IMAGEB, MIN, MAX, SCALEA, output_image_format=np.uint8
    )
    assert np.issubdtype(rgb_image.dtype, np.uint8)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_float64():
    """Test creating an RGB image with normalized float output format"""
    rgb_image = log_linear_rgb.make_log_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        MIN,
        MAX,
        SCALEA,
        output_image_format=np.float64,
    )
    assert np.issubdtype(rgb_image.dtype, float)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_linear_min_max():
    """Test using a min/max linear stretch determined from one image"""
    rgb_image = log_linear_rgb.make_linear_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        None,
        None,
        output_image_format=np.float64,
    )
    for i, (min_, max_, iref_) in enumerate(
        zip(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.08069534125307666, 0.032196043103128555, 0.032466842729915714],
        )
    ):
        assert_allclose(rgb_image[:, :, i].min(), min_)
        assert_allclose(rgb_image[:, :, i].max(), max_)
        assert_allclose(rgb_image[IX, IY, i], iref_)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_log_min_max():
    """Test using a min/max log stretch determined from one image"""
    rgb_image = log_linear_rgb.make_log_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        None,
        None,
        SCALEA,
        output_image_format=np.float64,
    )
    for i, (min_, max_, iref_) in enumerate(
        zip(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.6568837677677257, 0.5329319103684619, 0.5340539629318083],
        )
    ):
        assert_allclose(rgb_image[:, :, i].min(), min_)
        assert_allclose(rgb_image[:, :, i].max(), max_)
        assert_allclose(rgb_image[IX, IY, i], iref_)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_log_scalar_interval():
    """Test creating a black+white image using a linear stretch"""
    rgb_image = log_linear_rgb.make_log_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        MINSC,
        MAXSC,
        SCALEA,
        output_image_format=np.float64,
    )
    for i, (min_, max_, iref_) in enumerate(
        zip(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.8748719461075388],
            [0.6572779418489928, 0.5330153105260111, 0.41128606174423155],
        )
    ):
        assert_allclose(rgb_image[:, :, i].min(), min_)
        assert_allclose(rgb_image[:, :, i].max(), max_)
        assert_allclose(rgb_image[IX, IY, i], iref_)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_linear_bw():
    """Test creating a black+white image using a linear stretch"""
    rgb_image = log_linear_rgb.make_linear_rgb(
        IMAGER,
        IMAGER,
        IMAGER,
        MINSC,
        MAXSC,
        output_image_format=np.float64,
    )
    for i, (min_, max_, iref_) in enumerate(
        zip(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.08093024185629245, 0.08093024185629245, 0.08093024185629245],
        )
    ):
        assert_allclose(rgb_image[:, :, i].min(), min_)
        assert_allclose(rgb_image[:, :, i].max(), max_)
        assert_allclose(rgb_image[IX, IY, i], iref_)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_log_bw():
    """Test creating a black+white image using a log stretch"""
    rgb_image = log_linear_rgb.make_log_rgb(
        IMAGER,
        IMAGER,
        IMAGER,
        MINSC,
        MAXSC,
        SCALEA,
        output_image_format=np.float64,
    )
    for i, (min_, max_, iref_) in enumerate(
        zip(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.6572779418489928, 0.6572779418489928, 0.6572779418489928],
        )
    ):
        assert_allclose(rgb_image[:, :, i].min(), min_)
        assert_allclose(rgb_image[:, :, i].max(), max_)
        assert_allclose(rgb_image[IX, IY, i], iref_)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="requires matplotlib")
def test_make_log_rgb_file():
    """Test the function that does it all"""
    with tempfile.NamedTemporaryFile(suffix=".png") as temp:
        red = IMAGER
        green = IMAGEG
        blue = IMAGEB
        log_linear_rgb.make_log_rgb(red, green, blue, MIN, MAX, SCALEA, filename=temp)
        assert os.path.exists(temp.name)


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="requires matplotlib")
def test_make_linear_rgb_file():
    """Test the function that does it all"""
    with tempfile.NamedTemporaryFile(suffix=".png") as temp:
        red = IMAGER
        green = IMAGEG
        blue = IMAGEB
        log_linear_rgb.make_linear_rgb(red, green, blue, MIN, MAX, filename=temp)
        assert os.path.exists(temp.name)


def test_different_shapes_asserts():
    """Test the different shape assertion"""
    with pytest.raises(ValueError, match=r"shapes must match"):
        # just swap the dimensions to get a differently-shaped 'r'
        image_r = IMAGER.reshape(SHAPE[1], SHAPE[0])
        log_linear_rgb.make_log_rgb(image_r, IMAGEG, IMAGEB)


def test_incorrect_min_length():
    """Test incorrect input minimum array length"""
    with pytest.raises(ValueError, match=r"or 3 values for minimum"):
        log_linear_rgb.make_log_rgb(IMAGER, IMAGEG, IMAGEB, [MINSC, MINSC], MAX, SCALEA)


def test_incorrect_max_length():
    """Test incorrect input maximum array length"""
    with pytest.raises(ValueError, match=r"or 3 values for maximum"):
        log_linear_rgb.make_log_rgb(IMAGER, IMAGEG, IMAGEB, MAX, [MINSC, MINSC], SCALEA)


@pytest.mark.parametrize(("out_format"), INCORRECT_OUTPUT_TYPES)
def test_invalid_output_image_format(out_format):
    """Test incorrect output image format"""
    with pytest.raises(ValueError, match=r"'output_image_format' must be one"):
        log_linear_rgb.make_log_rgb(
            IMAGER,
            IMAGEG,
            IMAGEB,
            MIN,
            MAX,
            SCALEA,
            output_image_format=out_format,
        )

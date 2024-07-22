# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import tempfile

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB, HAS_PLT
from astropy.visualization import basic_rgb
from astropy.visualization.interval import ManualInterval
from astropy.visualization.stretch import LinearStretch, LogStretch

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
    stretch = LinearStretch()
    interval = []
    for i in range(3):
        interval.append(ManualInterval(vmin=MIN[i], vmax=MAX[i]))
    map_ = basic_rgb.RGBImageMapping(stretch=stretch, interval=interval)
    rgb_image = map_.make_rgb_image(IMAGER, IMAGEG, IMAGEB, output_dtype=np.float64)
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
    interval = []
    for i in range(3):
        interval.append(ManualInterval(vmin=MIN[i], vmax=MAX[i]))
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        interval=interval,
        output_dtype=np.float64,
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
    interval = []
    for i in range(3):
        interval.append(ManualInterval(vmin=MIN[i], vmax=MAX[i]))
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        interval=interval,
        stretch=LogStretch(a=SCALEA),
        output_dtype=np.float64,
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
    interval = []
    for i in range(3):
        interval.append(ManualInterval(vmin=MIN[i], vmax=MAX[i]))
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        interval=interval,
        stretch=LogStretch(a=SCALEA),
        output_dtype=np.uint8,
    )
    assert np.issubdtype(rgb_image.dtype, np.uint8)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_float64():
    """Test creating an RGB image with normalized float output format"""
    interval = []
    for i in range(3):
        interval.append(ManualInterval(vmin=MIN[i], vmax=MAX[i]))
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        interval=interval,
        stretch=LogStretch(a=SCALEA),
        output_dtype=np.float64,
    )
    assert np.issubdtype(rgb_image.dtype, float)
    if DISPLAY:
        _display_rgb(rgb_image, title=sys._getframe().f_code.co_name)


def test_linear_min_max():
    """Test using a min/max linear stretch determined from one image"""
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        interval=ManualInterval(vmin=None, vmax=None),
        output_dtype=np.float64,
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
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        interval=ManualInterval(vmin=None, vmax=None),
        stretch=LogStretch(a=SCALEA),
        output_dtype=np.float64,
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
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGEG,
        IMAGEB,
        interval=ManualInterval(vmin=MINSC, vmax=MAXSC),
        stretch=LogStretch(a=SCALEA),
        output_dtype=np.float64,
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
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGER,
        IMAGER,
        interval=ManualInterval(vmin=MINSC, vmax=MAXSC),
        output_dtype=np.float64,
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
    rgb_image = basic_rgb.make_rgb(
        IMAGER,
        IMAGER,
        IMAGER,
        interval=ManualInterval(vmin=MINSC, vmax=MAXSC),
        stretch=LogStretch(a=SCALEA),
        output_dtype=np.float64,
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
    interval = []
    for i in range(3):
        interval.append(ManualInterval(vmin=MIN[i], vmax=MAX[i]))
    with tempfile.NamedTemporaryFile(suffix=".png") as temp:
        red = IMAGER
        green = IMAGEG
        blue = IMAGEB
        basic_rgb.make_rgb(
            red,
            green,
            blue,
            interval=interval,
            stretch=LogStretch(a=SCALEA),
            filename=temp,
        )
        assert os.path.exists(temp.name)


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="requires matplotlib")
def test_make_linear_rgb_file():
    """Test the function that does it all"""
    interval = []
    for i in range(3):
        interval.append(ManualInterval(vmin=MIN[i], vmax=MAX[i]))
    with tempfile.NamedTemporaryFile(suffix=".png") as temp:
        red = IMAGER
        green = IMAGEG
        blue = IMAGEB
        basic_rgb.make_rgb(red, green, blue, interval=interval, filename=temp)
        assert os.path.exists(temp.name)


def test_different_shapes_asserts():
    """Test the different shape assertion"""
    with pytest.raises(ValueError, match=r"shapes must match"):
        # just swap the dimensions to get a differently-shaped 'r'
        image_r = IMAGER.reshape(SHAPE[1], SHAPE[0])
        basic_rgb.make_rgb(image_r, IMAGEG, IMAGEB, stretch=LogStretch(a=SCALEA))


def test_incorrect_interval_length():
    """Test incorrect input interval array length"""
    with pytest.raises(ValueError, match=r"3 instances for interval."):
        interval = ManualInterval(vmin=MINSC, vmax=MAXSC)
        basic_rgb.make_rgb(
            IMAGER,
            IMAGEG,
            IMAGEB,
            interval=[interval, interval],
            stretch=LogStretch(a=SCALEA),
        )


@pytest.mark.parametrize(("out_format"), INCORRECT_OUTPUT_TYPES)
def test_invalid_output_dtype(out_format):
    """Test incorrect output image format"""
    interval = []
    for i in range(3):
        interval.append(ManualInterval(vmin=MIN[i], vmax=MAX[i]))
    with pytest.raises(ValueError, match=r"'output_dtype' must be one"):
        basic_rgb.make_rgb(
            IMAGER,
            IMAGEG,
            IMAGEB,
            interval=interval,
            stretch=LogStretch(a=SCALEA),
            output_dtype=out_format,
        )

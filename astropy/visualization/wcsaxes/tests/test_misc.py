# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings
from contextlib import nullcontext

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.contour import QuadContourSet
from packaging.version import Version

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization.wcsaxes.core import WCSAxes
from astropy.visualization.wcsaxes.frame import (
    EllipticalFrame,
    RectangularFrame,
    RectangularFrame1D,
)
from astropy.visualization.wcsaxes.ticklabels import TickLabels
from astropy.visualization.wcsaxes.transforms import CurvedTransform
from astropy.visualization.wcsaxes.utils import get_coord_meta
from astropy.wcs import WCS
from astropy.wcs.wcsapi import HighLevelWCSWrapper, SlicedLowLevelWCS

ft_version = Version(mpl.ft2font.__freetype_version__)
FREETYPE_261 = ft_version == Version("2.6.1")

# We cannot use matplotlib.checkdep_usetex() anymore, see
# https://github.com/matplotlib/matplotlib/issues/23244
TEX_UNAVAILABLE = True

MATPLOTLIB_LT_3_7 = Version(mpl.__version__) < Version("3.7")


def teardown_function(function):
    plt.close("all")


def test_grid_regression(ignore_matplotlibrc):
    # Regression test for a bug that meant that if the rc parameter
    # axes.grid was set to True, WCSAxes would crash upon initialization.
    plt.rc("axes", grid=True)
    fig = plt.figure(figsize=(3, 3))
    WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])


def test_format_coord_regression(ignore_matplotlibrc, tmp_path):
    # Regression test for a bug that meant that if format_coord was called by
    # Matplotlib before the axes were drawn, an error occurred.
    fig = plt.figure(figsize=(3, 3))
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)
    assert ax.format_coord(10, 10) == ""
    assert ax.coords[0].format_coord(10) == ""
    assert ax.coords[1].format_coord(10) == ""
    fig.savefig(tmp_path / "nothing")
    assert ax.format_coord(10, 10) == "10.0 10.0 (world)"
    assert ax.coords[0].format_coord(10) == "10.0"
    assert ax.coords[1].format_coord(10) == "10.0"


TARGET_HEADER = fits.Header.fromstring(
    """
NAXIS   =                    2
NAXIS1  =                  200
NAXIS2  =                  100
CTYPE1  = 'RA---MOL'
CRPIX1  =                  500
CRVAL1  =                180.0
CDELT1  =                 -0.4
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--MOL'
CRPIX2  =                  400
CRVAL2  =                  0.0
CDELT2  =                  0.4
CUNIT2  = 'deg     '
COORDSYS= 'icrs    '
""",
    sep="\n",
)


@pytest.mark.parametrize("grid_type", ["lines", "contours"])
def test_no_numpy_warnings(ignore_matplotlibrc, tmp_path, grid_type):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS(TARGET_HEADER))
    ax.imshow(np.zeros((100, 200)))
    ax.coords.grid(color="white", grid_type=grid_type)

    # There should be no warnings raised if some pixels are outside WCS
    # (since this is normal).
    # BUT our own catch_warning was ignoring some warnings before, so now we
    # have to catch it. Otherwise, the pytest filterwarnings=error
    # setting in pyproject.toml will fail this test.
    # There are actually multiple warnings but they are all similar.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=r".*converting a masked element to nan.*"
        )
        warnings.filterwarnings(
            "ignore", message=r".*No contour levels were found within the data range.*"
        )
        warnings.filterwarnings(
            "ignore", message=r".*np\.asscalar\(a\) is deprecated since NumPy v1\.16.*"
        )
        warnings.filterwarnings(
            "ignore", message=r".*PY_SSIZE_T_CLEAN will be required.*"
        )
        fig.savefig(tmp_path / "test.png")


def test_invalid_frame_overlay(ignore_matplotlibrc):
    # Make sure a nice error is returned if a frame doesn't exist
    ax = plt.subplot(1, 1, 1, projection=WCS(TARGET_HEADER))
    with pytest.raises(ValueError, match=r"Frame banana not found"):
        ax.get_coords_overlay("banana")

    with pytest.raises(ValueError, match=r"Unknown frame: banana"):
        get_coord_meta("banana")


def test_plot_coord_transform(ignore_matplotlibrc):
    twoMASS_k_header = get_pkg_data_filename("data/2MASS_k_header")
    twoMASS_k_header = fits.Header.fromtextfile(twoMASS_k_header)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes(
        [0.15, 0.15, 0.8, 0.8], projection=WCS(twoMASS_k_header), aspect="equal"
    )
    ax.set_xlim(-0.5, 720.5)
    ax.set_ylim(-0.5, 720.5)

    c = SkyCoord(359.76045223 * u.deg, 0.26876217 * u.deg)
    with pytest.raises(TypeError):
        ax.plot_coord(c, "o", transform=ax.get_transform("galactic"))


def test_scatter_coord_transform(ignore_matplotlibrc):
    twoMASS_k_header = get_pkg_data_filename("data/2MASS_k_header")
    twoMASS_k_header = fits.Header.fromtextfile(twoMASS_k_header)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes(
        [0.15, 0.15, 0.8, 0.8], projection=WCS(twoMASS_k_header), aspect="equal"
    )
    ax.set_xlim(-0.5, 720.5)
    ax.set_ylim(-0.5, 720.5)

    c = SkyCoord(359.76045223 * u.deg, 0.26876217 * u.deg)
    with pytest.raises(TypeError):
        ax.scatter_coord(c, marker="o", transform=ax.get_transform("galactic"))


def test_set_label_properties(ignore_matplotlibrc):
    # Regression test to make sure that arguments passed to
    # set_xlabel/set_ylabel are passed to the underlying coordinate helpers

    ax = plt.subplot(1, 1, 1, projection=WCS(TARGET_HEADER))

    ax.set_xlabel("Test x label", labelpad=2, color="red")
    ax.set_ylabel("Test y label", labelpad=3, color="green")

    assert ax.coords[0].axislabels.get_text() == "Test x label"
    assert ax.coords[0].axislabels.get_minpad("b") == 2
    assert ax.coords[0].axislabels.get_color() == "red"

    assert ax.coords[1].axislabels.get_text() == "Test y label"
    assert ax.coords[1].axislabels.get_minpad("l") == 3
    assert ax.coords[1].axislabels.get_color() == "green"

    assert ax.get_xlabel() == "Test x label"
    assert ax.get_ylabel() == "Test y label"


GAL_HEADER = fits.Header.fromstring(
    """
SIMPLE  =                    T / conforms to FITS standard
BITPIX  =                  -32 / array data type
NAXIS   =                    3 / number of array dimensions
NAXIS1  =                   31
NAXIS2  =                 2881
NAXIS3  =                  480
EXTEND  =                    T
CTYPE1  = 'DISTMOD '
CRVAL1  =                  3.5
CDELT1  =                  0.5
CRPIX1  =                  1.0
CTYPE2  = 'GLON-CAR'
CRVAL2  =                180.0
CDELT2  =               -0.125
CRPIX2  =                  1.0
CTYPE3  = 'GLAT-CAR'
CRVAL3  =                  0.0
CDELT3  =                0.125
CRPIX3  =                241.0
""",
    sep="\n",
)


def test_slicing_warnings(ignore_matplotlibrc, tmp_path):
    # Regression test to make sure that no warnings are emitted by the tick
    # locator for the sliced axis when slicing a cube.

    # Scalar case

    wcs3d = WCS(naxis=3)
    wcs3d.wcs.ctype = ["x", "y", "z"]
    wcs3d.wcs.cunit = ["deg", "deg", "km/s"]
    wcs3d.wcs.crpix = [614.5, 856.5, 333]
    wcs3d.wcs.cdelt = [6.25, 6.25, 23]
    wcs3d.wcs.crval = [0.0, 0.0, 1.0]

    with warnings.catch_warnings():
        # https://github.com/astropy/astropy/issues/9690
        warnings.filterwarnings("ignore", message=r".*PY_SSIZE_T_CLEAN.*")
        plt.subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1))
        plt.savefig(tmp_path / "test.png")

    # Angle case

    wcs3d = WCS(GAL_HEADER)

    with warnings.catch_warnings():
        # https://github.com/astropy/astropy/issues/9690
        warnings.filterwarnings("ignore", message=r".*PY_SSIZE_T_CLEAN.*")
        plt.clf()
        plt.subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 2))
        plt.savefig(tmp_path / "test.png")


def test_plt_xlabel_ylabel(tmp_path):
    # Regression test for a bug that happened when using plt.xlabel
    # and plt.ylabel with Matplotlib 3.0

    plt.subplot(projection=WCS())
    plt.xlabel("Galactic Longitude")
    plt.ylabel("Galactic Latitude")
    plt.savefig(tmp_path / "test.png")


def test_grid_type_contours_transform(tmp_path):
    # Regression test for a bug that caused grid_type='contours' to not work
    # with custom transforms

    class CustomTransform(CurvedTransform):
        # We deliberately don't define the inverse, and has_inverse should
        # default to False.

        def transform(self, values):
            return values * 1.3

    transform = CustomTransform()
    coord_meta = {
        "type": ("scalar", "scalar"),
        "unit": (u.m, u.s),
        "wrap": (None, None),
        "name": ("x", "y"),
    }

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], transform=transform, coord_meta=coord_meta)
    fig.add_axes(ax)
    ax.grid(grid_type="contours")
    fig.savefig(tmp_path / "test.png")


def test_plt_imshow_origin():
    # Regression test for a bug that caused origin to be set to upper when
    # plt.imshow was called.

    ax = plt.subplot(projection=WCS())
    plt.imshow(np.ones((2, 2)))
    assert ax.get_xlim() == (-0.5, 1.5)
    assert ax.get_ylim() == (-0.5, 1.5)


def test_ax_imshow_origin():
    # Regression test for a bug that caused origin to be set to upper when
    # ax.imshow was called with no origin

    ax = plt.subplot(projection=WCS())
    ax.imshow(np.ones((2, 2)))
    assert ax.get_xlim() == (-0.5, 1.5)
    assert ax.get_ylim() == (-0.5, 1.5)


def test_grid_contour_large_spacing(tmp_path):
    # Regression test for a bug that caused a crash when grid was called and
    # didn't produce grid lines (due e.g. to too large spacing) and was then
    # called again.

    filename = tmp_path / "test.png"

    ax = plt.subplot(projection=WCS())
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.5, 1.5)
    ax.coords[0].set_ticks(values=[] * u.one)

    ax.coords[0].grid(grid_type="contours")
    plt.savefig(filename)

    ax.coords[0].grid(grid_type="contours")
    plt.savefig(filename)


def test_contour_return():
    # Regression test for a bug that caused contour and contourf to return None
    # instead of the contour object.

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)

    cset = ax.contour(np.arange(16).reshape(4, 4), transform=ax.get_transform("world"))
    assert isinstance(cset, QuadContourSet)

    cset = ax.contourf(np.arange(16).reshape(4, 4), transform=ax.get_transform("world"))
    assert isinstance(cset, QuadContourSet)


def test_contour_empty():
    # Regression test for a bug that caused contour to crash if no contours
    # were present.

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)

    if MATPLOTLIB_LT_3_7:
        ctx = pytest.warns(
            UserWarning, match="No contour levels were found within the data range"
        )
    else:
        ctx = nullcontext()

    with ctx:
        ax.contour(np.zeros((4, 4)), transform=ax.get_transform("world"))


def test_iterate_coords(ignore_matplotlibrc):
    # Regression test for a bug that caused ax.coords to return too few axes

    wcs3d = WCS(naxis=3)
    wcs3d.wcs.ctype = ["x", "y", "z"]
    wcs3d.wcs.cunit = ["deg", "deg", "km/s"]
    wcs3d.wcs.crpix = [614.5, 856.5, 333]
    wcs3d.wcs.cdelt = [6.25, 6.25, 23]
    wcs3d.wcs.crval = [0.0, 0.0, 1.0]

    ax = plt.subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1))

    x, y, z = ax.coords


def test_invalid_slices_errors(ignore_matplotlibrc):
    # Make sure that users get a clear message when specifying a WCS with
    # >2 dimensions without giving the 'slices' argument, or if the 'slices'
    # argument has too many/few elements.

    wcs3d = WCS(naxis=3)
    wcs3d.wcs.ctype = ["x", "y", "z"]

    plt.subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1))

    with pytest.raises(
        ValueError,
        match=r"WCS has more than 2 pixel dimensions, so 'slices' should be set",
    ):
        plt.subplot(1, 1, 1, projection=wcs3d)

    with pytest.raises(
        ValueError,
        match=(
            r"'slices' should have as many elements as WCS has pixel dimensions .should"
            r" be 3."
        ),
    ):
        plt.subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1, 2))

    wcs2d = WCS(naxis=2)
    wcs2d.wcs.ctype = ["x", "y"]

    plt.clf()
    ax = plt.subplot(1, 1, 1, projection=wcs2d)
    assert ax.frame_class is RectangularFrame
    plt.clf()
    ax = plt.subplot(1, 1, 1, projection=wcs2d, slices=("x", "y"))
    assert ax.frame_class is RectangularFrame
    plt.clf()
    ax = plt.subplot(1, 1, 1, projection=wcs2d, slices=("y", "x"))
    assert ax.frame_class is RectangularFrame
    plt.clf()
    ax = plt.subplot(1, 1, 1, projection=wcs2d, slices=["x", "y"])
    assert ax.frame_class is RectangularFrame
    plt.clf()
    ax = plt.subplot(1, 1, 1, projection=wcs2d, slices=(1, "x"))
    assert ax.frame_class is RectangularFrame1D

    wcs1d = WCS(naxis=1)
    wcs1d.wcs.ctype = ["x"]

    plt.clf()
    ax = plt.subplot(1, 1, 1, projection=wcs1d)
    assert ax.frame_class is RectangularFrame1D

    with pytest.raises(ValueError):
        plt.subplot(1, 1, 1, projection=wcs2d, slices=(1, "y"))


EXPECTED_REPR_1 = """
<CoordinatesMap with 3 world coordinates:

  index            aliases                type   ...    wrap   format_unit visible
  ----- ------------------------------ --------- ... --------- ----------- -------
      0                   distmod dist    scalar ...      None                  no
      1 pos.galactic.lon glon-car glon longitude ... 360.0 deg         deg     yes
      2 pos.galactic.lat glat-car glat  latitude ...      None         deg     yes

>
 """.strip()

EXPECTED_REPR_2 = """
<CoordinatesMap with 3 world coordinates:

  index            aliases                type   ...    wrap   format_unit visible
  ----- ------------------------------ --------- ... --------- ----------- -------
      0                   distmod dist    scalar ...      None                 yes
      1 pos.galactic.lon glon-car glon longitude ... 360.0 deg         deg     yes
      2 pos.galactic.lat glat-car glat  latitude ...      None         deg     yes

>
 """.strip()


def test_repr(ignore_matplotlibrc):
    # Unit test to make sure __repr__ looks as expected

    wcs3d = WCS(GAL_HEADER)

    # Cube header has world coordinates as distance, lon, lat, so start off
    # by slicing in a way that we select just lon,lat:

    ax = plt.subplot(1, 1, 1, projection=wcs3d, slices=(1, "x", "y"))
    assert repr(ax.coords) == EXPECTED_REPR_1

    # Now slice in a way that all world coordinates are still present:

    plt.clf()
    ax = plt.subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1))
    assert repr(ax.coords) == EXPECTED_REPR_2


@pytest.fixture
def time_spectral_wcs_2d():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["FREQ", "TIME"]
    wcs.wcs.set()
    return wcs


def test_time_wcs(time_spectral_wcs_2d):
    # Regression test for a bug that caused WCSAxes to error when using a WCS
    # with a time axis.

    plt.subplot(projection=time_spectral_wcs_2d)


@pytest.mark.skipif(TEX_UNAVAILABLE, reason="TeX is unavailable")
def test_simplify_labels_usetex(ignore_matplotlibrc, tmp_path):
    """Regression test for https://github.com/astropy/astropy/issues/8004."""
    plt.rc("text", usetex=True)

    header = {
        "NAXIS": 2,
        "NAXIS1": 360,
        "NAXIS2": 180,
        "CRPIX1": 180.5,
        "CRPIX2": 90.5,
        "CRVAL1": 180.0,
        "CRVAL2": 0.0,
        "CDELT1": -2 * np.sqrt(2) / np.pi,
        "CDELT2": 2 * np.sqrt(2) / np.pi,
        "CTYPE1": "RA---MOL",
        "CTYPE2": "DEC--MOL",
        "RADESYS": "ICRS",
    }

    wcs = WCS(header)
    fig, ax = plt.subplots(subplot_kw=dict(frame_class=EllipticalFrame, projection=wcs))
    ax.set_xlim(-0.5, header["NAXIS1"] - 0.5)
    ax.set_ylim(-0.5, header["NAXIS2"] - 0.5)
    ax.coords[0].set_ticklabel(exclude_overlapping=True)
    ax.coords[1].set_ticklabel(exclude_overlapping=True)
    ax.coords[0].set_ticks(spacing=45 * u.deg)
    ax.coords[1].set_ticks(spacing=30 * u.deg)
    ax.grid()

    fig.savefig(tmp_path / "plot.png")


@pytest.mark.parametrize(
    "usetex, unicode_minus, label_str",
    [
        (True, True, "$-{}$"),
        (True, False, "$-{}$"),
        (False, True, "\N{MINUS SIGN}{}"),
        (False, False, "-{}"),
    ],
)
def test_simplify_labels_minus_sign(
    ignore_matplotlibrc, usetex, unicode_minus, label_str
):
    # Ensure minus signs aren't removed from the front of labels across a grid of configuration possibilities
    if usetex and TEX_UNAVAILABLE:
        pytest.skip("TeX is unavailable")

    ticklabels = TickLabels(None)
    expected_labels = []
    for i in range(1, 6):
        label = label_str.format(i)
        ticklabels.add(
            axis="axis",
            world=0,
            angle=0,
            text=label,
            axis_displacement=0,
            data=(i, i),
        )
        expected_labels.append(label)

    with mpl.rc_context(
        rc={"text.usetex": usetex, "axes.unicode_minus": unicode_minus}
    ):
        ticklabels.simplify_labels()
    assert ticklabels.text["axis"] == expected_labels


@pytest.mark.parametrize("frame_class", [RectangularFrame, EllipticalFrame])
def test_set_labels_with_coords(ignore_matplotlibrc, frame_class):
    """Test if ``axis.set_xlabel()`` calls the correct ``coords[i]_set_axislabel()`` in a
    WCS plot. Regression test for https://github.com/astropy/astropy/issues/10435.
    """

    labels = ["RA", "Declination"]
    header = {
        "NAXIS": 2,
        "NAXIS1": 360,
        "NAXIS2": 180,
        "CRPIX1": 180.5,
        "CRPIX2": 90.5,
        "CRVAL1": 180.0,
        "CRVAL2": 0.0,
        "CDELT1": -2 * np.sqrt(2) / np.pi,
        "CDELT2": 2 * np.sqrt(2) / np.pi,
        "CTYPE1": "RA---AIT",
        "CTYPE2": "DEC--AIT",
    }

    wcs = WCS(header)
    fig, ax = plt.subplots(subplot_kw=dict(frame_class=frame_class, projection=wcs))
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])

    assert ax.get_xlabel() == labels[0]
    assert ax.get_ylabel() == labels[1]
    for i in range(2):
        assert ax.coords[i].get_axislabel() == labels[i]


@pytest.mark.parametrize("atol", [0.2, 1.0e-8])
def test_bbox_size(atol):
    # Test for the size of a WCSAxes bbox (only have Matplotlib >= 3.0 now)
    extents = [11.38888888888889, 3.5, 576.0, 432.0]

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)
    fig.canvas.draw()
    renderer = fig.canvas.renderer
    ax_bbox = ax.get_tightbbox(renderer)

    # Enforce strict test only with reference Freetype version
    if atol < 0.1 and not FREETYPE_261:
        pytest.xfail(
            "Exact BoundingBox dimensions are only ensured with FreeType 2.6.1"
        )
    assert np.allclose(ax_bbox.extents, extents, atol=atol)


def test_wcs_type_transform_regression():
    wcs = WCS(TARGET_HEADER)
    sliced_wcs = SlicedLowLevelWCS(wcs, np.s_[1:-1, 1:-1])
    ax = plt.subplot(1, 1, 1, projection=wcs)
    ax.get_transform(sliced_wcs)

    high_wcs = HighLevelWCSWrapper(sliced_wcs)
    ax.get_transform(sliced_wcs)


def test_multiple_draws_grid_contours(tmp_path):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS())
    ax.grid(color="black", grid_type="contours")
    fig.savefig(tmp_path / "plot.png")
    fig.savefig(tmp_path / "plot.png")


def test_get_coord_range_nan_regression():
    # Test to make sure there is no internal casting of NaN to integers
    # NumPy 1.24 raises a RuntimeWarning if a NaN is cast to an integer

    wcs = WCS(TARGET_HEADER)
    wcs.wcs.crval[0] = 0  # Re-position the longitude wrap to the middle
    ax = plt.subplot(1, 1, 1, projection=wcs)

    # Set the Y limits within valid latitudes/declinations
    ax.set_ylim(300, 500)

    # Set the X limits within valid longitudes/RAs, so the world coordinates have no NaNs
    ax.set_xlim(300, 700)
    assert np.allclose(
        ax.coords.get_coord_range(),
        np.array(
            [
                (-123.5219272110385, 122.49684897692201),
                (-44.02289164685554, 44.80732766607591),
            ]
        ),
    )

    # Extend the X limits to include invalid longitudes/RAs, so the world coordinates have NaNs
    ax.set_xlim(0, 700)
    assert np.allclose(
        ax.coords.get_coord_range(),
        np.array(
            [(-131.3193386797236, 180.0), (-44.02289164685554, 44.80732766607591)]
        ),
    )


def test_imshow_error():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS())
    with pytest.raises(ValueError, match="Cannot use images with origin='upper"):
        ax.imshow(np.ones(100).reshape(10, 10), origin="upper")


def test_label_setting():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS())
    # Check both xlabel and label kwargs work
    ax.set_xlabel(xlabel="label")
    ax.set_xlabel(label="label")
    # Check no label errors:
    with pytest.raises(
        TypeError, match=r"set_xlabel\(\) missing 1 required positional argument"
    ):
        ax.set_xlabel()

    # Check both xlabel and label kwargs work
    ax.set_ylabel(ylabel="label")
    ax.set_ylabel(label="label")
    # Check no label errors:
    with pytest.raises(
        TypeError, match=r"set_ylabel\(\) missing 1 required positional argument"
    ):
        ax.set_ylabel()


def test_invisible_bbox():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS())

    assert ax.get_tightbbox(fig.canvas.get_renderer()) is not None
    ax.set_visible(False)
    assert ax.get_tightbbox(fig.canvas.get_renderer()) is None

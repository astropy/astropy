# Licensed under a 3-clause BSD style license - see LICENSE.rst
import inspect
import warnings
from unittest.mock import MagicMock

import matplotlib as mpl
import numpy as np
import pytest
from matplotlib.backend_bases import KeyEvent
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.contour import QuadContourSet
from matplotlib.figure import Figure
from matplotlib.transforms import Affine2D
from numpy.testing import assert_allclose
from packaging.version import Version

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.visualization import quantity_support
from astropy.visualization.wcsaxes.core import WCSAxes
from astropy.visualization.wcsaxes.frame import (
    EllipticalFrame,
    RectangularFrame,
    RectangularFrame1D,
)
from astropy.visualization.wcsaxes.helpers import add_beam, add_scalebar
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


def test_grid_regression(ignore_matplotlibrc):
    # Regression test for a bug that meant that if the rc parameter
    # axes.grid was set to True, WCSAxes would crash upon initialization.
    mpl.rc("axes", grid=True)
    fig = Figure(figsize=(3, 3))
    WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])


def test_grid_visible_keyword(ignore_matplotlibrc):
    # Regression test for a bug that meant that ax.grid(visible=False),
    # using the parameter name from matplotlib.axes.Axes.grid, silently
    # turned the grid on instead of off.
    assert "visible" in inspect.signature(WCSAxes.grid).parameters

    fig = Figure(figsize=(3, 3))
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])

    ax.grid(visible=True)
    assert all(coord._grid_lines_kwargs["visible"] for coord in ax.coords)

    ax.grid(visible=False)
    assert not any(coord._grid_lines_kwargs["visible"] for coord in ax.coords)

    # Passing visible directly to CoordinateHelper.grid should also work
    ax.coords[0].grid(visible=True)
    assert ax.coords[0]._grid_lines_kwargs["visible"]
    ax.coords[0].grid(visible=False)
    assert not ax.coords[0]._grid_lines_kwargs["visible"]

    # The old parameter name should still work but emit a deprecation warning
    with pytest.warns(AstropyDeprecationWarning, match='"b" was deprecated'):
        ax.grid(b=False)
    assert not any(coord._grid_lines_kwargs["visible"] for coord in ax.coords)


def test_format_coord_regression(ignore_matplotlibrc, tmp_path):
    # Regression test for a bug that meant that if format_coord was called by
    # Matplotlib before the axes were drawn, an error occurred.
    fig = Figure(figsize=(3, 3))
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
    fig = Figure()
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
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS(TARGET_HEADER))
    with pytest.raises(ValueError, match=r"Frame banana not found"):
        ax.get_coords_overlay("banana")

    with pytest.raises(ValueError, match=r"Unknown frame: banana"):
        get_coord_meta("banana")


def test_plot_coord_transform(ignore_matplotlibrc):
    twoMASS_k_header = get_pkg_data_filename("data/2MASS_k_header")
    twoMASS_k_header = fits.Header.fromtextfile(twoMASS_k_header)
    fig = Figure(figsize=(6, 6))
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
    fig = Figure(figsize=(6, 6))
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

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS(TARGET_HEADER))

    ax.set_xlabel("Test x label", labelpad=2, color="red")
    ax.set_ylabel("Test y label", labelpad=3, color="green")

    assert ax.coords[0]._axislabels.get_text() == "Test x label"
    assert ax.coords[0]._axislabels.get_minpad("b") == 2
    assert ax.coords[0]._axislabels.get_color() == "red"

    assert ax.coords[1]._axislabels.get_text() == "Test y label"
    assert ax.coords[1]._axislabels.get_minpad("l") == 3
    assert ax.coords[1]._axislabels.get_color() == "green"

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

    fig = Figure()
    with warnings.catch_warnings():
        # https://github.com/astropy/astropy/issues/9690
        warnings.filterwarnings("ignore", message=r".*PY_SSIZE_T_CLEAN.*")
        fig.add_subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1))
        fig.savefig(tmp_path / "test.png")

    # Angle case

    wcs3d = WCS(GAL_HEADER)

    with warnings.catch_warnings():
        # https://github.com/astropy/astropy/issues/9690
        warnings.filterwarnings("ignore", message=r".*PY_SSIZE_T_CLEAN.*")
        fig.clear()
        fig.add_subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 2))
        fig.savefig(tmp_path / "test.png")


@pytest.fixture
def tmp_plt_figure():
    import matplotlib.pyplot as plt

    figure = plt.figure()
    yield figure
    plt.close(figure)


@pytest.mark.usefixtures("tmp_plt_figure")
def test_plt_xlabel_ylabel(tmp_path):
    # Regression test for a bug that happened when using plt.xlabel
    # and plt.ylabel with Matplotlib 3.0
    import matplotlib.pyplot as plt

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

    fig = Figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], transform=transform, coord_meta=coord_meta)
    fig.add_axes(ax)
    ax.grid(grid_type="contours")
    fig.savefig(tmp_path / "test.png")


@pytest.mark.usefixtures("tmp_plt_figure")
def test_plt_imshow_origin():
    # Regression test for a bug that caused origin to be set to upper when
    # plt.imshow was called.
    import matplotlib.pyplot as plt

    ax = plt.subplot(projection=WCS())
    plt.imshow(np.ones((2, 2)))
    assert ax.get_xlim() == (-0.5, 1.5)
    assert ax.get_ylim() == (-0.5, 1.5)


def test_ax_imshow_origin():
    # Regression test for a bug that caused origin to be set to upper when
    # ax.imshow was called with no origin
    fig = Figure()
    ax = fig.add_subplot(projection=WCS())
    ax.imshow(np.ones((2, 2)))
    assert ax.get_xlim() == (-0.5, 1.5)
    assert ax.get_ylim() == (-0.5, 1.5)


def test_grid_contour_large_spacing(tmp_path):
    # Regression test for a bug that caused a crash when grid was called and
    # didn't produce grid lines (due e.g. to too large spacing) and was then
    # called again.

    filename = tmp_path / "test.png"

    fig = Figure()
    ax = fig.add_subplot(projection=WCS())
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.5, 1.5)
    ax.coords[0].set_ticks(values=[] * u.one)

    ax.coords[0].grid(grid_type="contours")
    fig.savefig(filename)

    ax.coords[0].grid(grid_type="contours")
    fig.savefig(filename)


def test_contour_return():
    # Regression test for a bug that caused contour and contourf to return None
    # instead of the contour object.

    fig = Figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)

    cset = ax.contour(np.arange(16).reshape(4, 4), transform=ax.get_transform("world"))
    assert isinstance(cset, QuadContourSet)

    cset = ax.contourf(np.arange(16).reshape(4, 4), transform=ax.get_transform("world"))
    assert isinstance(cset, QuadContourSet)


def test_contour_empty():
    # Regression test for a bug that caused contour to crash if no contours
    # were present.

    fig = Figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)
    ax.contour(np.zeros((4, 4)), transform=ax.get_transform("world"))


def test_iterate_coords(ignore_matplotlibrc):
    # Regression test for a bug that caused ax.coords to return too few axes

    wcs3d = WCS(naxis=3)
    wcs3d.wcs.ctype = ["x", "y", "z"]
    wcs3d.wcs.cunit = ["deg", "deg", "km/s"]
    wcs3d.wcs.crpix = [614.5, 856.5, 333]
    wcs3d.wcs.cdelt = [6.25, 6.25, 23]
    wcs3d.wcs.crval = [0.0, 0.0, 1.0]

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1))

    x, y, z = ax.coords


def test_invalid_slices_errors(ignore_matplotlibrc):
    # Make sure that users get a clear message when specifying a WCS with
    # >2 dimensions without giving the 'slices' argument, or if the 'slices'
    # argument has too many/few elements.
    fig = Figure()
    wcs3d = WCS(naxis=3)
    wcs3d.wcs.ctype = ["x", "y", "z"]

    fig.add_subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1))

    with pytest.raises(
        ValueError,
        match=r"WCS has more than 2 pixel dimensions, so 'slices' should be set",
    ):
        fig.add_subplot(1, 1, 1, projection=wcs3d)

    with pytest.raises(
        ValueError,
        match=(
            r"'slices' should have as many elements as WCS has pixel dimensions .should"
            r" be 3."
        ),
    ):
        fig.add_subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1, 2))

    wcs2d = WCS(naxis=2)
    wcs2d.wcs.ctype = ["x", "y"]

    fig.clear()
    ax = fig.add_subplot(1, 1, 1, projection=wcs2d)
    assert ax.frame_class is RectangularFrame
    fig.clear()
    ax = fig.add_subplot(1, 1, 1, projection=wcs2d, slices=("x", "y"))
    assert ax.frame_class is RectangularFrame
    fig.clear()
    ax = fig.add_subplot(1, 1, 1, projection=wcs2d, slices=("y", "x"))
    assert ax.frame_class is RectangularFrame
    fig.clear()
    ax = fig.add_subplot(1, 1, 1, projection=wcs2d, slices=["x", "y"])
    assert ax.frame_class is RectangularFrame
    fig.clear()
    ax = fig.add_subplot(1, 1, 1, projection=wcs2d, slices=(1, "x"))
    assert ax.frame_class is RectangularFrame1D

    wcs1d = WCS(naxis=1)
    wcs1d.wcs.ctype = ["x"]

    fig.clear()
    ax = fig.add_subplot(1, 1, 1, projection=wcs1d)
    assert ax.frame_class is RectangularFrame1D

    with pytest.raises(ValueError):
        fig.add_subplot(1, 1, 1, projection=wcs2d, slices=(1, "y"))


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
    fig = Figure()
    wcs3d = WCS(GAL_HEADER)

    # Cube header has world coordinates as distance, lon, lat, so start off
    # by slicing in a way that we select just lon,lat:

    ax = fig.add_subplot(1, 1, 1, projection=wcs3d, slices=(1, "x", "y"))
    assert repr(ax.coords) == EXPECTED_REPR_1

    # Now slice in a way that all world coordinates are still present:

    fig.clear()
    ax = fig.add_subplot(1, 1, 1, projection=wcs3d, slices=("x", "y", 1))
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
    fig = Figure()
    fig.add_subplot(projection=time_spectral_wcs_2d)


@pytest.mark.skipif(TEX_UNAVAILABLE, reason="TeX is unavailable")
def test_simplify_labels_usetex(ignore_matplotlibrc, tmp_path):
    """Regression test for https://github.com/astropy/astropy/issues/8004."""
    mpl.rc("text", usetex=True)

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
    fig = Figure()
    ax = fig.add_subplot(frame_class=EllipticalFrame, projection=wcs)
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

    ticklabels = TickLabels(frame=MagicMock())
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
    fig = Figure()
    ax = fig.add_subplot(frame_class=frame_class, projection=wcs)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])

    assert ax.get_xlabel() == labels[0]
    assert ax.get_ylabel() == labels[1]
    for i in range(2):
        assert ax.coords[i].get_axislabel() == labels[i]


@pytest.mark.parametrize("atol", [0.2, 1.0e-8])
def test_bbox_size(atol):
    # Test for the size of a WCSAxes bbox (only have Matplotlib >= 3.0 now)
    extents = [11.38888888888889, 3.388888888888889, 576.0, 432.0]

    fig = Figure()
    canvas = FigureCanvasAgg(fig)
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
    canvas.figure.add_axes(ax)
    canvas.draw()
    renderer = fig.canvas.renderer
    ax_bbox = ax.get_tightbbox(renderer)

    # Enforce strict test only with reference Freetype version
    if atol < 0.1 and not FREETYPE_261:
        pytest.xfail(
            "Exact BoundingBox dimensions are only ensured with FreeType 2.6.1"
        )
    assert np.allclose(ax_bbox.extents, extents, atol=atol)


def test_wcs_type_transform_regression():
    fig = Figure()
    wcs = WCS(TARGET_HEADER)
    sliced_wcs = SlicedLowLevelWCS(wcs, np.s_[1:-1, 1:-1])
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    ax.get_transform(sliced_wcs)

    high_wcs = HighLevelWCSWrapper(sliced_wcs)
    ax.get_transform(sliced_wcs)


@pytest.mark.parametrize(("cunit", "cdelt"), [("deg", 10 / 3600), ("arcsec", 10)])
def test_get_transform_frame_with_non_degree_world_units(cunit, cdelt):
    # A transform from a coordinate frame must agree with the WCS it was built
    # from, despite different units in the WCS.
    wcs = WCS(naxis=2, preserve_units=True)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = [cunit, cunit]
    wcs.wcs.crpix = [50.5, 50.5]
    wcs.wcs.crval = [0, 0]
    wcs.wcs.cdelt = [cdelt, cdelt]
    wcs.wcs.set()

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # get_transform ends in display coordinates, removing transData leaves pixels
    world_to_pixel = ax.get_transform("icrs") - ax.transData

    coord = SkyCoord(0.05, 0, unit="deg")
    assert_allclose(
        world_to_pixel.transform(np.array([[0.05, 0.0]])),
        [wcs.world_to_pixel(coord)],
    )


def test_multiple_draws_grid_contours(tmp_path):
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS())
    ax.grid(color="black", grid_type="contours")
    fig.savefig(tmp_path / "plot.png")
    fig.savefig(tmp_path / "plot.png")


def test_get_coord_range_nan_regression():
    # Test to make sure there is no internal casting of NaN to integers
    # NumPy 1.24 raises a RuntimeWarning if a NaN is cast to an integer

    fig = Figure()
    wcs = WCS(TARGET_HEADER)
    wcs.wcs.crval[0] = 0  # Re-position the longitude wrap to the middle
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

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
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS())
    with pytest.raises(ValueError, match="Cannot use images with origin='upper"):
        ax.imshow(np.ones(100).reshape(10, 10), origin="upper")


def test_label_setting():
    fig = Figure()
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
    fig = Figure()
    _canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(1, 1, 1, projection=WCS())

    assert ax.get_tightbbox(fig.canvas.get_renderer()) is not None
    ax.set_visible(False)
    assert ax.get_tightbbox(fig.canvas.get_renderer()) is None


def test_tightbbox_pgf_measurement_renderer(tmp_path, ignore_matplotlibrc):
    # Regression test for a bug that caused get_tightbbox to emit actual
    # drawing commands, which fails with renderers that can only be used
    # for measurement, such as the PGF layout renderer (which has no
    # output stream), e.g. when saving a figure with layout='constrained'
    # to a PGF file.
    from matplotlib.backends.backend_pgf import FigureCanvasPgf

    fig = Figure(layout="constrained")
    FigureCanvasPgf(fig)
    ax = fig.add_subplot(projection=WCS())

    # Hide all text so that no TeX installation is needed
    for coord in ax.coords:
        coord.set_ticklabel_visible(False)
        coord.set_axislabel("")

    fig.savefig(tmp_path / "test.pgf")

    # When saving the plot above, get_tightbox should get called, but we can
    # double check that it also passes without failing if called explicitly
    assert ax.get_tightbbox(fig.canvas.get_renderer()) is not None


def test_get_axislabel_default():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = "RA---TAN", "DEC--TAN"

    fig = Figure()
    _canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    assert ax.coords[0].get_axislabel() == "pos.eq.ra"
    assert ax.coords[1].get_axislabel() == "pos.eq.dec"

    # Make sure that setting axis labels explicitly works, including to an
    # empty string
    ax.coords[0].set_axislabel("Right Ascension")
    ax.coords[1].set_axislabel("")
    assert ax.coords[0].get_axislabel() == "Right Ascension"
    assert ax.coords[1].get_axislabel() == ""


def test_plot_coord_slicing(ignore_matplotlibrc):
    # Regression test to ensure plot_coord does not raise a NotImplementedError
    # after slicing coordinates (issue #10254).

    cube_header = get_pkg_data_filename("data/cube_header")
    cube_header = fits.Header.fromtextfile(cube_header)

    fig = Figure(figsize=(6, 6))
    ax = fig.add_subplot(projection=WCS(cube_header), slices=("x", "y", 0))

    c = SkyCoord(52 * u.deg, 30.5 * u.deg)
    ax.plot_coord(c, "o")


SIMPLIFY_CASES = [
    (["13d14m15s", "13d14m15s"], ["13d14m15s", "15s"]),
    (["13d14m5s", "13d14m15s"], ["13d14m5s", "15s"]),
    (
        ["−29°25'55\"", "−29°25'54\"", "−29°25'53\"", "−29°25'52\"", "−29°25'51\""],
        ["−29°25'55\"", '54"', '53"', '52"', '51"'],
    ),
    (
        [
            r"0$\mathregular{^h}$",
            r"10$\mathregular{^h}$",
            r"10$\mathregular{^h}$",
            r"5$\mathregular{^h}$",
        ],
        [
            r"0$\mathregular{^h}$",
            r"10$\mathregular{^h}$",
            r"10$\mathregular{^h}$",
            r"5$\mathregular{^h}$",
        ],
    ),
    (
        [
            r"0$\mathregular{^h}$10$\mathregular{^m}$5$\mathregular{^m}$",
            r"0$\mathregular{^h}$10$\mathregular{^m}$10$\mathregular{^m}$",
            r"0$\mathregular{^h}$10$\mathregular{^m}$15$\mathregular{^m}$",
            r"0$\mathregular{^h}$10$\mathregular{^m}$15$\mathregular{^m}$",
        ],
        [
            r"0$\mathregular{^h}$10$\mathregular{^m}$5$\mathregular{^m}$",
            r"10$\mathregular{^m}$",
            r"15$\mathregular{^m}$",
            r"15$\mathregular{^m}$",
        ],
    ),
    (
        [
            "$17^{\\mathrm{h}}47^{\\mathrm{m}}53.8^{\\mathrm{s}}$",
            "$17^{\\mathrm{h}}47^{\\mathrm{m}}53.6^{\\mathrm{s}}$",
        ],
        [
            "$17^{\\mathrm{h}}47^{\\mathrm{m}}53.8^{\\mathrm{s}}$",
            "$53.6^{\\mathrm{s}}$",
        ],
    ),
]


@pytest.mark.parametrize(("before", "after"), SIMPLIFY_CASES)
def test_simplify_cases(before, after):
    # Regression test for a bug that caused the simplification to not work
    # correctly in the presence of dollar signs in LaTeX strings.

    ticklabels = TickLabels(frame=MagicMock())
    expected_labels = []

    for i, label in enumerate(before):
        ticklabels.add(
            axis="axis",
            world=0,
            angle=0,
            text=label,
            axis_displacement=i,
            data=(i, i),
        )

    ticklabels.simplify_labels()
    assert ticklabels.text["axis"] == after


def test_simplify_labels_prefix(ignore_matplotlibrc):
    # Regression test for crashes in simplify_labels (IndexError or
    # UnboundLocalError/NameError) when one tick label is a prefix of an
    # adjacent one, which can happen with tick labels of varying length.
    ticklabels = TickLabels(frame=MagicMock())
    for i, label in enumerate(["10d30m", "10d3"]):
        ticklabels.add(
            axis="axis", world=0, angle=0, text=label, axis_displacement=i, data=(i, i)
        )
    ticklabels.simplify_labels()
    assert ticklabels.text["axis"] == ["10d30m", "10d3"]


def test_get_transform_unit_mismatch():
    """
    Regression test for a bug that caused get_transform to ignore differences
    in WCS units.

    https://github.com/astropy/astropy/issues/18246
    """

    wcs_deg = WCS(naxis=2, preserve_units=True)
    wcs_deg.wcs.ctype = "RA---TAN", "DEC--TAN"
    wcs_deg.wcs.crval = 20, 30
    wcs_deg.wcs.cunit = "deg", "deg"
    wcs_deg.wcs.crpix = 1, 1
    wcs_deg.wcs.cdelt = 1 / 60, 1 / 60

    wcs_arcmin = WCS(naxis=2, preserve_units=True)
    wcs_arcmin.wcs.ctype = "RA---TAN", "DEC--TAN"
    wcs_arcmin.wcs.crval = 1200, 1800
    wcs_arcmin.wcs.cunit = "arcmin", "arcmin"
    wcs_arcmin.wcs.crpix = 1, 1
    wcs_arcmin.wcs.cdelt = 1, 1

    fig = Figure(figsize=(6, 6))
    ax = fig.add_subplot(projection=wcs_arcmin)
    transform1 = ax.get_transform(wcs_arcmin)
    transform2 = ax.get_transform(wcs_deg)

    # Since the two WCS are equivalent, the returned transforms should also
    # be equivalent

    rng = np.random.default_rng(12345)
    pixels = rng.uniform(0, 100, (2, 100))

    assert_allclose(transform1.transform(pixels), transform2.transform(pixels))


def test_get_coords_overlay_wcs_honors_coord_meta(ignore_matplotlibrc):
    # Regression test for a bug where get_coords_overlay silently discarded
    # a caller-supplied coord_meta when the frame passed in was a WCS,
    # unlike every other frame type.
    wcs = WCS(TARGET_HEADER)

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    custom_coord_meta = {
        "type": ("scalar", "scalar"),
        "wrap": (None, None),
        "unit": (u.deg, u.deg),
        "name": ("x", "y"),
    }

    overlay = ax.get_coords_overlay(WCS(TARGET_HEADER), coord_meta=custom_coord_meta)

    assert overlay[0].coord_type == "scalar"
    assert overlay[1].coord_type == "scalar"


def test_get_coords_overlay_elliptical_frame(ignore_matplotlibrc):
    """
    Regression test for calling get_coords_overlay on axes with a
    non-rectangular frame (e.g., EllipticalFrame).

    Previously this would emit an AstropyDeprecationWarning because 't'
    and 'r' are not valid spine names for EllipticalFrame.
    """
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["GLON-AIT", "GLAT-AIT"]
    wcs.wcs.crpix = [180.5, 90.5]
    wcs.wcs.cdelt = [-1, 1]
    wcs.wcs.crval = [0, 0]

    fig = Figure()
    ax = fig.add_subplot(projection=wcs, frame_class=EllipticalFrame)

    # This should not raise an AstropyDeprecationWarning
    overlay = ax.get_coords_overlay("icrs")
    assert overlay is not None


def test_reset_wcs_resets_display_coords_index(ignore_matplotlibrc):
    # Regression test for a bug where reset_wcs dropped overlays from
    # _all_coords without resetting _display_coords_index, so that if the
    # cursor display had been advanced to an overlay (by pressing "w"),
    # format_coord would later raise an IndexError.
    wcs = WCS(TARGET_HEADER)
    fig = Figure()
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    canvas.draw()

    ax.get_coords_overlay("icrs")

    # Simulate the user pressing "w" to advance the cursor display to the
    # overlay coordinates.
    event = KeyEvent("key_press_event", canvas, "w")
    canvas.callbacks.process("key_press_event", event)

    ax.reset_wcs(wcs=wcs)

    # This should not raise an IndexError.
    assert ax.format_coord(10, 10) != ""

    # Pressing "w" once more cycles past the world coordinates to the pixel
    # display mode, which is still valid after reset_wcs and should survive it.
    canvas.callbacks.process("key_press_event", event)
    assert "pixel" in ax.format_coord(10, 10)

    ax.reset_wcs(wcs=wcs)
    assert "pixel" in ax.format_coord(10, 10)


def test_auto_assign_coord_positions_no_consistent_option(
    ignore_matplotlibrc, tmp_path
):
    # Regression test for a bug where auto_assign_coord_positions crashed
    # with a TypeError if no assignment of spines to coordinates was
    # consistent with fixed tick positions (best_option stayed None).
    wcs = WCS(TARGET_HEADER)
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    # Coordinate 0 has fixed tick label position on 'b', so 'b' is excluded
    # from the spines available for automatic placement.
    ax.coords[0].set_ticklabel_position("b")
    ax.coords[0].set_ticks_position("b")

    # Coordinate 1 has fixed tick position on 'b' but automatic tick labels,
    # so every candidate spine assignment is inconsistent.
    ax.coords[1].set_ticks_position("b")
    ax.coords[1].set_ticklabel_position("#")

    fig.savefig(tmp_path / "plot.png")


@pytest.mark.parametrize("cunit", ["deg", "arcsec"])
def test_beam_scalebar_non_degree_world_units(cunit):
    # Regression test for a bug that caused beams and scale bars to have the
    # wrong size for WCS with world units other than degrees.
    wcs = WCS(naxis=2, preserve_units=True)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = [cunit, cunit]
    wcs.wcs.crval = [(10 * u.deg).to_value(cunit), (20 * u.deg).to_value(cunit)]
    wcs.wcs.cdelt = [(-0.01 * u.deg).to_value(cunit), (0.01 * u.deg).to_value(cunit)]

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    add_beam(ax, major=0.2 * u.deg, minor=0.1 * u.deg, angle=0 * u.deg)
    (ellipse,) = ax.artists[-1].get_child().get_children()
    assert_allclose([ellipse.width, ellipse.height], [10, 20])

    add_scalebar(ax, 0.3 * u.deg)
    (rectangle,) = ax.artists[-1].size_bar.get_children()
    assert_allclose(rectangle.get_width(), 30)


def test_beam_scalebar_mixed_world_units():
    # The two celestial axes can have different units, in which case the pixel
    # scale has to be converted to degrees before combining the world axes.
    wcs = WCS(naxis=2, preserve_units=True)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = ["deg", "arcsec"]
    wcs.wcs.crval = [10, (20 * u.deg).to_value(u.arcsec)]
    wcs.wcs.cdelt = [-0.01, (0.01 * u.deg).to_value(u.arcsec)]
    angle = np.radians(30)
    wcs.wcs.pc = [
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)],
    ]

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    add_beam(ax, major=0.2 * u.deg, minor=0.1 * u.deg, angle=0 * u.deg)
    (ellipse,) = ax.artists[-1].get_child().get_children()
    assert_allclose([ellipse.width, ellipse.height], [10, 20])

    add_scalebar(ax, 0.3 * u.deg)
    (rectangle,) = ax.artists[-1].size_bar.get_children()
    assert_allclose(rectangle.get_width(), 30)


def _celestial_wcs(cunit):
    wcs = WCS(naxis=2, preserve_units=True)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = [cunit, cunit]
    wcs.wcs.crpix = [50.5, 50.5]
    wcs.wcs.crval = [(10 * u.deg).to_value(cunit), (20 * u.deg).to_value(cunit)]
    wcs.wcs.cdelt = [(-0.01 * u.deg).to_value(cunit), (0.01 * u.deg).to_value(cunit)]
    wcs.wcs.set()
    return wcs


@pytest.fixture
def celestial_axes():
    fig = Figure()
    FigureCanvasAgg(fig)
    return fig.add_subplot(1, 1, 1, projection=_celestial_wcs("arcsec"))


@pytest.mark.parametrize("cunit", ["deg", "arcsec"])
@pytest.mark.parametrize("frame", ["world", "icrs", "galactic"])
@pytest.mark.parametrize("unit", [u.deg, u.arcsec])
def test_plot_quantity_with_transform(cunit, frame, unit):
    # Regression test for a bug that caused the unit of a Quantity plotted
    # with a world or frame transform to be ignored.
    # https://github.com/astropy/astropy/issues/20777
    wcs = _celestial_wcs(cunit)
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    transform = ax.get_transform(frame)
    world_to_pixel = transform - ax.transData

    coord = SkyCoord(10.1, 20.1, unit="deg")
    sph = coord.transform_to("icrs" if frame == "world" else frame).spherical
    lon, lat = sph.lon.to(unit), sph.lat.to(unit)

    expected = [wcs.world_to_pixel(coord)]
    (line,) = ax.plot(lon, lat, "o", transform=transform)
    assert_allclose(world_to_pixel.transform(line.get_xydata()), expected)
    points = ax.scatter(lon, lat, transform=transform)
    assert_allclose(world_to_pixel.transform(points.get_offsets()), expected)
    text = ax.text(lon, lat, "label", transform=transform)
    position = np.array([text.get_position()])
    assert_allclose(world_to_pixel.transform(position), expected)
    # Lists of Quantity scalars are also converted
    (line,) = ax.plot([lon], [lat], "o", transform=transform)
    assert_allclose(world_to_pixel.transform(line.get_xydata()), expected)
    points = ax.scatter([lon], [lat], transform=transform)
    assert_allclose(world_to_pixel.transform(points.get_offsets()), expected)


_X = [10.1, 10.2, 10.3] * u.deg
_Y = [20.1, 20.2, 20.4] * u.deg
_X2D, _Y2D = np.meshgrid(_X, _Y)
_Z2D = np.arange(9).reshape(3, 3)

# Plotting methods that accept a transform, called with the functions x and y
# applied to the coordinates
QUANTITY_PLOT_CALLS = {
    "plot": lambda ax, x, y, **kwargs: ax.plot(x(_X), y(_Y), **kwargs),
    "plot_scalars": lambda ax, x, y, **kwargs: ax.plot(
        x(_X[0]), y(_Y[0]), "o", **kwargs
    ),
    "plot_groups": lambda ax, x, y, **kwargs: ax.plot(
        x(_X), y(_Y), "o", x(_X), y(_Y[::-1]), x(_X[:2]), y(_Y[:2]), "-", **kwargs
    ),
    "step": lambda ax, x, y, **kwargs: ax.step(x(_X), y(_Y), **kwargs),
    "fill": lambda ax, x, y, **kwargs: ax.fill(x(_X), y(_Y), "r", **kwargs),
    "scatter": lambda ax, x, y, **kwargs: ax.scatter(x(_X), y(_Y), **kwargs),
    "scatter_keywords": lambda ax, x, y, **kwargs: ax.scatter(
        x=x(_X), y=y(_Y), **kwargs
    ),
    "text": lambda ax, x, y, **kwargs: ax.text(x(_X[0]), y(_Y[0]), "label", **kwargs),
    "errorbar": lambda ax, x, y, **kwargs: ax.errorbar(
        x(_X), y(_Y), xerr=x(_X / 500), yerr=y(_Y / 500), capsize=3, **kwargs
    ),
    "pcolor": lambda ax, x, y, **kwargs: ax.pcolor(x(_X2D), y(_Y2D), _Z2D, **kwargs),
    "pcolormesh": lambda ax, x, y, **kwargs: ax.pcolormesh(
        x(_X2D), y(_Y2D), _Z2D, **kwargs
    ),
    "contour": lambda ax, x, y, **kwargs: ax.contour(x(_X2D), y(_Y2D), _Z2D, **kwargs),
    "contourf": lambda ax, x, y, **kwargs: ax.contourf(
        x(_X2D), y(_Y2D), _Z2D, **kwargs
    ),
}


def _render(call, *args, frame=None):
    fig = Figure(figsize=(2, 2), dpi=50)
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(1, 1, 1, projection=_celestial_wcs("arcsec"))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    kwargs = {} if frame is None else {"transform": ax.get_transform(frame)}
    call(ax, *args, **kwargs)
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba()).copy()


@pytest.mark.parametrize("call", QUANTITY_PLOT_CALLS.values(), ids=QUANTITY_PLOT_CALLS)
@pytest.mark.parametrize(("frame", "unit"), [("world", u.arcsec), ("icrs", u.deg)])
def test_quantity_same_as_values_in_transform_units(call, frame, unit):
    # Quantities have to give the same result as plain values in the units
    # that the transform expects, which for an arcsecond WCS are arcseconds
    # for the world transform and degrees for a celestial frame transform.
    def to_arcmin(q):
        return q.to(u.arcmin)

    def to_value(q):
        return q.to_value(unit)

    expected = _render(call, to_value, to_value, frame=frame)
    assert np.array_equal(_render(call, to_arcmin, to_arcmin, frame=frame), expected)
    # Including when only one of the coordinates is a Quantity
    assert np.array_equal(_render(call, to_arcmin, to_value, frame=frame), expected)
    # Make sure that the comparison is sensitive to the units
    assert not np.array_equal(
        _render(call, lambda q: q.value, lambda q: q.value, frame="world"),
        _render(call, lambda q: q.value, lambda q: q.value, frame="icrs"),
    )


@pytest.mark.parametrize("call", QUANTITY_PLOT_CALLS.values(), ids=QUANTITY_PLOT_CALLS)
def test_quantity_without_transform(celestial_axes, call):
    # Without a transform, or with the pixel or data transforms, coordinates
    # are pixel coordinates, so quantities have to be in pixels
    for kwargs in [
        {},
        {"transform": celestial_axes.transData},
        {"transform": celestial_axes.get_transform("pixel")},
        {"transform": celestial_axes.get_transform(_celestial_wcs("deg"))},
    ]:
        with pytest.raises(u.UnitConversionError, match="pixel coordinates"):
            call(celestial_axes, lambda q: q, lambda q: q, **kwargs)
        call(celestial_axes, lambda q: q.value * u.pix, lambda q: q.value, **kwargs)


def test_pixel_quantity_values(celestial_axes):
    ax = celestial_axes
    (line,) = ax.plot([1, 2] * u.pix, [3, 4] * u.pix)
    assert_allclose(line.get_xydata(), [[1, 3], [2, 4]])
    # Including with prefixed pixel units, and for a single y argument
    (line,) = ax.plot([0.001, 0.002] * u.kpix, [3, 4])
    assert_allclose(line.get_xydata(), [[1, 3], [2, 4]])
    (line,) = ax.plot([3, 4] * u.pix, "o")
    assert_allclose(line.get_xydata(), [[0, 3], [1, 4]])
    points = ax.scatter([1, 2] * u.pix, [3, 4] * u.pix)
    assert_allclose(points.get_offsets(), [[1, 3], [2, 4]])


@pytest.mark.parametrize("unit", [u.kpix, u.mpix, u.Unit(2 * u.pix)])
@pytest.mark.parametrize("use_transform", [False, True])
def test_pixel_derived_quantity_values(celestial_axes, unit, use_transform):
    # Units that are derived from pixels are converted to pixels
    ax = celestial_axes
    kwargs = {"transform": ax.get_transform("pixel")} if use_transform else {}
    x = ([1, 2] * u.pix).to(unit)
    y = ([3, 4] * u.pix).to(unit)
    (line,) = ax.plot(x, y, **kwargs)
    assert_allclose(line.get_xydata(), [[1, 3], [2, 4]])
    points = ax.scatter(x, y, **kwargs)
    assert_allclose(points.get_offsets(), [[1, 3], [2, 4]])
    text = ax.text(x[0], y[0], "label", **kwargs)
    assert_allclose(text.get_position(), [1, 3])


def test_plot_quantity_incompatible_unit(celestial_axes):
    ax = celestial_axes
    with pytest.raises(u.UnitConversionError, match="expects units of 'arcsec'"):
        ax.plot(1 * u.m, 1 * u.deg, transform=ax.get_transform("world"))
    with pytest.raises(u.UnitConversionError, match="expects units of 'deg'"):
        ax.scatter(1 * u.deg, 1 * u.pix, transform=ax.get_transform("icrs"))


def test_plot_quantity_unknown_transform(celestial_axes):
    # The units that arbitrary transforms expect are not known
    ax = celestial_axes
    for transform in [
        ax.transAxes,
        Affine2D().translate(1, 1) + ax.get_transform("world"),
        ax.get_transform(Affine2D()),
    ]:
        with pytest.raises(TypeError, match="units that other transforms expect"):
            ax.plot([0.1, 0.2] * u.deg, [0.1, 0.2] * u.deg, transform=transform)
        ax.plot([0.1, 0.2], [0.1, 0.2], transform=transform)


def test_plot_quantity_axes_without_wcs():
    # The world units are not known for axes created with a transform
    fig = Figure()
    coord_meta = {
        "type": ("longitude", "latitude"),
        "wrap": (None, None),
        "unit": (u.deg, u.deg),
        "name": ("lon", "lat"),
    }
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], transform=Affine2D(), coord_meta=coord_meta)
    fig.add_axes(ax)
    transform = ax.get_transform("world")
    with pytest.raises(TypeError, match="units that other transforms expect"):
        ax.plot([1, 2] * u.deg, [1, 2] * u.deg, transform=transform)
    ax.plot([1, 2], [1, 2], transform=transform)
    # Pixel coordinates are still known
    (line,) = ax.plot([1, 2] * u.pix, [3, 4] * u.pix)
    assert_allclose(line.get_xydata(), [[1, 3], [2, 4]])


def test_plot_quantity_with_quantity_support(celestial_axes):
    # Quantities are converted before Matplotlib sees them, so the result does
    # not depend on whether quantity_support is enabled
    ax = celestial_axes
    transform = ax.get_transform("world")
    with quantity_support():
        (line,) = ax.plot(10.1 * u.deg, 20.1 * u.deg, "o", transform=transform)
        with pytest.raises(u.UnitConversionError, match="pixel coordinates"):
            ax.plot(10.1 * u.deg, 20.1 * u.deg, "o")
    assert_allclose(line.get_xydata(), [[36360, 72360]])
    assert ax.xaxis.get_units() is None


def test_plot_quantity_1d():
    # For 1D frames, only the x axis is a world coordinate
    wcs = WCS(naxis=1)
    wcs.wcs.ctype = ["WAVE"]
    wcs.wcs.cunit = ["um"]
    wcs.wcs.crval = [1]
    wcs.wcs.crpix = [1]
    wcs.wcs.cdelt = [0.01]
    wcs.wcs.set()

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    transform = ax.get_transform("world")
    (line,) = ax.plot(15000 * u.AA, 3 * u.Jy, "o", transform=transform)
    assert_allclose((transform - ax.transData).transform(line.get_xydata()), [[50, 3]])
    # Without a transform, x is in pixels and y is not a coordinate
    (line,) = ax.plot([1, 2] * u.pix, [3, 4] * u.Jy)
    assert_allclose(line.get_xydata(), [[1, 3], [2, 4]])
    (line,) = ax.plot([3, 4] * u.Jy)
    assert_allclose(line.get_xydata(), [[0, 3], [1, 4]])
    with pytest.raises(u.UnitConversionError, match="pixel coordinates"):
        ax.plot([1, 2] * u.um, [3, 4] * u.Jy)


def test_plot_quantity_1d_slice_of_2d():
    # For a 1D slice of a celestial WCS, the y axis also shows data values
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=WCS(TARGET_HEADER), slices=("x", 0))
    (line,) = ax.plot([1, 2] * u.pix, [3, 4] * u.Jy)
    assert_allclose(line.get_xydata(), [[1, 3], [2, 4]])


def test_reset_wcs_keeps_quantity_positions():
    # Quantities are converted when plotted, to the units of the WCS that the
    # transform uses, so resetting the WCS to one with different units does
    # not move existing artists.
    def get_pixel_position(line):
        world_to_pixel = line.get_transform() - line.axes.transData
        return world_to_pixel.transform(line.get_xydata())

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1, projection=_celestial_wcs("deg"))
    (line,) = ax.plot(10.1 * u.deg, 20.2 * u.deg, transform=ax.get_transform("world"))
    expected = get_pixel_position(line)

    ax.reset_wcs(_celestial_wcs("arcsec"))
    assert_allclose(get_pixel_position(line), expected)

    # New artists use the units of the new WCS
    (line,) = ax.plot(10.1 * u.deg, 20.2 * u.deg, transform=ax.get_transform("world"))
    assert_allclose(line.get_xydata(), [[36360, 72720]])
    assert_allclose(get_pixel_position(line), expected, rtol=1e-5)

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import matplotlib.lines
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib import rc_context
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Rectangle

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.tests.figures import figure_test
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyUserWarning
from astropy.visualization.wcsaxes import WCSAxes, add_beam, add_scalebar
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.visualization.wcsaxes.patches import Quadrangle, SphericalCircle
from astropy.wcs import WCS


class BaseImageTests:
    @classmethod
    def setup_class(cls):
        msx_header = get_pkg_data_filename("data/msx_header")
        cls.msx_header = fits.Header.fromtextfile(msx_header)

        rosat_header = get_pkg_data_filename("data/rosat_header")
        cls.rosat_header = fits.Header.fromtextfile(rosat_header)

        twoMASS_k_header = get_pkg_data_filename("data/2MASS_k_header")
        cls.twoMASS_k_header = fits.Header.fromtextfile(twoMASS_k_header)

        cube_header = get_pkg_data_filename("data/cube_header")
        cls.cube_header = fits.Header.fromtextfile(cube_header)

        slice_header = get_pkg_data_filename("data/slice_header")
        cls.slice_header = fits.Header.fromtextfile(slice_header)

    def teardown_method(self, method):
        plt.close("all")


class TestBasic(BaseImageTests):
    @figure_test
    def test_tight_layout(self):
        # Check that tight_layout works on a WCSAxes.
        fig = plt.figure(figsize=(8, 6))
        for i in (1, 2):
            fig.add_subplot(2, 1, i, projection=WCS(self.msx_header))
        fig.tight_layout()
        return fig

    @figure_test
    def test_image_plot(self):
        # Test for plotting image and also setting values of ticks
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8], projection=WCS(self.msx_header), aspect="equal"
        )
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        ax.coords[0].set_ticks([-0.30, 0.0, 0.20] * u.degree, size=5, width=1)
        return fig

    @figure_test
    def test_axes_off(self):
        # Test for turning the axes off
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=WCS(self.msx_header))
        ax.imshow(np.arange(12).reshape((3, 4)))
        ax.set_axis_off()
        return fig

    @figure_test
    @pytest.mark.parametrize("axisbelow", [True, False, "line"])
    def test_axisbelow(self, axisbelow):
        # Test that tick marks, labels, and gridlines are drawn with the
        # correct zorder controlled by the axisbelow property.
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8], projection=WCS(self.msx_header), aspect="equal"
        )
        ax.set_axisbelow(axisbelow)
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        ax.coords[0].set_ticks([-0.30, 0.0, 0.20] * u.degree, size=5, width=1)
        ax.grid()
        ax.coords[0].set_auto_axislabel(False)
        ax.coords[1].set_auto_axislabel(False)

        # Add an image (default zorder=0).
        ax.imshow(np.zeros((64, 64)))

        # Add a patch (default zorder=1).
        r = Rectangle((30.0, 50.0), 60.0, 50.0, facecolor="green", edgecolor="red")
        ax.add_patch(r)

        # Add a line (default zorder=2).
        ax.plot([32, 128], [32, 128], linewidth=10)

        return fig

    @figure_test
    def test_contour_overlay(self):
        # Test for overlaying contours on images
        path = get_pkg_data_filename("galactic_center/gc_msx_e.fits")
        with fits.open(path) as pf:
            data = pf[0].data

        wcs_msx = WCS(self.msx_header)

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.15, 0.15, 0.8, 0.8],
            projection=WCS(self.twoMASS_k_header),
            aspect="equal",
        )
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        # Overplot contour
        ax.contour(
            data,
            transform=ax.get_transform(wcs_msx),
            colors="orange",
            levels=[2.5e-5, 5e-5, 1.0e-4],
        )
        ax.coords[0].set_ticks(size=5, width=1)
        ax.coords[1].set_ticks(size=5, width=1)
        ax.set_xlim(0.0, 720.0)
        ax.set_ylim(0.0, 720.0)

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @figure_test
    def test_contourf_overlay(self):
        # Test for overlaying contours on images
        path = get_pkg_data_filename("galactic_center/gc_msx_e.fits")
        with fits.open(path) as pf:
            data = pf[0].data

        wcs_msx = WCS(self.msx_header)

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.15, 0.15, 0.8, 0.8],
            projection=WCS(self.twoMASS_k_header),
            aspect="equal",
        )
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        # Overplot contour
        ax.contourf(
            data, transform=ax.get_transform(wcs_msx), levels=[2.5e-5, 5e-5, 1.0e-4]
        )
        ax.coords[0].set_ticks(size=5, width=1)
        ax.coords[1].set_ticks(size=5, width=1)
        ax.set_xlim(0.0, 720.0)
        ax.set_ylim(0.0, 720.0)

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @figure_test
    def test_overlay_features_image(self):
        # Test for overlaying grid, changing format of ticks, setting spacing
        # and number of ticks

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.25, 0.25, 0.65, 0.65], projection=WCS(self.msx_header), aspect="equal"
        )

        # Change the format of the ticks
        ax.coords[0].set_major_formatter("dd:mm:ss")
        ax.coords[1].set_major_formatter("dd:mm:ss.ssss")

        # Overlay grid on image
        ax.grid(color="red", alpha=1.0, lw=1, linestyle="dashed")

        # Set the spacing of ticks on the 'glon' axis to 4 arcsec
        ax.coords["glon"].set_ticks(spacing=4 * u.arcsec, size=5, width=1)

        # Set the number of ticks on the 'glat' axis to 9
        ax.coords["glat"].set_ticks(number=9, size=5, width=1)

        # Set labels on axes
        ax.coords["glon"].set_axislabel("Galactic Longitude", minpad=1.6)
        ax.coords["glat"].set_axislabel("Galactic Latitude", minpad=-0.75)

        # Change the frame linewidth and color
        ax.coords.frame.set_color("red")
        ax.coords.frame.set_linewidth(2)

        assert ax.coords.frame.get_color() == "red"
        assert ax.coords.frame.get_linewidth() == 2

        return fig

    @figure_test
    def test_curvilinear_grid_patches_image(self):
        # Overlay curvilinear grid and patches on image

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8], projection=WCS(self.rosat_header), aspect="equal"
        )

        ax.set_xlim(-0.5, 479.5)
        ax.set_ylim(-0.5, 239.5)

        ax.grid(color="black", alpha=1.0, lw=1, linestyle="dashed")

        p = Circle((300, 100), radius=40, ec="yellow", fc="none")
        ax.add_patch(p)

        p = Circle(
            (30.0, 20.0),
            radius=20.0,
            ec="orange",
            fc="none",
            transform=ax.get_transform("world"),
        )
        ax.add_patch(p)

        p = Circle(
            (60.0, 50.0),
            radius=20.0,
            ec="red",
            fc="none",
            transform=ax.get_transform("fk5"),
        )
        ax.add_patch(p)

        p = Circle(
            (40.0, 60.0),
            radius=20.0,
            ec="green",
            fc="none",
            transform=ax.get_transform("galactic"),
        )
        ax.add_patch(p)

        return fig

    @figure_test
    def test_cube_slice_image(self):
        # Test for cube slicing

        fig = plt.figure()
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8],
            projection=WCS(self.cube_header),
            slices=(50, "y", "x"),
            aspect="equal",
        )

        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)

        ax.coords[2].set_axislabel("Velocity m/s")

        ax.coords[1].set_ticks(spacing=0.2 * u.deg, width=1)
        ax.coords[2].set_ticks(spacing=400 * u.m / u.s, width=1)

        ax.coords[1].set_ticklabel(exclude_overlapping=True)
        ax.coords[2].set_ticklabel(exclude_overlapping=True)

        ax.coords[0].grid(grid_type="contours", color="purple", linestyle="solid")
        ax.coords[1].grid(grid_type="contours", color="orange", linestyle="solid")
        ax.coords[2].grid(grid_type="contours", color="red", linestyle="solid")

        return fig

    @figure_test
    def test_cube_slice_image_lonlat(self):
        # Test for cube slicing. Here we test with longitude and latitude since
        # there is some longitude-specific code in _update_grid_contour.

        fig = plt.figure()
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8],
            projection=WCS(self.cube_header),
            slices=("x", "y", 50),
            aspect="equal",
        )

        ax.set_xlim(-0.5, 106.5)
        ax.set_ylim(-0.5, 106.5)

        ax.coords[0].grid(grid_type="contours", color="blue", linestyle="solid")
        ax.coords[1].grid(grid_type="contours", color="red", linestyle="solid")

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @figure_test
    def test_plot_coord(self):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.15, 0.15, 0.8, 0.8],
            projection=WCS(self.twoMASS_k_header),
            aspect="equal",
        )
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        c = SkyCoord(266 * u.deg, -29 * u.deg)
        lines = ax.plot_coord(c, "o")

        # Test that plot_coord returns the results from ax.plot
        assert isinstance(lines, list)
        assert isinstance(lines[0], matplotlib.lines.Line2D)

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @figure_test
    def test_scatter_coord(self):
        from matplotlib.collections import PathCollection

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.15, 0.15, 0.8, 0.8],
            projection=WCS(self.twoMASS_k_header),
            aspect="equal",
        )
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        c = SkyCoord(266 * u.deg, -29 * u.deg)
        sc = ax.scatter_coord(c, marker="o")

        # Test that plot_coord returns the results from ax.plot
        assert isinstance(sc, PathCollection)

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @figure_test
    def test_plot_line(self):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.15, 0.15, 0.8, 0.8],
            projection=WCS(self.twoMASS_k_header),
            aspect="equal",
        )
        ax.set_xlim(-0.5, 720.5)
        ax.set_ylim(-0.5, 720.5)

        c = SkyCoord([266, 266.8] * u.deg, [-29, -28.9] * u.deg)
        ax.plot_coord(c)

        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)

        return fig

    @figure_test
    def test_changed_axis_units(self):
        # Test to see if changing the units of axis works
        fig = plt.figure()
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8],
            projection=WCS(self.cube_header),
            slices=(50, "y", "x"),
            aspect="equal",
        )
        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)
        ax.coords[0].set_ticks_position("")
        ax.coords[0].set_ticklabel_position("")
        ax.coords[0].set_axislabel_position("")
        ax.coords[1].set_ticks_position("lr")
        ax.coords[1].set_ticklabel_position("l")
        ax.coords[1].set_axislabel_position("l")
        ax.coords[2].set_ticks_position("bt")
        ax.coords[2].set_ticklabel_position("b")
        ax.coords[2].set_axislabel_position("b")
        ax.coords[2].set_major_formatter("x.xx")
        ax.coords[2].set_format_unit(u.km / u.s)
        ax.coords[2].set_axislabel("Velocity km/s")
        ax.coords[1].set_ticks(width=1)
        ax.coords[2].set_ticks(width=1)
        ax.coords[1].set_ticklabel(exclude_overlapping=True)
        ax.coords[2].set_ticklabel(exclude_overlapping=True)

        return fig

    @figure_test
    def test_minor_ticks(self):
        # Test for drawing minor ticks
        fig = plt.figure()
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8],
            projection=WCS(self.cube_header),
            slices=(50, "y", "x"),
            aspect="equal",
        )
        ax.set_xlim(-0.5, 52.5)
        ax.set_ylim(-0.5, 106.5)
        ax.coords[0].set_ticks_position("")
        ax.coords[0].set_ticklabel_position("")
        ax.coords[0].set_axislabel_position("")
        ax.coords[1].set_ticks_position("lr")
        ax.coords[1].set_ticklabel_position("l")
        ax.coords[1].set_axislabel_position("l")
        ax.coords[2].set_ticks_position("bt")
        ax.coords[2].set_ticklabel_position("b")
        ax.coords[2].set_axislabel_position("b")
        ax.coords[2].set_ticklabel(exclude_overlapping=True)
        ax.coords[1].set_ticklabel(exclude_overlapping=True)
        ax.coords[2].display_minor_ticks(True)
        ax.coords[1].display_minor_ticks(True)
        ax.coords[2].set_minor_frequency(3)
        ax.coords[1].set_minor_frequency(10)

        return fig

    @figure_test
    def test_ticks_labels(self):
        fig = plt.figure(figsize=(6, 6))
        ax = WCSAxes(fig, [0.1, 0.1, 0.7, 0.7], wcs=None)
        fig.add_axes(ax)
        ax.set_xlim(-0.5, 2)
        ax.set_ylim(-0.5, 2)
        ax.coords[0].set_ticks(size=10, color="blue", alpha=0.2, width=1)
        ax.coords[1].set_ticks(size=20, color="red", alpha=0.9, width=1)
        ax.coords[0].set_ticks_position("all")
        ax.coords[1].set_ticks_position("all")
        ax.coords[0].set_axislabel("X-axis", size=20)
        ax.coords[1].set_axislabel(
            "Y-axis",
            color="green",
            size=25,
            weight="regular",
            style="normal",
            family="cmtt10",
        )
        ax.coords[0].set_axislabel_position("t")
        ax.coords[1].set_axislabel_position("r")
        ax.coords[0].set_ticklabel(
            color="purple",
            size=15,
            alpha=1,
            weight="light",
            style="normal",
            family="cmss10",
        )
        ax.coords[1].set_ticklabel(
            color="black", size=18, alpha=0.9, weight="bold", family="cmr10"
        )
        ax.coords[0].set_ticklabel_position("all")
        ax.coords[1].set_ticklabel_position("r")

        return fig

    @figure_test
    def test_no_ticks(self):
        # Check that setting no ticks works
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8], projection=WCS(self.msx_header), aspect="equal"
        )
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        ax.coords[0].set_ticks(number=0)
        ax.coords[0].grid(True)
        return fig

    @figure_test
    def test_rcparams(self):
        # Test custom rcParams

        with rc_context(
            {
                "axes.labelcolor": "purple",
                "axes.labelsize": 14,
                "axes.labelweight": "bold",
                "axes.linewidth": 3,
                "axes.facecolor": "0.5",
                "axes.edgecolor": "green",
                "xtick.color": "red",
                "xtick.labelsize": 8,
                "xtick.direction": "in",
                "xtick.minor.visible": True,
                "xtick.minor.size": 5,
                "xtick.major.size": 20,
                "xtick.major.width": 3,
                "xtick.major.pad": 10,
                "grid.color": "blue",
                "grid.linestyle": ":",
                "grid.linewidth": 1,
                "grid.alpha": 0.5,
            }
        ):
            fig = plt.figure(figsize=(6, 6))
            ax = WCSAxes(fig, [0.15, 0.1, 0.7, 0.7], wcs=None)
            fig.add_axes(ax)
            ax.set_xlim(-0.5, 2)
            ax.set_ylim(-0.5, 2)
            ax.grid()
            ax.set_xlabel("X label")
            ax.set_ylabel("Y label")
            ax.coords[0].set_ticklabel(exclude_overlapping=True)
            ax.coords[1].set_ticklabel(exclude_overlapping=True)
            return fig

    @figure_test
    def test_tick_angles(self):
        # Test that tick marks point in the correct direction, even when the
        # axes limits extend only over a few FITS pixels. Addresses #45, #46.
        w = WCS()
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.crval = [90, 70]
        w.wcs.cdelt = [16, 16]
        w.wcs.crpix = [1, 1]
        w.wcs.radesys = "ICRS"
        w.wcs.equinox = 2000.0
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=w)
        ax.set_xlim(1, -1)
        ax.set_ylim(-1, 1)
        ax.grid(color="gray", alpha=0.5, linestyle="solid")
        ax.coords["ra"].set_ticks(color="red", size=20)
        ax.coords["dec"].set_ticks(color="red", size=20)
        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)
        return fig

    @figure_test
    def test_tick_angles_non_square_axes(self):
        # Test that tick marks point in the correct direction, even when the
        # axes limits extend only over a few FITS pixels, and the axes are
        # non-square.
        w = WCS()
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.crval = [90, 70]
        w.wcs.cdelt = [16, 16]
        w.wcs.crpix = [1, 1]
        w.wcs.radesys = "ICRS"
        w.wcs.equinox = 2000.0
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=w)
        ax.set_xlim(1, -1)
        ax.set_ylim(-1, 1)
        ax.grid(color="gray", alpha=0.5, linestyle="solid")
        ax.coords["ra"].set_ticks(color="red", size=20)
        ax.coords["dec"].set_ticks(color="red", size=20)
        # In previous versions, all angle axes defaulted to being displayed in
        # degrees. We now automatically show RA axes in hour angle units, but
        # for backward-compatibility with previous reference images we
        # explicitly use degrees here.
        ax.coords[0].set_format_unit(u.degree)
        return fig

    @figure_test
    def test_set_coord_type(self):
        # Test for setting coord_type
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes(
            [0.2, 0.2, 0.6, 0.6], projection=WCS(self.msx_header), aspect="equal"
        )
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        ax.coords[0].set_coord_type("scalar")
        ax.coords[1].set_coord_type("scalar")
        ax.coords[0].set_major_formatter("x.xxx")
        ax.coords[1].set_major_formatter("x.xxx")
        ax.coords[0].set_ticklabel(exclude_overlapping=True)
        ax.coords[1].set_ticklabel(exclude_overlapping=True)
        return fig

    @figure_test
    def test_ticks_regression(self):
        # Regression test for a bug that caused ticks aligned exactly with a
        # sampled frame point to not appear. This also checks that tick labels
        # don't get added more than once, and that no error occurs when e.g.
        # the top part of the frame is all at the same coordinate as one of the
        # potential ticks (which causes the tick angle calculation to return
        # NaN).
        wcs = WCS(self.slice_header)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.25, 0.25, 0.5, 0.5], projection=wcs, aspect="auto")
        limits = wcs.wcs_world2pix([0, 0], [35e3, 80e3], 0)[1]
        ax.set_ylim(*limits)
        ax.coords[0].set_ticks(spacing=0.002 * u.deg)
        ax.coords[1].set_ticks(spacing=5 * u.km / u.s)
        ax.coords[0].set_ticklabel(alpha=0.5)  # to see multiple labels
        ax.coords[1].set_ticklabel(alpha=0.5)
        ax.coords[0].set_ticklabel_position("all")
        ax.coords[1].set_ticklabel_position("all")
        return fig

    @figure_test
    def test_axislabels_regression(self):
        # Regression test for a bug that meant that if tick labels were made
        # invisible with ``set_visible(False)``, they were still added to the
        # list of bounding boxes for tick labels, but with default values of 0
        # to 1, which caused issues.
        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.25, 0.25, 0.5, 0.5], projection=wcs, aspect="auto")
        ax.coords[0].set_axislabel("Label 1")
        ax.coords[1].set_axislabel("Label 2")
        ax.coords[1].set_axislabel_visibility_rule("always")
        ax.coords[1].ticklabels.set_visible(False)
        return fig

    @figure_test(savefig_kwargs={"bbox_inches": "tight"})
    def test_noncelestial_angular(self, tmp_path):
        # Regression test for a bug that meant that when passing a WCS that had
        # angular axes and using set_coord_type to set the coordinates to
        # longitude/latitude, but where the WCS wasn't recognized as celestial,
        # the WCS units are not converted to deg, so we can't assume that
        # transform will always return degrees.

        wcs = WCS(naxis=2)

        wcs.wcs.ctype = ["solar-x", "solar-y"]
        wcs.wcs.cunit = ["arcsec", "arcsec"]

        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(1, 1, 1, projection=wcs)

        ax.imshow(np.zeros([1024, 1024]), origin="lower")

        ax.coords[0].set_coord_type("longitude", coord_wrap=180 * u.deg)
        ax.coords[1].set_coord_type("latitude")

        ax.coords[0].set_major_formatter("s.s")
        ax.coords[1].set_major_formatter("s.s")

        ax.coords[0].set_format_unit(u.arcsec, show_decimal_unit=False)
        ax.coords[1].set_format_unit(u.arcsec, show_decimal_unit=False)

        ax.grid(color="white", ls="solid")

        # Force drawing (needed for format_coord)
        fig.savefig(tmp_path / "nothing")

        assert ax.format_coord(512, 512) == "513.0 513.0 (world)"

        return fig

    @figure_test
    def test_patches_distortion(self, tmp_path):
        # Check how patches get distorted (and make sure that scatter markers
        # and SphericalCircle don't)

        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.25, 0.25, 0.5, 0.5], projection=wcs, aspect="equal")

        # Pixel coordinates
        r = Rectangle((30.0, 50.0), 60.0, 50.0, edgecolor="green", facecolor="none")
        ax.add_patch(r)

        # FK5 coordinates
        r = Rectangle(
            (266.4, -28.9),
            0.3,
            0.3,
            edgecolor="cyan",
            facecolor="none",
            transform=ax.get_transform("fk5"),
        )
        ax.add_patch(r)

        # FK5 coordinates
        c = Circle(
            (266.4, -29.1),
            0.15,
            edgecolor="magenta",
            facecolor="none",
            transform=ax.get_transform("fk5"),
        )
        ax.add_patch(c)

        # Pixel coordinates
        ax.scatter(
            [40, 100, 130],
            [30, 130, 60],
            s=100,
            edgecolor="red",
            facecolor=(1, 0, 0, 0.5),
        )

        # World coordinates (should not be distorted)
        ax.scatter(
            266.78238,
            -28.769255,
            transform=ax.get_transform("fk5"),
            s=300,
            edgecolor="red",
            facecolor="none",
        )

        # World coordinates (should not be distorted)
        r1 = SphericalCircle(
            (266.4 * u.deg, -29.1 * u.deg),
            0.15 * u.degree,
            edgecolor="purple",
            facecolor="none",
            transform=ax.get_transform("fk5"),
        )
        ax.add_patch(r1)

        r2 = SphericalCircle(
            SkyCoord(266.4 * u.deg, -29.1 * u.deg),
            0.15 * u.degree,
            edgecolor="purple",
            facecolor="none",
            transform=ax.get_transform("fk5"),
        )

        with pytest.warns(
            AstropyUserWarning,
            match="Received `center` of representation type "
            "<class 'astropy.coordinates.representation.CartesianRepresentation'> "
            "will be converted to SphericalRepresentation",
        ):
            r3 = SphericalCircle(
                SkyCoord(
                    x=-0.05486461,
                    y=-0.87204803,
                    z=-0.48633538,
                    representation_type="cartesian",
                ),
                0.15 * u.degree,
                edgecolor="purple",
                facecolor="none",
                transform=ax.get_transform("fk5"),
            )

        ax.coords[0].set_ticklabel_visible(False)
        ax.coords[1].set_ticklabel_visible(False)

        # Test to verify that SphericalCircle works irrespective of whether
        # the input(center) is a tuple or a SkyCoord object.
        assert (r1.get_xy() == r2.get_xy()).all()
        assert np.allclose(r1.get_xy(), r3.get_xy())
        assert np.allclose(r2.get_xy()[0], [266.4, -29.25])

        return fig

    @figure_test
    def test_quadrangle(self, tmp_path):
        # Test that Quadrangle can have curved edges while Rectangle does not
        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes([0.25, 0.25, 0.5, 0.5], projection=wcs, aspect="equal")
        ax.set_xlim(0, 10000)
        ax.set_ylim(-10000, 0)

        # Add a quadrangle patch (100 degrees by 20 degrees)
        q = Quadrangle(
            (255, -70) * u.deg,
            100 * u.deg,
            20 * u.deg,
            label="Quadrangle",
            edgecolor="blue",
            facecolor="none",
            transform=ax.get_transform("icrs"),
        )
        ax.add_patch(q)

        # Add a rectangle patch (100 degrees by 20 degrees)
        r = Rectangle(
            (255, -70),
            100,
            20,
            label="Rectangle",
            edgecolor="red",
            facecolor="none",
            linestyle="--",
            transform=ax.get_transform("icrs"),
        )
        ax.add_patch(r)

        ax.coords[0].set_ticklabel_visible(False)
        ax.coords[1].set_ticklabel_visible(False)

        return fig

    @figure_test
    def test_beam_shape_from_args(self, tmp_path):
        # Test for adding the beam shape with the beam parameters as arguments
        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_axes([0.2, 0.2, 0.6, 0.6], projection=wcs, aspect="equal")
        ax.set_xlim(-10, 10)
        ax.set_ylim(-10, 10)

        add_beam(
            ax,
            major=2 * u.arcmin,
            minor=1 * u.arcmin,
            angle=-30 * u.degree,
            corner="bottom right",
            frame=True,
            borderpad=0.0,
            pad=1.0,
            color="black",
        )

        return fig

    @figure_test
    def test_beam_shape_from_header(self, tmp_path):
        # Test for adding the beam shape with the beam parameters from a header
        hdr = self.msx_header
        hdr["BMAJ"] = (2 * u.arcmin).to(u.degree).value
        hdr["BMIN"] = (1 * u.arcmin).to(u.degree).value
        hdr["BPA"] = 30.0

        wcs = WCS(hdr)
        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_axes([0.2, 0.2, 0.6, 0.6], projection=wcs, aspect="equal")
        ax.set_xlim(-10, 10)
        ax.set_ylim(-10, 10)

        add_beam(ax, header=hdr)

        return fig

    @figure_test
    def test_scalebar(self, tmp_path):
        # Test for adding a scale bar
        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_axes([0.2, 0.2, 0.6, 0.6], projection=wcs, aspect="equal")
        ax.set_xlim(-10, 10)
        ax.set_ylim(-10, 10)

        add_scalebar(
            ax,
            2 * u.arcmin,
            label="2'",
            corner="top right",
            borderpad=1.0,
            label_top=True,
        )

        return fig

    @figure_test
    def test_elliptical_frame(self):
        # Regression test for a bug (astropy/astropy#6063) that caused labels to
        # be incorrectly simplified.

        wcs = WCS(self.msx_header)
        fig = plt.figure(figsize=(5, 3))
        fig.add_axes([0.2, 0.2, 0.6, 0.6], projection=wcs, frame_class=EllipticalFrame)
        return fig

    @figure_test
    def test_hms_labels(self):
        # This tests the apparance of the hms superscripts in tick labels
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes(
            [0.3, 0.2, 0.65, 0.6], projection=WCS(self.twoMASS_k_header), aspect="equal"
        )
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.5, 0.5)
        ax.coords[0].set_ticks(spacing=0.2 * 15 * u.arcsec)
        return fig

    @figure_test(style={"text.usetex": True})
    def test_latex_labels(self):
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_axes(
            [0.3, 0.2, 0.65, 0.6], projection=WCS(self.twoMASS_k_header), aspect="equal"
        )
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.5, 0.5)
        ax.coords[0].set_ticks(spacing=0.2 * 15 * u.arcsec)
        return fig

    @figure_test
    def test_tick_params(self):
        # This is a test to make sure that tick_params works correctly. We try
        # and test as much as possible with a single reference image.

        wcs = WCS()
        wcs.wcs.ctype = ["lon", "lat"]

        fig = plt.figure(figsize=(6, 6))

        # The first subplot tests:
        # - that plt.tick_params works
        # - that by default both axes are changed
        # - changing the tick direction and appearance, the label appearance and padding
        ax = fig.add_subplot(2, 2, 1, projection=wcs)
        plt.tick_params(
            direction="in",
            length=20,
            width=5,
            pad=6,
            labelsize=6,
            color="red",
            labelcolor="blue",
        )

        ax.coords[0].set_auto_axislabel(False)
        ax.coords[1].set_auto_axislabel(False)

        # The second subplot tests:
        # - that specifying grid parameters doesn't actually cause the grid to
        #   be shown (as expected)
        # - that axis= can be given integer coordinates or their string name
        # - that the tick positioning works (bottom/left/top/right)
        # Make sure that we can pass things that can index coords
        ax = fig.add_subplot(2, 2, 2, projection=wcs)
        plt.tick_params(
            axis=0,
            direction="in",
            length=20,
            width=5,
            pad=4,
            labelsize=6,
            color="red",
            labelcolor="blue",
            bottom=True,
            grid_color="purple",
        )
        plt.tick_params(
            axis="lat",
            direction="out",
            labelsize=8,
            color="blue",
            labelcolor="purple",
            left=True,
            right=True,
            grid_color="red",
        )

        ax.coords[0].set_auto_axislabel(False)
        ax.coords[1].set_auto_axislabel(False)

        # The third subplot tests:
        # - that ax.tick_params works
        # - that the grid has the correct settings once shown explicitly
        # - that we can use axis='x' and axis='y'
        ax = fig.add_subplot(2, 2, 3, projection=wcs)
        ax.tick_params(
            axis="x",
            direction="in",
            length=20,
            width=5,
            pad=20,
            labelsize=6,
            color="red",
            labelcolor="blue",
            bottom=True,
            grid_color="purple",
        )
        ax.tick_params(
            axis="y",
            direction="out",
            labelsize=8,
            color="blue",
            labelcolor="purple",
            left=True,
            right=True,
            grid_color="red",
        )
        plt.grid()

        ax.coords[0].set_auto_axislabel(False)
        ax.coords[1].set_auto_axislabel(False)

        # The final subplot tests:
        # - that we can use tick_params on a specific coordinate
        # - that the label positioning can be customized
        # - that the colors argument works
        # - that which='minor' works
        ax = fig.add_subplot(2, 2, 4, projection=wcs)
        ax.coords[0].tick_params(
            length=4,
            pad=2,
            colors="orange",
            labelbottom=True,
            labeltop=True,
            labelsize=10,
        )
        ax.coords[1].display_minor_ticks(True)
        ax.coords[1].tick_params(which="minor", length=6)

        ax.coords[0].set_auto_axislabel(False)
        ax.coords[1].set_auto_axislabel(False)

        return fig


@pytest.fixture
def wave_wcs_1d():
    wcs = WCS(naxis=1)
    wcs.wcs.ctype = ["WAVE"]
    wcs.wcs.cunit = ["m"]
    wcs.wcs.crpix = [1]
    wcs.wcs.cdelt = [5]
    wcs.wcs.crval = [45]
    wcs.wcs.set()
    return wcs


@figure_test
def test_1d_plot_1d_wcs(wave_wcs_1d):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=wave_wcs_1d)
    (lines,) = ax.plot([10, 12, 14, 12, 10])

    ax.set_xlabel("this is the x-axis")
    ax.set_ylabel("this is the y-axis")

    return fig


@figure_test
def test_1d_plot_1d_wcs_format_unit(wave_wcs_1d):
    """
    This test ensures that the format unit is updated and displayed for both
    the axis ticks and default axis labels.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=wave_wcs_1d)
    (lines,) = ax.plot([10, 12, 14, 12, 10])

    ax.coords[0].set_format_unit("nm")

    return fig


@pytest.fixture
def spatial_wcs_2d():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]
    wcs.wcs.crpix = [3.0] * 2
    wcs.wcs.cdelt = [15] * 2
    wcs.wcs.crval = [50.0] * 2
    wcs.wcs.set()
    return wcs


@figure_test
def test_1d_plot_2d_wcs_correlated(spatial_wcs_2d):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=spatial_wcs_2d, slices=("x", 0))
    (lines,) = ax.plot([10, 12, 14, 12, 10], "-o", color="orange")

    ax.coords["glon"].set_ticks(color="red")
    ax.coords["glon"].set_ticklabel(color="red")
    ax.coords["glon"].grid(color="red")

    ax.coords["glat"].set_ticks(color="blue")
    ax.coords["glat"].set_ticklabel(color="blue")
    ax.coords["glat"].grid(color="blue")

    return fig


@pytest.fixture
def spatial_wcs_2d_small_angle():
    """
    This WCS has an almost linear correlation between the pixel and world axes
    close to the reference pixel.
    """
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["HPLN-TAN", "HPLT-TAN"]
    wcs.wcs.crpix = [3.0] * 2
    wcs.wcs.cdelt = [10 / 3600, 5 / 3600]
    wcs.wcs.crval = [0] * 2
    wcs.wcs.set()
    return wcs


@pytest.mark.parametrize(
    "slices, bottom_axis",
    [
        # Remember SLLWCS takes slices in array order
        (np.s_[0, :], "custom:pos.helioprojective.lon"),
        (np.s_[:, 0], "custom:pos.helioprojective.lat"),
    ],
)
@figure_test
def test_1d_plot_1d_sliced_low_level_wcs(
    spatial_wcs_2d_small_angle, slices, bottom_axis
):
    """
    Test that a SLLWCS through a coupled 2D WCS plots as line OK.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=spatial_wcs_2d_small_angle[slices])
    (lines,) = ax.plot([10, 12, 14, 12, 10], "-o", color="orange")

    # Draw to trigger rendering the ticks.
    plt.draw()

    assert ax.coords[bottom_axis].ticks.get_visible_axes() == ["b"]

    return fig


@pytest.mark.parametrize(
    "slices, bottom_axis", [(("x", 0), "hpln"), ((0, "x"), "hplt")]
)
@figure_test
def test_1d_plot_put_varying_axis_on_bottom_lon(
    spatial_wcs_2d_small_angle, slices, bottom_axis
):
    """
    When we plot a 1D slice through spatial axes, we want to put the axis which
    actually changes on the bottom.

    For example an aligned wcs, pixel grid where you plot a lon slice through a
    lat axis, you would end up with no ticks on the bottom as the lon doesn't
    change, and a set of lat ticks on the top because it does but it's the
    correlated axis not the actual one you are plotting against.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=spatial_wcs_2d_small_angle, slices=slices)
    (lines,) = ax.plot([10, 12, 14, 12, 10], "-o", color="orange")

    # Draw to trigger rendering the ticks.
    plt.draw()

    assert ax.coords[bottom_axis].ticks.get_visible_axes() == ["b"]

    return fig


@figure_test
def test_allsky_labels_wrap():
    # Regression test for a bug that caused some tick labels to not be shown
    # when looking at all-sky maps in the case where coord_wrap < 360

    fig = plt.figure(figsize=(4, 4))

    icen = 0

    for ctype in [("GLON-CAR", "GLAT-CAR"), ("HGLN-CAR", "HGLT-CAR")]:
        for cen in [0, 90, 180, 270]:
            icen += 1

            wcs = WCS(naxis=2)
            wcs.wcs.ctype = ctype
            wcs.wcs.crval = cen, 0
            wcs.wcs.crpix = 360.5, 180.5
            wcs.wcs.cdelt = -0.5, 0.5

            ax = fig.add_subplot(8, 1, icen, projection=wcs)
            ax.set_xlim(-0.5, 719.5)
            ax.coords[0].set_ticks(spacing=50 * u.deg)
            ax.coords[0].set_ticks_position("b")
            ax.coords[0].set_auto_axislabel(False)
            ax.coords[1].set_auto_axislabel(False)
            ax.coords[1].set_ticklabel_visible(False)
            ax.coords[1].set_ticks_visible(False)

    fig.subplots_adjust(hspace=2, left=0.05, right=0.95, bottom=0.1, top=0.95)

    return fig


@figure_test
def test_tickable_gridlines():
    wcs = WCS(
        {
            "naxis": 2,
            "naxis1": 360,
            "naxis2": 180,
            "crpix1": 180.5,
            "crpix2": 90.5,
            "cdelt1": -1,
            "cdelt2": 1,
            "ctype1": "RA---CAR",
            "ctype2": "DEC--CAR",
        }
    )

    fig = Figure()
    ax = fig.add_subplot(projection=wcs)
    ax.set_xlim(-0.5, 360 - 0.5)
    ax.set_ylim(-0.5, 150 - 0.5)

    lon, lat = ax.coords
    lon.grid()
    lat.grid()

    overlay = ax.get_coords_overlay("galactic")
    overlay[0].set_ticks(spacing=30 * u.deg)
    overlay[1].set_ticks(spacing=30 * u.deg)

    # Test both single-character and multi-character names
    overlay[1].add_tickable_gridline("g", -30 * u.deg)
    overlay[0].add_tickable_gridline("const-glon", 30 * u.deg)

    overlay[0].grid(color="magenta")
    overlay[0].set_ticklabel_position("gt")
    overlay[0].set_ticklabel(color="magenta")
    overlay[0].set_axislabel("Galactic longitude", color="magenta")

    overlay[1].grid(color="blue")
    overlay[1].set_ticklabel_position(("const-glon", "r"))
    overlay[1].set_ticklabel(color="blue")
    overlay[1].set_axislabel("Galactic latitude", color="blue")

    return fig

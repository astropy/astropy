# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from matplotlib.figure import Figure

from astropy.tests.figures import figure_test
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization.wcsaxes.frame import BaseFrame
from astropy.wcs import WCS

from .test_images import BaseImageTests


class HexagonalFrame(BaseFrame):
    spine_names = "abcdef"

    def update_spines(self):
        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        ymid = 0.5 * (ymin + ymax)
        xmid1 = (xmin + xmax) / 4.0
        xmid2 = (xmin + xmax) * 3.0 / 4.0

        self["a"].data = np.array(([xmid1, ymin], [xmid2, ymin]))
        self["b"].data = np.array(([xmid2, ymin], [xmax, ymid]))
        self["c"].data = np.array(([xmax, ymid], [xmid2, ymax]))
        self["d"].data = np.array(([xmid2, ymax], [xmid1, ymax]))
        self["e"].data = np.array(([xmid1, ymax], [xmin, ymid]))
        self["f"].data = np.array(([xmin, ymid], [xmid1, ymin]))


class TestFrame(BaseImageTests):
    @figure_test(tolerance=0.5)
    def test_custom_frame(self):
        wcs = WCS(self.msx_header)

        fig = Figure(figsize=(4, 4))
        ax = WCSAxes(fig, [0.15, 0.15, 0.7, 0.7], wcs=wcs, frame_class=HexagonalFrame)
        fig.add_axes(ax)

        ax.coords.grid(color="white")

        im = ax.imshow(
            np.ones((149, 149)),
            vmin=0.0,
            vmax=2.0,
            origin="lower",
            cmap="gist_heat",
        )

        minpad = {}
        minpad["a"] = minpad["d"] = 1
        minpad["b"] = minpad["c"] = minpad["e"] = minpad["f"] = 2.75

        ax.coords["glon"].set_axislabel("Longitude", minpad=minpad)
        ax.coords["glon"].set_axislabel_position("ad")

        ax.coords["glat"].set_axislabel("Latitude", minpad=minpad)
        ax.coords["glat"].set_axislabel_position("bcef")

        ax.coords["glon"].set_ticklabel_position("ad")
        ax.coords["glat"].set_ticklabel_position("bcef")

        # Set limits so that no labels overlap
        ax.set_xlim(5.5, 100.5)
        ax.set_ylim(5.5, 110.5)

        # Clip the image to the frame
        im.set_clip_path(ax.coords.frame.patch)

        return fig

    @figure_test
    def test_update_clip_path_rectangular(self, tmp_path):
        fig = Figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], aspect="equal")

        fig.add_axes(ax)

        ax.set_xlim(0.0, 2.0)
        ax.set_ylim(0.0, 2.0)

        # Force drawing, which freezes the clip path returned by WCSAxes
        fig.savefig(tmp_path / "nothing")

        ax.imshow(np.zeros((12, 4)))

        ax.set_xlim(-0.5, 3.5)
        ax.set_ylim(-0.5, 11.5)

        ax.coords[0].set_auto_axislabel(False)
        ax.coords[1].set_auto_axislabel(False)

        return fig

    @figure_test
    def test_update_clip_path_nonrectangular(self, tmp_path):
        fig = Figure()
        ax = WCSAxes(
            fig, [0.1, 0.1, 0.8, 0.8], aspect="equal", frame_class=HexagonalFrame
        )

        fig.add_axes(ax)

        ax.set_xlim(0.0, 2.0)
        ax.set_ylim(0.0, 2.0)

        # Force drawing, which freezes the clip path returned by WCSAxes
        fig.savefig(tmp_path / "nothing")

        ax.imshow(np.zeros((12, 4)))

        ax.set_xlim(-0.5, 3.5)
        ax.set_ylim(-0.5, 11.5)

        return fig

    @figure_test
    def test_update_clip_path_change_wcs(self, tmp_path):
        # When WCS is changed, a new frame is created, so we need to make sure
        # that the path is carried over to the new frame.

        fig = Figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], aspect="equal")

        fig.add_axes(ax)

        ax.set_xlim(0.0, 2.0)
        ax.set_ylim(0.0, 2.0)

        # Force drawing, which freezes the clip path returned by WCSAxes
        fig.savefig(tmp_path / "nothing")

        ax.reset_wcs()

        ax.imshow(np.zeros((12, 4)))

        ax.set_xlim(-0.5, 3.5)
        ax.set_ylim(-0.5, 11.5)

        ax.coords[0].set_auto_axislabel(False)
        ax.coords[1].set_auto_axislabel(False)

        return fig

    def test_copy_frame_properties_change_wcs(self):
        # When WCS is changed, a new frame is created, so we need to make sure
        # that the color and linewidth are transferred over

        fig = Figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
        fig.add_axes(ax)
        ax.coords.frame.set_linewidth(5)
        ax.coords.frame.set_color("purple")
        ax.reset_wcs()
        assert ax.coords.frame.get_linewidth() == 5
        assert ax.coords.frame.get_color() == "purple"

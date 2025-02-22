# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np
from matplotlib.figure import Figure

from astropy import units as u
from astropy.tests.figures import figure_test
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization.wcsaxes.transforms import CurvedTransform
from astropy.wcs import WCS

from .test_images import BaseImageTests

# Create fake transforms that roughly mimic a polar projection


class DistanceToLonLat(CurvedTransform):
    has_inverse = True

    def __init__(self, R=6e3):
        super().__init__()
        self.R = R

    def transform(self, xy):
        x, y = xy[:, 0], xy[:, 1]
        lam = np.degrees(np.arctan2(y, x))
        phi = 90.0 - np.degrees(np.hypot(x, y) / self.R)
        return np.array((lam, phi)).transpose()

    transform_non_affine = transform

    def inverted(self):
        return LonLatToDistance(R=self.R)


class LonLatToDistance(CurvedTransform):
    def __init__(self, R=6e3):
        super().__init__()
        self.R = R

    def transform(self, lamphi):
        lam, phi = lamphi[:, 0], lamphi[:, 1]
        r = np.radians(90 - phi) * self.R
        x = r * np.cos(np.radians(lam))
        y = r * np.sin(np.radians(lam))
        return np.array((x, y)).transpose()

    transform_non_affine = transform

    def inverted(self):
        return DistanceToLonLat(R=self.R)


class TestTransformCoordMeta(BaseImageTests):
    @figure_test
    def test_coords_overlay(self):
        # Set up a simple WCS that maps pixels to non-projected distances
        wcs = WCS(naxis=2)
        wcs.wcs.ctype = ["x", "y"]
        wcs.wcs.cunit = ["km", "km"]
        wcs.wcs.crpix = [614.5, 856.5]
        wcs.wcs.cdelt = [6.25, 6.25]
        wcs.wcs.crval = [0.0, 0.0]

        fig = Figure(figsize=(4, 4))

        ax = WCSAxes(fig, [0.15, 0.15, 0.7, 0.7], wcs=wcs)
        fig.add_axes(ax)

        s = DistanceToLonLat(R=6378.273)

        ax.coords["x"].set_ticklabel_position("")
        ax.coords["y"].set_ticklabel_position("")

        coord_meta = {}
        coord_meta["type"] = ("longitude", "latitude")
        coord_meta["wrap"] = (360.0 * u.deg, None)
        coord_meta["unit"] = (u.deg, u.deg)
        coord_meta["name"] = "lon", "lat"

        overlay = ax.get_coords_overlay(s, coord_meta=coord_meta)

        overlay.grid(color="red")
        overlay["lon"].grid(color="red", linestyle="solid", alpha=0.3)
        overlay["lat"].grid(color="blue", linestyle="solid", alpha=0.3)

        overlay["lon"].set_ticklabel(size=7, exclude_overlapping=True)
        overlay["lat"].set_ticklabel(size=7, exclude_overlapping=True)

        overlay["lon"].set_ticklabel_position("brtl")
        overlay["lat"].set_ticklabel_position("brtl")

        overlay["lon"].set_ticks(spacing=10.0 * u.deg)
        overlay["lat"].set_ticks(spacing=10.0 * u.deg)

        ax.set_xlim(-0.5, 1215.5)
        ax.set_ylim(-0.5, 1791.5)

        return fig

    @figure_test
    def test_coords_overlay_auto_coord_meta(self):
        fig = Figure(figsize=(4, 4))

        ax = WCSAxes(fig, [0.15, 0.15, 0.7, 0.7], wcs=WCS(self.msx_header))
        fig.add_axes(ax)

        ax.grid(color="red", alpha=0.5, linestyle="solid")

        overlay = ax.get_coords_overlay("fk5")  # automatically sets coord_meta

        overlay.grid(color="black", alpha=0.5, linestyle="solid")

        overlay["ra"].set_ticks(color="black")
        overlay["dec"].set_ticks(color="black")

        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)

        return fig

    @figure_test
    def test_direct_init(self):
        s = DistanceToLonLat(R=6378.273)

        coord_meta = {}
        coord_meta["type"] = ("longitude", "latitude")
        coord_meta["wrap"] = (360.0 * u.deg, None)
        coord_meta["unit"] = (u.deg, u.deg)
        coord_meta["name"] = "lon", "lat"
        fig = Figure(figsize=(4, 4))

        ax = WCSAxes(fig, [0.15, 0.15, 0.7, 0.7], transform=s, coord_meta=coord_meta)
        fig.add_axes(ax)

        ax.coords["lon"].grid(color="red", linestyle="solid", alpha=0.3)
        ax.coords["lat"].grid(color="blue", linestyle="solid", alpha=0.3)

        ax.coords["lon"].set_auto_axislabel(False)
        ax.coords["lat"].set_auto_axislabel(False)

        ax.coords["lon"].set_ticklabel(size=7, exclude_overlapping=True)
        ax.coords["lat"].set_ticklabel(size=7, exclude_overlapping=True)

        ax.coords["lon"].set_ticklabel_position("brtl")
        ax.coords["lat"].set_ticklabel_position("brtl")

        ax.coords["lon"].set_ticks(spacing=10.0 * u.deg)
        ax.coords["lat"].set_ticks(spacing=10.0 * u.deg)

        ax.set_xlim(-400.0, 500.0)
        ax.set_ylim(-300.0, 400.0)

        return fig

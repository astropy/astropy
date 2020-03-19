# Licensed under a 3-clause BSD style license - see LICENSE.rst
import matplotlib.pyplot as plt
import pytest
from matplotlib.backend_bases import KeyEvent

import numpy as np

import astropy.units as u
from astropy.coordinates import FK5, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization.wcsaxes.core import WCSAxes
from astropy.wcs import WCS
from astropy.coordinates import galactocentric_frame_defaults

from .test_images import BaseImageTests


class TestDisplayWorldCoordinate(BaseImageTests):

    def teardown_method(self, method):
        plt.close('all')

    def test_overlay_coords(self, ignore_matplotlibrc, tmpdir):
        wcs = WCS(self.msx_header)

        fig = plt.figure(figsize=(4, 4))
        canvas = fig.canvas

        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
        fig.add_axes(ax)

        # On some systems, fig.canvas.draw is not enough to force a draw, so we
        # save to a temporary file.
        fig.savefig(tmpdir.join('test1.png').strpath)

        # Testing default displayed world coordinates
        string_world = ax._display_world_coords(0.523412, 0.518311)
        assert string_world == '0\xb029\'45" -0\xb029\'20" (world)'

        # Test pixel coordinates
        event1 = KeyEvent('test_pixel_coords', canvas, 'w')
        fig.canvas.key_press_event(event1.key, guiEvent=event1)
        string_pixel = ax._display_world_coords(0.523412, 0.523412)
        assert string_pixel == "0.523412 0.523412 (pixel)"

        event3 = KeyEvent('test_pixel_coords', canvas, 'w')
        fig.canvas.key_press_event(event3.key, guiEvent=event3)
        # Test that it still displays world coords when there are no overlay coords
        string_world2 = ax._display_world_coords(0.523412, 0.518311)
        assert string_world2 == '0\xb029\'45" -0\xb029\'20" (world)'

        overlay = ax.get_coords_overlay('fk5')

        # Regression test for bug that caused format to always be taken from
        # main world coordinates.
        overlay[0].set_major_formatter('d.ddd')

        # On some systems, fig.canvas.draw is not enough to force a draw, so we
        # save to a temporary file.
        fig.savefig(tmpdir.join('test2.png').strpath)

        event4 = KeyEvent('test_pixel_coords', canvas, 'w')
        fig.canvas.key_press_event(event4.key, guiEvent=event4)
        # Test that it displays the overlay world coordinates
        string_world3 = ax._display_world_coords(0.523412, 0.518311)

        assert string_world3 == '267.176\xb0 -28\xb045\'56" (world, overlay 1)'

        overlay = ax.get_coords_overlay(FK5())

        # Regression test for bug that caused format to always be taken from
        # main world coordinates.
        overlay[0].set_major_formatter('d.ddd')

        # On some systems, fig.canvas.draw is not enough to force a draw, so we
        # save to a temporary file.
        fig.savefig(tmpdir.join('test3.png').strpath)

        event5 = KeyEvent('test_pixel_coords', canvas, 'w')
        fig.canvas.key_press_event(event4.key, guiEvent=event4)
        # Test that it displays the overlay world coordinates
        string_world4 = ax._display_world_coords(0.523412, 0.518311)

        assert string_world4 == '267.176\xb0 -28\xb045\'56" (world, overlay 2)'

        overlay = ax.get_coords_overlay(FK5(equinox=Time("J2030")))

        # Regression test for bug that caused format to always be taken from
        # main world coordinates.
        overlay[0].set_major_formatter('d.ddd')

        # On some systems, fig.canvas.draw is not enough to force a draw, so we
        # save to a temporary file.
        fig.savefig(tmpdir.join('test4.png').strpath)

        event6 = KeyEvent('test_pixel_coords', canvas, 'w')
        fig.canvas.key_press_event(event5.key, guiEvent=event6)
        # Test that it displays the overlay world coordinates
        string_world5 = ax._display_world_coords(0.523412, 0.518311)

        assert string_world5 == '267.652\xb0 -28\xb046\'23" (world, overlay 3)'

    def test_cube_coords(self, ignore_matplotlibrc, tmpdir):
        wcs = WCS(self.cube_header)

        fig = plt.figure(figsize=(4, 4))
        canvas = fig.canvas

        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs, slices=('y', 50, 'x'))
        fig.add_axes(ax)

        # On some systems, fig.canvas.draw is not enough to force a draw, so we
        # save to a temporary file.
        fig.savefig(tmpdir.join('test.png').strpath)

        # Testing default displayed world coordinates
        string_world = ax._display_world_coords(0.523412, 0.518311)
        assert string_world == '3h26m52.0s 30\xb037\'17\" 2563 (world)'

        # Test pixel coordinates
        event1 = KeyEvent('test_pixel_coords', canvas, 'w')
        fig.canvas.key_press_event(event1.key, guiEvent=event1)
        string_pixel = ax._display_world_coords(0.523412, 0.523412)
        assert string_pixel == "0.523412 0.523412 (pixel)"

    def test_cube_coords_uncorr_slicing(self, ignore_matplotlibrc, tmpdir):

        # Regression test for a bug that occurred with coordinate formatting if
        # some dimensions were uncorrelated and sliced out.

        wcs = WCS(self.cube_header)

        fig = plt.figure(figsize=(4, 4))
        canvas = fig.canvas

        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs, slices=('x', 'y', 2))
        fig.add_axes(ax)

        # On some systems, fig.canvas.draw is not enough to force a draw, so we
        # save to a temporary file.
        fig.savefig(tmpdir.join('test.png').strpath)

        # Testing default displayed world coordinates
        string_world = ax._display_world_coords(0.523412, 0.518311)
        assert string_world == '3h26m56.6s 30\xb018\'19\" (world)'

        # Test pixel coordinates
        event1 = KeyEvent('test_pixel_coords', canvas, 'w')
        fig.canvas.key_press_event(event1.key, guiEvent=event1)
        string_pixel = ax._display_world_coords(0.523412, 0.523412)
        assert string_pixel == "0.523412 0.523412 (pixel)"

    def test_plot_coord_3d_transform(self):
        wcs = WCS(self.msx_header)

        with galactocentric_frame_defaults.set('latest'):
            coord = SkyCoord(0 * u.kpc, 0 * u.kpc, 0 * u.kpc, frame='galactocentric')

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        point, = ax.plot_coord(coord, 'ro')

        np.testing.assert_allclose(point.get_xydata()[0], [0, 0], atol=1e-4)

# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ..core import WCSAxes
import matplotlib.pyplot as plt
from matplotlib.backend_bases import KeyEvent

from astropy.wcs import WCS
from astropy.extern import six
from astropy.tests.helper import pytest

from .test_images import BaseImageTests


class TestDisplayWorldCoordinate(BaseImageTests):

    def test_overlay_coords(self):
        wcs = WCS(self.msx_header)

        fig = plt.figure(figsize=(4,4))
        canvas = fig.canvas

        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
        fig.add_axes(ax)
        fig.canvas.draw()

        # Testing default displayed world coordinates
        string_world = ax._display_world_coords(0.523412, 0.518311)
        assert string_world == six.u('0\xb029\'45" -0\xb029\'20" (world)')

        # Test pixel coordinates
        event1 = KeyEvent('test_pixel_coords', canvas, 'p')
        fig.canvas.key_press_event(event1.key, guiEvent=event1)
        string_pixel = ax._display_world_coords(0.523412, 0.523412)
        assert string_pixel == "x=0.523412 y=0.523412"

        event2 = KeyEvent('test_pixel_coords', canvas, 'p')
        fig.canvas.key_press_event(event2.key, guiEvent=event2)

        event3 = KeyEvent('test_pixel_coords', canvas, 'o')
        fig.canvas.key_press_event(event3.key, guiEvent=event3)
        # Test that it still displays world coords when there are no overlay coords
        string_world2 = ax._display_world_coords(0.523412, 0.518311)
        assert string_world2 == six.u('0\xb029\'45" -0\xb029\'20" (world)')

        overlay = ax.get_coords_overlay('fk5')
        fig.canvas.draw()
        event4 = KeyEvent('test_pixel_coords', canvas, 'o')
        fig.canvas.key_press_event(event4.key, guiEvent=event4)
        # Test that it displays the overlay world coordinates
        string_world3 = ax._display_world_coords(0.523412, 0.518311)
        assert string_world3 == six.u('0\xb029\'45" -0\xb029\'20" (world), 267\xb010\'32" 331\xb014\'04" (overlay coords)')

    def test_cube_coords(self):
        wcs = WCS(self.cube_header)

        fig = plt.figure(figsize=(4,4))
        canvas = fig.canvas

        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs, slices=('y', 50, 'x'))
        fig.add_axes(ax)
        fig.canvas.draw()

        # Testing default displayed world coordinates
        string_world = ax._display_world_coords(0.523412, 0.518311)
        assert string_world == six.u('2563 51\xb043\'01" (world)')

        # Test pixel coordinates
        event1 = KeyEvent('test_pixel_coords', canvas, 'p')
        fig.canvas.key_press_event(event1.key, guiEvent=event1)
        string_pixel = ax._display_world_coords(0.523412, 0.523412)
        assert string_pixel == "x=0.523412 y=0.523412"

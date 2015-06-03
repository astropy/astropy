# Licensed under a 3-clause BSD style license - see LICENSE.rst

from matplotlib.axes import Axes, subplot_class_factory
from matplotlib.transforms import Affine2D, Bbox, Transform

from astropy.wcs import WCS

from .transforms import (WCSPixel2WorldTransform, WCSWorld2PixelTransform,
                         CoordinateTransform)
from .coordinates_map import CoordinatesMap
from .utils import get_coord_meta
from .wcs_utils import wcs_to_celestial_frame
from .frame import RectangularFrame
import numpy as np

__all__ = ['WCSAxes', 'WCSAxesSubplot']

VISUAL_PROPERTIES = ['facecolor', 'edgecolor', 'linewidth', 'alpha', 'linestyle']

IDENTITY = WCS(naxis=2)
IDENTITY.wcs.ctype = ["X", "Y"]
IDENTITY.wcs.crval = [1., 1.]
IDENTITY.wcs.crpix = [1., 1.]
IDENTITY.wcs.cdelt = [1., 1.]


class WCSAxes(Axes):

    def __init__(self, fig, rect, wcs=None, transform=None, coord_meta=None,
                 transData=None, slices=None, frame_class=RectangularFrame,
                 **kwargs):

        super(WCSAxes, self).__init__(fig, rect, **kwargs)
        self._bboxes = []

        self.frame_class = frame_class

        if not (transData is None):
            # User wants to override the transform for the final
            # data->pixel mapping
            self.transData = transData

        self.reset_wcs(wcs=wcs, slices=slices, transform=transform, coord_meta=coord_meta)
        self._hide_parent_artists()
        self.format_coord = self._display_world_coords
        self._display_coords_index = 0
        fig.canvas.mpl_connect('key_press_event', self._set_cursor_prefs)
        self.patch = self.coords.frame.patch

    def _display_world_coords(self, x, y):

        if self._display_coords_index == -1:
            return "%s %s (pixel)" % (x, y)

        pixel = np.array([x, y])

        coords = self._all_coords[self._display_coords_index]

        world = coords._transform.transform(np.array([pixel]))[0]

        xw = coords[self._x_index].format_coord(world[self._x_index])
        yw = coords[self._y_index].format_coord(world[self._y_index])

        if self._display_coords_index == 0:
            system = "world"
        else:
            system = "world, overlay {0}".format(self._display_coords_index)

        coord_string = "%s %s (%s)" % (xw, yw, system)

        return coord_string

    def _set_cursor_prefs(self, event, **kwargs):
        if event.key == 'w':
            self._display_coords_index += 1
            if self._display_coords_index + 1 > len(self._all_coords):
                self._display_coords_index = -1

    def _hide_parent_artists(self):
        # Turn off spines and current axes
        for s in self.spines.values():
            s.set_visible(False)

        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)

    def reset_wcs(self, wcs=None, slices=None, transform=None, coord_meta=None):
        """
        Reset the current Axes, to use a new WCS object.
        """

        # Here determine all the coordinate axes that should be shown.
        if wcs is None and transform is None:
            self.wcs = IDENTITY
        else:
            self.wcs = wcs

        # If we are making a new WCS, we need to preserve the path object since
        # it may already be used by objects that have been plotted, and we need
        # to continue updating it. CoordinatesMap will create a new frame
        # instance, but we can tell that instance to keep using the old path.
        if hasattr(self, 'coords'):
            previous_frame_path = self.coords.frame._path
        else:
            previous_frame_path = None

        self.coords = CoordinatesMap(self, wcs=self.wcs, slice=slices,
                                     transform=transform, coord_meta=coord_meta,
                                     frame_class=self.frame_class,
                                     previous_frame_path=previous_frame_path)

        self._all_coords = [self.coords]

        if slices is None:
            self.slices = ('x', 'y')
            self._x_index = 0
            self._y_index = 1
        else:
            self.slices = slices
            self._x_index = self.slices.index('x')
            self._y_index = self.slices.index('y')

        # Common default settings for Rectangular Frame
        if self.frame_class is RectangularFrame:
            for coord_index in range(len(self.slices)):
                if self.slices[coord_index] == 'x':
                    self.coords[coord_index].set_axislabel_position('b')
                    self.coords[coord_index].set_ticklabel_position('b')
                elif self.slices[coord_index] == 'y':
                    self.coords[coord_index].set_axislabel_position('l')
                    self.coords[coord_index].set_ticklabel_position('l')
                else:
                    self.coords[coord_index].set_axislabel_position('')
                    self.coords[coord_index].set_ticklabel_position('')
                    self.coords[coord_index].set_ticks_position('')

    def draw(self, renderer, inframe=False):

        # We need to make sure that that frame path is up to date
        self.coords.frame._update_patch_path()

        super(WCSAxes, self).draw(renderer, inframe)

        # Here need to find out range of all coordinates, and update range for
        # each coordinate axis. For now, just assume it covers the whole sky.

        self._bboxes = []
        self._ticklabels_bbox = []
        visible_ticks = []

        for coords in self._all_coords:

            coords.frame.update()
            for coord in coords:
                coord._draw(renderer, bboxes=self._bboxes,
                            ticklabels_bbox=self._ticklabels_bbox)
                visible_ticks.extend(coord.ticklabels.get_visible_axes())

        for coords in self._all_coords:

            for coord in coords:
                coord._draw_axislabels(renderer, bboxes=self._bboxes,
                                       ticklabels_bbox=self._ticklabels_bbox,
                                       visible_ticks=visible_ticks)

        self.coords.frame.draw(renderer)

    def set_xlabel(self, label):
        self.coords[self._x_index].set_axislabel(label)

    def set_ylabel(self, label):
        self.coords[self._y_index].set_axislabel(label)

    def get_xlabel(self):
        return self.coords[self._x_index].get_axislabel()

    def get_ylabel(self):
        return self.coords[self._y_index].get_axislabel()

    def get_coords_overlay(self, frame, coord_meta=None):

        # Here we can't use get_transform because that deals with
        # pixel-to-pixel transformations when passing a WCS object.
        if isinstance(frame, WCS):
            coords = CoordinatesMap(self, frame, frame_class=self.frame_class)
        else:
            if coord_meta is None:
                coord_meta = get_coord_meta(frame)
            transform = self._get_transform_no_transdata(frame)
            coords = CoordinatesMap(self, transform=transform, coord_meta=coord_meta, frame_class=self.frame_class)

        self._all_coords.append(coords)

        # Common settings for overlay
        coords[0].set_axislabel_position('t')
        coords[1].set_axislabel_position('r')
        coords[0].set_ticklabel_position('t')
        coords[1].set_ticklabel_position('r')

        self.overlay_coords = coords

        return coords

    def get_transform(self, frame):
        """
        Return a transform from the specified frame to display coordinates.

        This does not include the transData transformation

        Parameters
        ----------
        frame : :class:`~astropy.wcs.WCS` or :class:`~matplotlib.transforms.Transform` or str
            The ``frame`` parameter can have several possible types:
                * :class:`~astropy.wcs.WCS` instance: assumed to be a
                  transformation from pixel to world coordinates, where the
                  world coordinates are the same as those in the WCS
                  transformation used for this ``WCSAxes`` instance. This is
                  used for example to show contours, since this involves
                  plotting an array in pixel coordinates that are not the
                  final data coordinate and have to be transformed to the
                  common world coordinate system first.
                * :class:`~matplotlib.transforms.Transform` instance: it is
                  assumed to be a transform to the world coordinates that are
                  part of the WCS used to instantiate this ``WCSAxes``
                  instance.
                * ``'pixel'`` or ``'world'``: return a transformation that
                  allows users to plot in pixel/data coordinates (essentially
                  an identity transform) and ``world`` (the default
                  world-to-pixel transformation used to instantiate the
                  ``WCSAxes`` instance).
                * ``'fk5'`` or ``'galactic'``: return a transformation from
                  the specified frame to the pixel/data coordinates.
                * :class:`~astropy.coordinates.BaseCoordinateFrame` instance.
        """
        return self._get_transform_no_transdata(frame).inverted() + self.transData

    def _get_transform_no_transdata(self, frame):
        """
        Return a transform from data to the specified frame
        """

        if self.wcs is None and frame != 'pixel':
            raise ValueError('No WCS specified, so only pixel coordinates are available')

        if isinstance(frame, WCS):

            coord_in = wcs_to_celestial_frame(self.wcs)
            coord_out = wcs_to_celestial_frame(frame)

            if coord_in == coord_out:

                return (WCSPixel2WorldTransform(self.wcs, slice=self.slices)
                        + WCSWorld2PixelTransform(frame))

            else:

                return (WCSPixel2WorldTransform(self.wcs, slice=self.slices)
                        + CoordinateTransform(self.wcs, frame)
                        + WCSWorld2PixelTransform(frame))

        elif frame == 'pixel':

            return Affine2D()

        elif isinstance(frame, Transform):

            pixel2world = WCSPixel2WorldTransform(self.wcs, slice=self.slices)

            return pixel2world + frame

        else:

            pixel2world = WCSPixel2WorldTransform(self.wcs, slice=self.slices)

            if frame == 'world':

                return pixel2world

            else:
                coordinate_transform = CoordinateTransform(self.wcs, frame)

                if coordinate_transform.same_frames:
                    return pixel2world
                else:
                    return pixel2world + CoordinateTransform(self.wcs, frame)

    def get_tightbbox(self, renderer):

        if not self.get_visible():
            return

        bb = [b for b in self._bboxes if b and (b.width != 0 or b.height != 0)]

        if bb:
            _bbox = Bbox.union(bb)
            return _bbox
        else:
            return self.get_window_extent(renderer)

    def grid(self, draw_grid=True, **kwargs):
        """
        Plot gridlines for both coordinates.

        Standard matplotlib appearance options (color, alpha, etc.) can be
        passed as keyword arguments.

        Parameters
        ----------
        draw_grid : bool
            Whether to show the gridlines
        """
        if draw_grid and hasattr(self, 'coords'):
            self.coords.grid(draw_grid=draw_grid, **kwargs)

# In the following, we put the generated subplot class in a temporary class and
# we then inherit it - if we don't do this, the generated class appears to
# belong in matplotlib, not in WCSAxes, from the API's point of view.


class WCSAxesSubplot(subplot_class_factory(WCSAxes)):
    pass

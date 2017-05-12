# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function, division, absolute_import

import numpy as np

from matplotlib.axes import Axes, subplot_class_factory
from matplotlib.transforms import Affine2D, Bbox, Transform

from ...coordinates import SkyCoord, BaseCoordinateFrame
from ...wcs import WCS
from ...wcs.utils import wcs_to_celestial_frame
from ...extern import six

from .transforms import (WCSPixel2WorldTransform, WCSWorld2PixelTransform,
                         CoordinateTransform)
from .coordinates_map import CoordinatesMap
from .utils import get_coord_meta
from .frame import EllipticalFrame, RectangularFrame

__all__ = ['WCSAxes', 'WCSAxesSubplot']

VISUAL_PROPERTIES = ['facecolor', 'edgecolor', 'linewidth', 'alpha', 'linestyle']

IDENTITY = WCS(naxis=2)
IDENTITY.wcs.ctype = ["X", "Y"]
IDENTITY.wcs.crval = [1., 1.]
IDENTITY.wcs.crpix = [1., 1.]
IDENTITY.wcs.cdelt = [1., 1.]


class WCSAxes(Axes):
    """
    The main axes class that can be used to show world coordinates from a WCS.

    Parameters
    ----------
    fig : `~matplotlib.figure.Figure`
        The figure to add the axes to
    rect : list
        The position of the axes in the figure in relative units. Should be
        given as ``[left, bottom, width, height]``.
    wcs : :class:`~astropy.wcs.WCS`, optional
        The WCS for the data. If this is specified, ``transform`` cannot be
        specified.
    transform : `~matplotlib.transforms.Transform`, optional
        The transform for the data. If this is specified, ``wcs`` cannot be
        specified.
    coord_meta : dict, optional
        A dictionary providing additional metadata when ``transform`` is
        specified. This should include the keys ``type``, ``wrap``, and
        ``unit``. Each of these should be a list with as many items as the
        dimension of the WCS. The ``type`` entries should be one of
        ``longitude``, ``latitude``, or ``scalar``, the ``wrap`` entries should
        give, for the longitude, the angle at which the coordinate wraps (and
        `None` otherwise), and the ``unit`` should give the unit of the
        coordinates as :class:`~astropy.units.Unit` instances.
    transData : `~matplotlib.transforms.Transform`, optional
        Can be used to override the default data -> pixel mapping.
    slices : tuple, optional
        For WCS transformations with more than two dimensions, we need to
        choose which dimensions are being shown in the 2D image. The slice
        should contain one ``x`` entry, one ``y`` entry, and the rest of the
        values should be integers indicating the slice through the data. The
        order of the items in the slice should be the same as the order of the
        dimensions in the :class:`~astropy.wcs.WCS`, and the opposite of the
        order of the dimensions in Numpy. For example, ``(50, 'x', 'y')`` means
        that the first WCS dimension (last Numpy dimension) will be sliced at
        an index of 50, the second WCS and Numpy dimension will be shown on the
        x axis, and the final WCS dimension (first Numpy dimension) will be
        shown on the y-axis (and therefore the data will be plotted using
        ``data[:, :, 50].transpose()``)
    frame_class : type, optional
        The class for the frame, which should be a subclass of
        :class:`~astropy.visualization.wcsaxes.frame.BaseFrame`. The default is to use a
        :class:`~astropy.visualization.wcsaxes.frame.RectangularFrame`
    """

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
        self._drawn = False

    def _display_world_coords(self, x, y):

        if not self._drawn:
            return ""

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

    # We now overload ``imshow`` because we need to make sure that origin is
    # set to ``lower`` for all images, which means that we need to flip RGB
    # images.
    def imshow(self, X, *args, **kwargs):
        """
        Wrapper to Matplotlib's :meth:`~matplotlib.axes.Axes.imshow`.

        If an RGB image is passed as a PIL object, it will be flipped
        vertically and ``origin`` will be set to ``lower``, since WCS
        transformations - like FITS files - assume that the origin is the lower
        left pixel of the image (whereas RGB images have the origin in the top
        left).

        All arguments are passed to :meth:`~matplotlib.axes.Axes.imshow`.
        """

        origin = kwargs.get('origin', None)

        if origin == 'upper':
            raise ValueError("Cannot use images with origin='upper' in WCSAxes.")

        # To check whether the image is a PIL image we can check if the data
        # has a 'getpixel' attribute - this is what Matplotlib's AxesImage does

        try:
            from PIL.Image import Image, FLIP_TOP_BOTTOM
        except ImportError:
            # We don't need to worry since PIL is not installed, so user cannot
            # have passed RGB image.
            pass
        else:
            if isinstance(X, Image) or hasattr(X, 'getpixel'):
                X = X.transpose(FLIP_TOP_BOTTOM)
                kwargs['origin'] = 'lower'

        return super(WCSAxes, self).imshow(X, *args, **kwargs)

    def plot_coord(self, *args, **kwargs):
        """
        Plot `~astropy.coordinates.SkyCoord` or
        `~astropy.coordinates.BaseCoordinateFrame` objects onto the axes.

        The first argument to
        :meth:`~astropy.visualization.wcsaxes.WCSAxes.plot_coord` should be a
        coordinate, which will then be converted to the first two parameters to
        `matplotlib.axes.Axes.plot`. All other arguments are the same as
        `matplotlib.axes.Axes.plot`. If not specified a ``transform`` keyword
        argument will be created based on the coordinate.

        Parameters
        ----------
        coordinate : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate object to plot on the axes. This is converted to the
            first two arguments to `matplotlib.axes.Axes.plot`.

        See Also
        --------

        matplotlib.axes.Axes.plot : This method is called from this function with all arguments passed to it.

        """

        if isinstance(args[0], (SkyCoord, BaseCoordinateFrame)):

            # Extract the frame from the first argument.
            frame0 = args[0]
            if isinstance(frame0, SkyCoord):
                frame0 = frame0.frame

            plot_data = []
            for coord in self.coords:
                if coord.coord_type == 'longitude':
                    plot_data.append(frame0.data.lon.to(coord.coord_unit).value)
                elif coord.coord_type == 'latitude':
                    plot_data.append(frame0.data.lat.to(coord.coord_unit).value)
                else:
                    raise NotImplementedError("Coordinates cannot be plotted with this "
                                              "method because the WCS does not represent longitude/latitude.")

            if 'transform' in kwargs.keys():
                raise TypeError("The 'transform' keyword argument is not allowed,"
                                " as it is automatically determined by the input coordinate frame.")

            transform = self.get_transform(frame0)
            kwargs.update({'transform': transform})

            args = tuple(plot_data) + args[1:]

        super(WCSAxes, self).plot(*args, **kwargs)

    def reset_wcs(self, wcs=None, slices=None, transform=None, coord_meta=None):
        """
        Reset the current Axes, to use a new WCS object.
        """

        # Here determine all the coordinate axes that should be shown.
        if wcs is None and transform is None:

            self.wcs = IDENTITY

        else:

            # We now force call 'set', which ensures the WCS object is
            # consistent, which will only be important if the WCS has been set
            # by hand. For example if the user sets a celestial WCS by hand and
            # forgets to set the units, WCS.wcs.set() will do this.
            if wcs is not None:
                wcs.wcs.set()

            self.wcs = wcs

        # If we are making a new WCS, we need to preserve the path object since
        # it may already be used by objects that have been plotted, and we need
        # to continue updating it. CoordinatesMap will create a new frame
        # instance, but we can tell that instance to keep using the old path.
        if hasattr(self, 'coords'):
            previous_frame = {'path': self.coords.frame._path,
                              'color': self.coords.frame.get_color(),
                              'linewidth': self.coords.frame.get_linewidth()}
        else:
            previous_frame = {'path': None}

        self.coords = CoordinatesMap(self, wcs=self.wcs, slice=slices,
                                     transform=transform, coord_meta=coord_meta,
                                     frame_class=self.frame_class,
                                     previous_frame_path=previous_frame['path'])

        if previous_frame['path'] is not None:
            self.coords.frame.set_color(previous_frame['color'])
            self.coords.frame.set_linewidth(previous_frame['linewidth'])

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
        # Common default settings for Elliptical Frame
        elif self.frame_class is EllipticalFrame:
            for coord_index in range(len(self.slices)):
                if self.slices[coord_index] == 'x':
                    self.coords[coord_index].set_axislabel_position('h')
                    self.coords[coord_index].set_ticklabel_position('h')
                    self.coords[coord_index].set_ticks_position('h')
                elif self.slices[coord_index] == 'y':
                    self.coords[coord_index].set_ticks_position('c')
                    self.coords[coord_index].set_axislabel_position('c')
                    self.coords[coord_index].set_ticklabel_position('c')
                else:
                    self.coords[coord_index].set_axislabel_position('')
                    self.coords[coord_index].set_ticklabel_position('')
                    self.coords[coord_index].set_ticks_position('')

    def draw(self, renderer, inframe=False):

        # In Axes.draw, the following code can result in the xlim and ylim
        # values changing, so we need to force call this here to make sure that
        # the limits are correct before we update the patch.
        locator = self.get_axes_locator()
        if locator:
            pos = locator(self, renderer)
            self.apply_aspect(pos)
        else:
            self.apply_aspect()

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

        self._drawn = True

    def set_xlabel(self, label, labelpad=1, **kwargs):
        self.coords[self._x_index].set_axislabel(label, minpad=labelpad, **kwargs)

    def set_ylabel(self, label, labelpad=1, **kwargs):
        self.coords[self._y_index].set_axislabel(label, minpad=labelpad, **kwargs)

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
            coords = CoordinatesMap(self, transform=transform,
                                    coord_meta=coord_meta,
                                    frame_class=self.frame_class)

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

                return (WCSPixel2WorldTransform(self.wcs, slice=self.slices) +
                        WCSWorld2PixelTransform(frame))

            else:

                return (WCSPixel2WorldTransform(self.wcs, slice=self.slices) +
                        CoordinateTransform(self.wcs, frame) +
                        WCSWorld2PixelTransform(frame))

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

    def grid(self, b=None, axis='both', **kwargs):
        """
        Plot gridlines for both coordinates.

        Standard matplotlib appearance options (color, alpha, etc.) can be
        passed as keyword arguments. This behaves like `matplotlib.axes.Axes`
        except that if no arguments are specified, the grid is shown rather
        than toggled.

        Parameters
        ----------
        b : bool
            Whether to show the gridlines.
        """

        if not hasattr(self, 'coords'):
            return

        which = kwargs.pop('which', 'major')
        if which != 'major':
            raise NotImplementedError('Plotting the grid for the minor ticks is '
                                      'not supported.')

        if axis == 'both':
            self.coords.grid(draw_grid=b, **kwargs)
        elif axis == 'x':
            self.coords[0].grid(draw_grid=b, **kwargs)
        elif axis == 'y':
            self.coords[1].grid(draw_grid=b, **kwargs)
        else:
            raise ValueError('axis should be one of x/y/both')

# In the following, we put the generated subplot class in a temporary class and
# we then inherit it - if we don't do this, the generated class appears to
# belong in matplotlib, not in WCSAxes, from the API's point of view.


class WCSAxesSubplot(subplot_class_factory(WCSAxes)):
    """
    A subclass class for WCSAxes
    """
    pass

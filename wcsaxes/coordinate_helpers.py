"""
This file defines the classes used to represent a 'coordinate', which includes
axes, ticks, tick labels, and grid lines.
"""

import abc

import numpy as np

from matplotlib.ticker import Formatter
from matplotlib.transforms import Affine2D, ScaledTranslation
from matplotlib.collections import PathCollection
from matplotlib.text import Text

from .formatter_locator import AngleFormatterLocator, ScalarFormatterLocator
from .ticks import Ticks
from .ticklabels import TickLabels
from .grid_paths import get_lon_lat_path, get_gridline_path
from .utils import get_pixels_to_data_scales

from . import six


# TODO: very little here is sky-specific. Only the wrapping. Maybe we can avoid the use of an empty base class.


class BaseCoordinateHelper(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self, parent_axes=None, transform=None, coord_index=None):

        super(BaseCoordinateHelper, self).__init__()

        # Keep a reference to the parent axes and the transform
        self.parent_axes = parent_axes
        self.transform = transform
        self.coord_index = coord_index

        # Initialize ticks
        self.dpi_transform = Affine2D()
        self.offset_transform = ScaledTranslation(0, 0, self.dpi_transform)
        self.ticks = Ticks(transform=self.parent_axes.transData +
                           self.offset_transform)

        self.ticklabels = TickLabels(transform=self.parent_axes.transData,
                                     figure=self.parent_axes.get_figure())

        # Initialize formatter/locator
        self._formatter_locator = self._formatter_locator_class()

        # Initialize container for the grid lines
        self.grid_lines = PathCollection([],
                                         transform=self.parent_axes.transData,
                                         facecolor='none', visible=False,
                                         figure=self.parent_axes.get_figure())

    def grid(self, draw_grid=True, **kwargs):
        """
        Plot gridlines for this coordinate.

        Standard matplotlib appearance options (color, alpha, etc.) can be
        passed as keyword arguments.

        Parameters
        ----------
        draw_grid : bool
            Whether to show the gridlines
        """

        if 'color' in kwargs:
            kwargs['edgecolor'] = kwargs.pop('color')

        self.grid_lines.set(**kwargs)

        if self.grid_lines.get_visible():
            if not draw_grid:
                self.grid_lines.set_visible(False)
        else:
            self.grid_lines.set_visible(True)

    def set_major_formatter(self, formatter):
        """
        Set the formatter to use for the major tick labels.

        Parameters
        ----------
        formatter : str or Formatter
            The format or formatter to use.
        """
        if isinstance(formatter, Formatter):
            raise NotImplementedError()  # figure out how to swap out formatter
        elif isinstance(formatter, six.string_types):
            self._formatter_locator.format = formatter
        else:
            raise TypeError("formatter should be a string or a Formatter "
                            "instance")

    def set_ticks(self, values=None, spacing=None, number=None):
        """
        Set the location of the ticks.

        Only one of the options from ``values``, ``spacing``, or ``number``
        should be specified.

        Parameters
        ----------
        values : iterable, optional
            The coordinate values at which to show the ticks.
        spacing : float, optional
            The spacing between ticks.
        number : float, optional
            The approximate number of ticks shown.
        """
        if values is not None:
            self._formatter_locator.values = values
        elif spacing is not None:
            self._formatter_locator.spacing = spacing
        elif number is not None:
            self._formatter_locator.number = number
        else:
            raise ValueError("one of values, spacing, or number should be "
                             "specified")

    @property
    def locator(self):
        return _formatter_locator.locator

    @property
    def formatter(self):
        return _formatter_locator.formatter

    def draw(self, renderer):

        renderer.open_group('coordinate_axis')

        self._update_ticks(renderer)
        self.ticks.draw(renderer)
        self.ticklabels.draw(renderer)
        self.grid_lines.draw(renderer)

        renderer.close_group('coordinate_axis')

    def set_ticks_position(self):
        """
        Set the axes on which the ticks for this coordinate should
        appear. Should be a string containing zero or more of ``b``, ``t``,
        ``l``, ``r``.
        """
        raise NotImplementedError()

    def set_ticklabel_position(self, position):
        """
        Set the axes on which the ticklabels for this coordinate should
        appear. Should be a string containing zero or more of ``b``, ``t``,
        ``l``, ``r``.
        """
        raise NotImplementedError()

    def set_ticks_color(self, color):
        self.ticks.set_color(color)

    def set_ticks_size(self, color):
        self.ticks.set_ticksize(size)

    def set_ticklabel_size(self, size):
        self.ticklabels.set_size(size)

    def set_ticklabel_color(self, color):
        self.ticklabels.set_color(color)

    def _update_grid(self):

        # For 3-d WCS with a correlated third axis, the *proper* way of
        # drawing a grid should be to find the world coordinates of all pixels
        # and drawing contours. What we are doing here assumes that we can
        # define the grid lines with just two of the coordinates (and
        # therefore assumes that the other coordinates are fixed and set to
        # the value in the slice). Here we basically assume that if the WCS
        # had a third axis, it has been abstracted away in the transformation.

        coord_range = self.parent_axes.get_coord_range()

        tick_world_coordinates = self._formatter_locator.locator(*coord_range[self.coord_index])

        paths = []
        for w in tick_world_coordinates:
            if self.coord_index == 0:
                x_world = np.repeat(w, 1000)
                y_world = np.linspace(coord_range[1][0], coord_range[1][1], 1000)
            else:
                x_world = np.linspace(coord_range[0][0], coord_range[0][1], 1000)
                y_world = np.repeat(w, 1000)
            xy_world = np.vstack([x_world, y_world]).transpose()
            paths.append(self._get_gridline(xy_world))

        self.grid_lines.set_paths(paths)


class AngleCoordinateHelper(BaseCoordinateHelper):

    _formatter_locator_class = AngleFormatterLocator

    is_angle = True

    def _update_ticks(self, renderer):

        # Here we should determine the location and rotation of all the ticks.
        # For each axis, we can check the intersections for the specific
        # coordinate and once we have the tick positions, we can use the WCS to
        # determine the rotations.

        # Find the range of coordinates in all directions
        coord_range = self.parent_axes.get_coord_range()

        # First find the ticks we want to show
        tick_world_coordinates = self._formatter_locator.locator(*coord_range[self.coord_index]) % 360.

        # We want to allow non-standard rectangular frames, so we just rely on
        # the parent axes to tell us what the bounding frame is.
        x_pix, y_pix = self.parent_axes._sample_bounding_frame(1000)

        # TODO: the above could be abstracted to work with any line, which
        # could then be re-used for floating axes1

        # Find angle normal to border and inwards
        dx = x_pix - np.roll(x_pix, 1)
        dy = y_pix - np.roll(y_pix, 1)
        normal_angle = np.degrees(np.arctan2(dx, -dy))

        # Transform to world coordinates
        world = self.transform.inverted().transform(np.vstack([x_pix, y_pix]).transpose())[:, self.coord_index]

        # Let's just code up the algorithm with simple loops and we can then
        # see whether to optimize it array-wise, or just cythonize it.

        # We find for each interval the starting and ending coordinate,
        # ensuring that we take wrapping into account correctly.
        w1 = world % 360.
        w2 = np.roll(world, -1) % 360.
        w1[w2 - w1 > 180.] += 360
        w2[w1 - w2 > 180.] += 360

        # Need to check ticks as well as ticks + 360, since the above can
        # produce pairs such as 359 to 361 or 0.5 to 1.5, both of which would
        # match a tick at 0.75.
        check_ticks = np.hstack([tick_world_coordinates, tick_world_coordinates + 360.])

        self.ticks.clear()

        for t in check_ticks:

            # Find steps where a tick is present
            intersections = np.nonzero(((t > w1) & (t < w2)) | ((t < w1) & (t > w2)))[0]

            # Loop over ticks, and find exact pixel coordinates by linear
            # interpolation
            for imin in intersections:
                imax = (imin + 1) % len(world)
                frac = (t - w1[imin]) / (w2[imin] - w1[imin])
                x_pix_i = x_pix[imin] + frac * (x_pix[imax] - x_pix[imin])
                y_pix_i = y_pix[imin] + frac * (y_pix[imax] - y_pix[imin])

                self.ticks.add(pixel=(x_pix_i, y_pix_i),
                               world=t % 360.,
                               angle=normal_angle[imin])

                self.ticklabels.add(pixel=(x_pix_i, y_pix_i),
                                    world=t % 360.,
                                    angle=normal_angle[imin],
                                    text=self._formatter_locator.formatter([t % 360.])[0])

        xscale, yscale = get_pixels_to_data_scales(self.parent_axes)

        self.ticklabels.set_pixel_to_data_scaling(xscale, yscale)

        self._update_grid()

    def _get_gridline(self, xy_world):
        return get_lon_lat_path(self.parent_axes, self.transform, xy_world)


class ScalarCoordinateHelper(BaseCoordinateHelper):

    _formatter_locator_class = ScalarFormatterLocator

    is_angle = False

    def _update_ticks(self, renderer):

        # Here we should determine the location and rotation of all the ticks.
        # For each axis, we can check the intersections for the specific
        # coordinate and once we have the tick positions, we can use the WCS to
        # determine the rotations.

        # Find the range of coordinates in all directions
        coord_range = self.parent_axes.get_coord_range()

        # First find the ticks we want to show
        tick_world_coordinates = self._formatter_locator.locator(*coord_range[self.coord_index])

        # We want to allow non-standard rectangular frames, so we just rely on
        # the parent axes to tell us what the bounding frame is.
        x_pix, y_pix = self.parent_axes._sample_bounding_frame(1000)

        # TODO: the above could be abstracted to work with any line, which
        # could then be re-used for floating axes1

        # Find angle normal to border and inwards
        dx = x_pix - np.roll(x_pix, 1)
        dy = y_pix - np.roll(y_pix, 1)
        normal_angle = np.degrees(np.arctan2(dx, -dy))

        # Transform to world coordinates
        world = self.transform.inverted().transform(np.vstack([x_pix, y_pix]).transpose())[:, self.coord_index]

        # Let's just code up the algorithm with simple loops and we can then
        # see whether to optimize it array-wise, or just cythonize it.

        # We find for each interval the starting and ending coordinate,
        # ensuring that we take wrapping into account correctly.
        w1 = world
        w2 = np.roll(world, -1)

        # Need to check ticks as well as ticks + 360, since the above can
        # produce pairs such as 359 to 361 or 0.5 to 1.5, both of which would
        # match a tick at 0.75.
        check_ticks = tick_world_coordinates

        self.ticks.clear()

        for t in check_ticks:

            # Find steps where a tick is present
            intersections = np.nonzero(((t > w1) & (t < w2)) | ((t < w1) & (t > w2)))[0]

            # Loop over ticks, and find exact pixel coordinates by linear
            # interpolation
            for imin in intersections:
                imax = (imin + 1) % len(world)
                frac = (t - w1[imin]) / (w2[imin] - w1[imin])
                x_pix_i = x_pix[imin] + frac * (x_pix[imax] - x_pix[imin])
                y_pix_i = y_pix[imin] + frac * (y_pix[imax] - y_pix[imin])

                self.ticks.add(pixel=(x_pix_i, y_pix_i),
                               world=t,
                               angle=normal_angle[imin])

                self.ticklabels.add(pixel=(x_pix_i, y_pix_i),
                                    world=t,
                                    angle=normal_angle[imin],
                                    text=self._formatter_locator.formatter([t])[0])

        xscale, yscale = get_pixels_to_data_scales(self.parent_axes)

        self.ticklabels.set_pixel_to_data_scaling(xscale, yscale)

        self._update_grid()

    def _get_gridline(self, xy_world):
        return get_gridline_path(self.parent_axes, self.transform, xy_world)

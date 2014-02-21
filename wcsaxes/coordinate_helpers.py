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

from .formatter_locator import AngleFormatterLocator
from .ticks import Ticks
from .grid_paths import get_lon_lat_path

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

        # Initialize formatter/locator
        self._formatter_locator = self._formatter_locator_class()

        # Initialize container for the grid lines
        self.grid_lines = PathCollection([],
                                         transform=self.parent_axes.transData,
                                         facecolor='none', visible=False,
                                         figure=self.parent_axes.get_figure())

        # Default parameters for tick labels
        self.set_ticklabel_size('medium')
        self.set_ticklabel_color('black')

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
        self._ticklabel_size = size

    def set_ticklabel_color(self, color):
        self._ticklabel_color = color


class SkyCoordinateHelper(BaseCoordinateHelper):

    _formatter_locator_class = AngleFormatterLocator

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
        world = self.transform.inverted().transform(np.vstack([x_pix, y_pix]).transpose())[:,self.coord_index]

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

        # Add labels

        tick_normal = np.array(self.ticks.angle) % 360

        for i in range(len(self.ticks)):

            coordinate = self.ticks.world[i]
            x, y = self.ticks.pixel[i]
            angle = self.ticks.angle[i]
            label_text = self._formatter_locator.formatter([coordinate])[0]

            t = Text(text=label_text, color=self._ticklabel_color, size=self._ticklabel_size, clip_on=False)
            text_size = renderer.points_to_pixels(t.get_size())

            # TODO: there must be a better way to figure out the text size in
            # data units
            axmin, aymin, axmax, aymax = self.parent_axes.get_window_extent().bounds
            xmin, xmax = self.parent_axes.get_xlim()
            ymin, ymax = self.parent_axes.get_ylim()
            xscale = abs((xmax - xmin) / (axmax - axmin))
            yscale = abs((ymax - ymin) / (aymax - aymin))

            pad = text_size * 0.4

            # TODO: do something smarter for arbitrary directions
            if np.abs(tick_normal[i]) < 45.:
                ha = 'right'
                va = 'center'
                dx = -pad
                dy = 0.
            elif np.abs(tick_normal[i] - 90.) < 45:
                ha = 'center'
                va = 'top'
                dx = 0
                dy = -pad
            elif np.abs(tick_normal[i] - 180.) < 45:
                ha = 'left'
                va = 'center'
                dx = pad
                dy = 0
            else:
                ha = 'center'
                va = 'bottom'
                dx = 0
                dy = pad

            dx *= xscale
            dy *= yscale

            t.set_position((x+dx, y+dy))
            t.set_ha(ha)
            t.set_va(va)

            t.set_figure(self.parent_axes.get_figure())
            t.set_transform(self.parent_axes.transData)

            t.draw(renderer)

        self._update_grid()

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
                lon = np.repeat(w, 1000)
                lat = np.linspace(coord_range[1][0], coord_range[1][1], 1000)
            else:
                lon = np.linspace(coord_range[0][0], coord_range[0][1], 1000)
                lat = np.repeat(w, 1000)
            lon_lat = np.vstack([lon, lat]).transpose()
            paths.append(get_lon_lat_path(self.parent_axes, self.transform, lon_lat))

        self.grid_lines.set_paths(paths)




class ScalarCoordinateHelper(BaseCoordinateHelper):
    pass

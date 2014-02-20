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

    @abc.abstractmethod
    def set_major_formatter(self, formatter):
        """
        Set the major formatter for the ticks - should be either a
        ``Formatter`` instance for world coordinates, or a string such as
        ``dd:mm:ss.s``.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def set_ticks(self, spacing=None, number=None):
        """
        Set the spacing/value of the ticks. This can take:

        * A list of tick position
        * A spacing with the ``spacing=`` option (can take quantities)
        * An approximate number of tick marks with the ``number=`` option
        """
        raise NotImplementedError()

    def grid(self):
        """
        Draw grid lines for just this coordinate. Should return a
        ``LineCollection`` instance.
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def locator(self):
        raise NotImplementedError()

    @abc.abstractproperty
    def formatter(self):
        raise NotImplementedError()


class SkyCoordinateHelper(BaseCoordinateHelper):

    def __init__(self, parent_axes=None, transform=None, coord_index=None):

        super(SkyCoordinateHelper, self).__init__()

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
        self._formatter_locator_helper = AngleFormatterLocator()

        self.tick_positions_world = []

        # Initialize container for the grid lines
        self.grid_lines = PathCollection([])
        self.parent_axes.add_collection(self.grid_lines)
        self.grid_lines.set_visible(False)

        self.text_labels = []

    def _update_ticks(self, renderer):

        # Here we should determine the location and rotation of all the ticks.
        # For each axis, we can check the intersections for the specific
        # coordinate and once we have the tick positions, we can use the WCS to
        # determine the rotations.

        coord_range = self.parent_axes.get_coord_range()

        # First find the ticks we want to show
        tick_world_coordinates = self._formatter_locator_helper.locator(*coord_range[self.coord_index])

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

        # Smooth out discontinuities
        for i in range(len(world) - 1):
            if(abs(world[i] - world[i + 1]) > 180.):
                if(world[i] > world[i + 1]):
                    world[i + 1:] = world[i + 1:] + 360.
                else:
                    world[i + 1:] = world[i + 1:] - 360.

        # Search for intersections
        # TODO: can be made more accurate/robust
        tick_world_coordinates_2 = []
        tick_normal = []
        tick_pixel_coordinates = []
        tick_angles = []
        for tick_position in tick_world_coordinates:
            inter = np.where(((world[:-1] <= tick_position) & (world[1:] > tick_position)) |
                             ((world[:-1] > tick_position) & (world[1:] <= tick_position)))[0]
            for i in inter:
                tick_pixel_coordinates.append((x_pix[i], y_pix[i]))
                tick_angles.append(normal_angle[i])
                tick_world_coordinates_2.append(tick_position)
                tick_normal.append(normal_angle[i])

        # And finally we set the locations and angles
        self.ticks._set_positions(world=tick_world_coordinates,
                                 pixel=tick_pixel_coordinates,
                                 angles=tick_angles)

        # Add labels
        for label in self.text_labels:
            label.remove()

        tick_normal = np.array(tick_normal) % 360

        for i in range(len(tick_world_coordinates_2)):

            coordinate = tick_world_coordinates_2[i]
            x, y = tick_pixel_coordinates[i]
            angle = tick_angles[i]
            label_text = self._formatter_locator_helper.formatter([coordinate])[0]

            t = Text(text=label_text, color='green', clip_on=False)
            text_size = renderer.points_to_pixels(t.get_size())

            pad = text_size * 0.4

            # In future, do something smarter for arbitrary directions
            if np.abs(tick_normal[i] - 90.) < 45:
                ha = 'center'
                va = 'bottom'
                dx = 0
                dy = -text_size - pad
            elif np.abs(tick_normal[i] - 180.) < 45:
                ha = 'left'
                va = 'center'
                dx = pad
                dy = 0
            elif np.abs(tick_normal[i] - 270.) < 45:
                ha = 'center'
                va = 'bottom'
                dx = 0
                dy = pad
            else:
                ha = 'right'
                va = 'center'
                dx = -pad
                dy = 0.

            t.set_position((x+dx, y+dy))
            t.set_ha(ha)
            t.set_va(va)

            self.parent_axes.add_artist(t)

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

        paths = []
        for w in self.ticks.get_world_coordinates():
            if self.coord_index == 0:
                lon = np.repeat(w, 1000)
                lat = np.linspace(coord_range[1][0], coord_range[1][1], 1000)
            else:
                lon = np.linspace(coord_range[0][0], coord_range[0][1], 1000)
                lat = np.repeat(w, 1000)
            lon_lat = np.vstack([lon, lat]).transpose()
            paths.append(get_lon_lat_path(self.parent_axes, self.transform, lon_lat))

        self.grid_lines.set_paths(paths)

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
            self._formatter_locator_helper.format = formatter
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
            self._formatter_locator_helper.values = values
        elif spacing is not None:
            self._formatter_locator_helper.spacing = spacing
        elif number is not None:
            self._formatter_locator_helper.number = number
        else:
            raise ValueError("one of values, spacing, or number should be "
                             "specified")

    @property
    def locator(self):
        return _formatter_locator_helper.locator

    @property
    def formatter(self):
        return _formatter_locator_helper.formatter

    def draw(self, renderer):

        renderer.open_group('coordinate_axis')

        self._update_ticks(renderer)
        self.ticks.draw(renderer)

        renderer.close_group('coordinate_axis')


class ScalarCoordinateHelper(BaseCoordinateHelper):
    pass

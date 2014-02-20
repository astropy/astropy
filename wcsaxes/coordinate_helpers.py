import abc

import numpy as np

from matplotlib.ticker import Formatter
from matplotlib.transforms import Affine2D, ScaledTranslation
from matplotlib.collections import PathCollection

from .formatter_locator import AngleFormatterLocator
from . import six
from .ticks import Ticks
from .grid_paths import get_lon_lat_path

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

    def _update_ticks(self, coord_range):

        # Here we should determine the location and rotation of all the ticks.
        # For each axis, we can check the intersections for the specific
        # coordinate and once we have the tick positions, we can use the WCS to
        # determine the rotations.

        # First find the ticks we want to show
        self.tick_positions_world = self._formatter_locator_helper.locator(*coord_range)

        # We want to allow non-standard rectangular frames, so we just rely on
        # the parent axes to tell us what the bounding frame is.
        x_pix, y_pix = self.parent_axes._sample_bounding_frame(1000)

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
        locations = []
        angles = []
        for tick_position in self.tick_positions_world:
            inter = np.where(((world[:-1] <= tick_position) & (world[1:] > tick_position)) |
                             ((world[:-1] > tick_position) & (world[1:] <= tick_position)))[0]
            for i in inter:
                locations.append((x_pix[i], y_pix[i]))
                angles.append(normal_angle[i])

        print(angles)
        # And finally we set the locations and angles
        self.ticks.set_locs_angles(zip(locations, angles))

        self.grid(True)

    def grid(self, draw_grid):
        paths = []
        for w in self.tick_positions_world:
            if self.coord_index == 0:
                lon = np.repeat(w, 1000)
                lat = np.linspace(-89.999, 89.999, 1000)
            else:
                lon = np.linspace(-179.999, 179.999, 1000)
                lat = np.repeat(w, 1000)
            lon_lat = np.vstack([lon, lat]).transpose()
            paths.append(get_lon_lat_path(self.parent_axes, self.transform, lon_lat))
        self.parent_axes.add_collection(PathCollection(paths, edgecolors='b', facecolors='none', alpha=0.4))

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

        # self._update_ticks()
        self.ticks.draw(renderer)

        # if self.grid_lines.get_visible():
        #     self.grid_lines.draw(renderer)

        renderer.close_group('coordinate_axis')


class ScalarCoordinateHelper(BaseCoordinateHelper):
    pass

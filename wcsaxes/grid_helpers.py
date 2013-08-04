import abc

from mpl_toolkits.axisartist import angle_helper, GridHelperCurveLinear

from .coordinate_helpers import SkyCoordinateHelper
from .transforms import WCSWorld2PixelTransform


class BaseGridHelper(object):

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def grid(self):
        """
        Draw grid lines for both coordinates.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def __getitem__(self):
        """
        Access the coordinates by index or by name
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def set_visible(self, visibility):
        """
        Set whether the coordinate system is visible
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def enable_offset_mode(self, reference_coordinates):
        """
        Enable the offset mode given a set of reference cooridnates
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def disable_offset_mode(self):
        """
        Disable the offset mode
        """
        raise NotImplementedError()


class SkyGridHelper(BaseGridHelper):

    def __init__(self, axes, wcs):

        self._axes = axes
        self._wcs = wcs

        # The following code was taken from the demo_floating_axis example in
        # Matplotlib. We make use of helpers in the mpl_toolkits package.

        self._extreme_finder = angle_helper.ExtremeFinderCycle(20, 20,
                                                               lon_cycle=360,
                                                               lat_cycle=None,
                                                               lon_minmax=(-180., 180),
                                                               lat_minmax = (-90., 90.)
                                                               )

        # Set up transform
        self._transform = WCSWorld2PixelTransform(self._wcs)

        # Set up grid helper
        self._grid_helper = GridHelperCurveLinear(self._transform, extreme_finder=self._extreme_finder)

        # Set up coordinates
        self._coords = {}
        self._coords[0] = SkyCoordinateHelper(parent=self, index=1)
        self._coords[1] = SkyCoordinateHelper(parent=self, index=2)

        # Set up aliases for coordinates
        name_1 = self._wcs.wcs.ctype[0][:4]
        self._coords[name_1.lower()] = self._coords[0]
        name_2 = self._wcs.wcs.ctype[1][:4]
        self._coords[name_2.lower()] = self._coords[1]

    @property
    def grid_helper(self):
        return self._grid_helper

    def __getitem__(self, item):
        if isinstance(item, basestring):
            return self._coords[item.lower()]
        else:
            return self._coords[item]

    def grid(self, *args, **kwargs):
        self._axes.grid(*args, **kwargs)

    def set_visible(self, visibility):
        raise NotImplementedError()

    def enable_offset_mode(self, reference_coordinates):
        raise NotImplementedError()

    def disable_offset_mode(self):
        raise NotImplementedError()


class ScalarGridHelper(BaseGridHelper):
    pass

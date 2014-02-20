import abc

from .coordinate_helpers import SkyCoordinateHelper
from .transforms import WCSWorld2PixelTransform
from . import six


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


class SkyCoordinatesMap(object):

    def __init__(self, axes, wcs):

        # Keep track of parent axes and WCS
        self._axes = axes
        self._wcs = wcs

        # Set up transform
        self._transform = WCSWorld2PixelTransform(self._wcs)

        # Set up coordinates
        self._coords = {}
        self._coords[0] = SkyCoordinateHelper(parent_axes=axes,
                                              transform=self._transform, coord_index=0)
        self._coords[1] = SkyCoordinateHelper(parent_axes=axes,
                                              transform=self._transform, coord_index=1)

        # Set up aliases for coordinates
        name_1 = self._wcs.wcs.ctype[0][:4]
        self._coords[name_1.lower()] = self._coords[0]
        name_2 = self._wcs.wcs.ctype[1][:4]
        self._coords[name_2.lower()] = self._coords[1]

    def __getitem__(self, item):
        if isinstance(item, six.string_types):
            return self._coords[item.lower()]
        else:
            return self._coords[item]

    def set_visible(self, visibility):
        raise NotImplementedError()

    def enable_offset_mode(self, reference_coordinates):
        raise NotImplementedError()

    def disable_offset_mode(self):
        raise NotImplementedError()


class ScalarGridHelper(BaseGridHelper):
    pass

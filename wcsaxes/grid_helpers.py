from .coordinate_helpers import AngleCoordinateHelper, ScalarCoordinateHelper
from .transforms import WCSWorld2PixelTransform
from .utils import coord_type_from_ctype

from . import six


class CoordinatesMap(object):

    def __init__(self, axes, wcs):

        # Keep track of parent axes and WCS
        self._axes = axes
        self._wcs = wcs

        # Set up transform
        self._transform = WCSWorld2PixelTransform(self._wcs)

        # Set up coordinates
        self._coords = {}
        for coord_index in [0, 1]:
            coord_type = coord_type_from_ctype(wcs.wcs.ctype[coord_index])
            if coord_type in ['longitude', 'latitude']:
                self._coords[coord_index] = AngleCoordinateHelper(parent_axes=axes,
                                                                  transform=self._transform,
                                                                  coord_index=coord_index,
                                                                  coord_type=coord_type)
            else:
                self._coords[coord_index] = ScalarCoordinateHelper(parent_axes=axes,
                                                                   transform=self._transform,
                                                                   coord_index=coord_index,
                                                                   coord_type=coord_type)

        # Set up aliases for coordinates
        name_1 = self._wcs.wcs.ctype[0][:4].replace('-', '')
        self._coords[name_1.lower()] = self._coords[0]
        name_2 = self._wcs.wcs.ctype[1][:4].replace('-', '')
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

    def grid(self, show_grid=True, **kwrags):
        """
        Plot gridlines for both coordinates.

        Standard matplotlib appearance options (color, alpha, etc.) can be
        passed as keyword arguments.

        Parameters
        ----------
        draw_grid : bool
            Whether to show the gridlines
        """
        self[0].grid(show_grid=show_grid, **kwargs)
        self[1].grid(show_grid=show_grid, **kwargs)

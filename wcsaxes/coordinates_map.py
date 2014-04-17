from .coordinate_helpers import CoordinateHelper
from .transforms import WCSPixel2WorldTransform
from .utils import coord_type_from_ctype
from .frame import RectangularFrame

from . import six


class CoordinatesMap(object):

    def __init__(self, axes, wcs, transform=None, slice=None):

        # Keep track of parent axes and WCS
        self._axes = axes
        self._wcs = wcs

        # Set up transform
        if transform is None:
            self._transform = WCSPixel2WorldTransform(self._wcs, slice=slice)
        else:
            self._transform = transform

        self.frame = RectangularFrame(axes, self._transform)

        # Set up coordinates
        self._coords = {}
        for coord_index in range(self._wcs.wcs.naxis):
            coord_type = coord_type_from_ctype(wcs.wcs.ctype[coord_index])
            self._coords[coord_index] = CoordinateHelper(parent_axes=axes,
                                                         transform=self._transform,
                                                         coord_index=coord_index,
                                                         coord_type=coord_type,
                                                         frame=self.frame)


        # Set up aliases for coordinates
        name_1 = self._wcs.wcs.ctype[0][:4].replace('-', '')
        self._coords[name_1.lower()] = self._coords[0]
        name_2 = self._wcs.wcs.ctype[1][:4].replace('-', '')
        self._coords[name_2.lower()] = self._coords[1]

        # Set up default labels if possible

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

    def __iter__(self):
        for coord in range(self._wcs.wcs.naxis):
            yield self._coords[coord]

    def grid(self, draw_grid=True, grid_type='lines', **kwargs):
        """
        Plot gridlines for both coordinates.

        Standard matplotlib appearance options (color, alpha, etc.) can be
        passed as keyword arguments.

        Parameters
        ----------
        draw_grid : bool
            Whether to show the gridlines
        grid_type : { 'lines' | 'contours' }
            Whether to plot the contours by determining the grid lines in
            world coordinates and then plotting them in world coordinates
            (``'lines'``) or by determining the world coordinates at many
            positions in the image and then drawing contours
            (``'contours'``). The first is recommended for 2-d images, while
            for 3-d (or higher dimensional) cubes, the ``'contours'`` option
            is recommended.
        """
        for coord in self:
            coord.grid(draw_grid=draw_grid, grid_type=grid_type, **kwargs)

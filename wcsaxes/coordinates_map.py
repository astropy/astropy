# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy.extern import six
from .coordinate_helpers import CoordinateHelper
from .transforms import WCSPixel2WorldTransform
from .utils import coord_type_from_ctype
from .frame import RectangularFrame
from .coordinate_range import find_coordinate_range


class CoordinatesMap(object):

    def __init__(self, axes, wcs=None, transform=None, coord_meta=None,
                 slice=None, frame_class=RectangularFrame,
                 previous_frame_path=None):

        # Keep track of parent axes and WCS
        self._axes = axes

        if wcs is None:
            if transform is None:
                raise ValueError("Either `wcs` or `transform` are required")
            if coord_meta is None:
                raise ValueError("`coord_meta` is required when "
                                 "`transform` is passed")
            self._transform = transform
            naxis = 2
        else:
            if transform is not None:
                raise ValueError("Cannot specify both `wcs` and `transform`")
            if coord_meta is not None:
                raise ValueError("Cannot pass `coord_meta` if passing `wcs`")
            self._transform = WCSPixel2WorldTransform(wcs, slice=slice)
            naxis = wcs.wcs.naxis

        self.frame = frame_class(axes, self._transform, path=previous_frame_path)

        # Set up coordinates
        self._coords = []
        self._aliases = {}

        for coord_index in range(naxis):

            # Extract coordinate metadata from WCS object or transform
            if wcs is not None:
                coord_type, coord_wrap = coord_type_from_ctype(wcs.wcs.ctype[coord_index])
                coord_unit = wcs.wcs.cunit[coord_index]
                name = wcs.wcs.ctype[coord_index][:4].replace('-', '')
            else:
                try:
                    coord_type = coord_meta['type'][coord_index]
                    coord_wrap = coord_meta['wrap'][coord_index]
                    coord_unit = coord_meta['unit'][coord_index]
                    name = coord_meta['name'][coord_index]
                except IndexError:
                    raise ValueError("coord_meta items should have a length of {0}".format(len(wcs.wcs.naxis)))

            self._coords.append(CoordinateHelper(parent_axes=axes,
                                                 parent_map=self,
                                                 transform=self._transform,
                                                 coord_index=coord_index,
                                                 coord_type=coord_type,
                                                 coord_wrap=coord_wrap,
                                                 coord_unit=coord_unit,
                                                 frame=self.frame))

            # Set up aliases for coordinates
            self._aliases[name.lower()] = coord_index

    def __getitem__(self, item):
        if isinstance(item, six.string_types):
            return self._coords[self._aliases[item.lower()]]
        else:
            return self._coords[item]

    def set_visible(self, visibility):
        raise NotImplementedError()

    def enable_offset_mode(self, reference_coordinates):
        raise NotImplementedError()

    def disable_offset_mode(self):
        raise NotImplementedError()

    def __iter__(self):
        for coord in self._coords:
            yield coord

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

    def get_coord_range(self):
        xmin, xmax = self._axes.get_xlim()
        ymin, ymax = self._axes.get_ylim()
        return find_coordinate_range(self._transform,
                                     [xmin, xmax, ymin, ymax],
                                     [coord.coord_type for coord in self],
                                     [coord.coord_unit for coord in self])

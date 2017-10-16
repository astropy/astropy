# Licensed under a 3-clause BSD style license - see LICENSE.rst



from .coordinate_helpers import CoordinateHelper
from .transforms import WCSPixel2WorldTransform
from .utils import coord_type_from_ctype
from .frame import RectangularFrame
from .coordinate_range import find_coordinate_range


class CoordinatesMap:
    """
    A container for coordinate helpers that represents a coordinate system.

    This object can be used to access coordinate helpers by index (like a list)
    or by name (like a dictionary).

    Parameters
    ----------
    axes : :class:`~astropy.visualization.wcsaxes.WCSAxes`
        The axes the coordinate map belongs to.
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
    slice : tuple, optional
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
    previous_frame_path : `~matplotlib.path.Path`, optional
        When changing the WCS of the axes, the frame instance will change but
        we might want to keep re-using the same underlying matplotlib
        `~matplotlib.path.Path` - in that case, this can be passed to this
        keyword argument.
    """

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
        if isinstance(item, str):
            return self._coords[self._aliases[item.lower()]]
        else:
            return self._coords[item]

    def set_visible(self, visibility):
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

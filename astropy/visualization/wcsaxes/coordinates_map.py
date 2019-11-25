# Licensed under a 3-clause BSD style license - see LICENSE.rst

from textwrap import indent
from collections import OrderedDict

from .coordinate_helpers import CoordinateHelper
from .frame import RectangularFrame, RectangularFrame1D
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
    transform : `~matplotlib.transforms.Transform`, optional
        The transform for the data.
    coord_meta : dict, optional
        A dictionary providing additional metadata. This should include the keys
        ``type``, ``wrap``, and ``unit``. Each of these should be a list with as
        many items as the dimension of the coordinate system. The ``type``
        entries should be one of ``longitude``, ``latitude``, or ``scalar``, the
        ``wrap`` entries should give, for the longitude, the angle at which the
        coordinate wraps (and `None` otherwise), and the ``unit`` should give
        the unit of the coordinates as :class:`~astropy.units.Unit` instances.
        This can optionally also include a ``format_unit`` entry giving the
        units to use for the tick labels (if not specified, this defaults to
        ``unit``).
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

    def __init__(self, axes, transform=None, coord_meta=None,
                 frame_class=RectangularFrame, previous_frame_path=None):

        self._axes = axes
        self._transform = transform

        self.frame = frame_class(axes, self._transform, path=previous_frame_path)

        # Set up coordinates
        self._coords = []
        self._aliases = {}

        visible_count = 0

        for index in range(len(coord_meta['type'])):

            # Extract coordinate metadata
            coord_type = coord_meta['type'][index]
            coord_wrap = coord_meta['wrap'][index]
            coord_unit = coord_meta['unit'][index]
            name = coord_meta['name'][index]

            visible = True
            if 'visible' in coord_meta:
                visible = coord_meta['visible'][index]

            format_unit = None
            if 'format_unit' in coord_meta:
                format_unit = coord_meta['format_unit'][index]

            default_label = name[0] if isinstance(name, (tuple, list)) else name
            if 'default_axis_label' in coord_meta:
                default_label = coord_meta['default_axis_label'][index]

            coord_index = None
            if visible:
                visible_count += 1
                coord_index = visible_count - 1

            self._coords.append(CoordinateHelper(parent_axes=axes,
                                                 parent_map=self,
                                                 transform=self._transform,
                                                 coord_index=coord_index,
                                                 coord_type=coord_type,
                                                 coord_wrap=coord_wrap,
                                                 coord_unit=coord_unit,
                                                 format_unit=format_unit,
                                                 frame=self.frame,
                                                 default_label=default_label))

            # Set up aliases for coordinates
            if isinstance(name, tuple):
                for nm in name:
                    self._aliases[nm] = index
            else:
                self._aliases[name.lower()] = index

    def __getitem__(self, item):
        if isinstance(item, str):
            return self._coords[self._aliases[item.lower()]]
        else:
            return self._coords[item]

    def __contains__(self, item):
        if isinstance(item, str):
            return item.lower() in self._aliases
        else:
            return 0 <= item < len(self._coords)

    def set_visible(self, visibility):
        raise NotImplementedError()

    def __iter__(self):
        for coord in self._coords:
            yield coord

    def grid(self, draw_grid=True, grid_type=None, **kwargs):
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
            is recommended. By default, 'lines' is used if the transform has
            an inverse, otherwise 'contours' is used.
        """
        for coord in self:
            coord.grid(draw_grid=draw_grid, grid_type=grid_type, **kwargs)

    def get_coord_range(self):
        xmin, xmax = self._axes.get_xlim()

        if isinstance(self.frame, RectangularFrame1D):
            extent = [xmin, xmax]
        else:
            ymin, ymax = self._axes.get_ylim()
            extent = [xmin, xmax, ymin, ymax]

        return find_coordinate_range(self._transform,
                                     extent,
                                     [coord.coord_type for coord in self if coord.coord_index is not None],
                                     [coord.coord_unit for coord in self if coord.coord_index is not None],
                                     [coord.coord_wrap for coord in self if coord.coord_index is not None])

    def _as_table(self):

        # Import Table here to avoid importing the astropy.table package
        # every time astropy.visualization.wcsaxes is imported.
        from astropy.table import Table  # noqa

        rows = []
        for icoord, coord in enumerate(self._coords):
            aliases = [key for key, value in self._aliases.items() if value == icoord]
            row = OrderedDict([('index', icoord), ('aliases', ' '.join(aliases)),
                               ('type', coord.coord_type), ('unit', coord.coord_unit),
                               ('wrap', coord.coord_wrap), ('format_unit', coord.get_format_unit()),
                               ('visible', 'no' if coord.coord_index is None else 'yes')])
            rows.append(row)
        return Table(rows=rows)

    def __repr__(self):
        s = f'<CoordinatesMap with {len(self._coords)} world coordinates:\n\n'
        table = indent(str(self._as_table()), '  ')
        return s + table + '\n\n>'

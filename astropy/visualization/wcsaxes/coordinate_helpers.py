# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function, division, absolute_import

"""
This file defines the classes used to represent a 'coordinate', which includes
axes, ticks, tick labels, and grid lines.
"""

import numpy as np

from matplotlib.ticker import Formatter
from matplotlib.transforms import Affine2D, ScaledTranslation
from matplotlib.patches import PathPatch
from matplotlib import rcParams

from ... import units as u
from ...extern import six

from .formatter_locator import AngleFormatterLocator, ScalarFormatterLocator
from .ticks import Ticks
from .ticklabels import TickLabels
from .axislabels import AxisLabels
from .grid_paths import get_lon_lat_path, get_gridline_path

__all__ = ['CoordinateHelper']


def wrap_angle_at(values, coord_wrap):
    # On ARM processors, np.mod emits warnings if there are NaN values in the
    # array, although this doesn't seem to happen on other processors.
    with np.errstate(invalid='ignore'):
        return np.mod(values - coord_wrap, 360.) - (360. - coord_wrap)


class CoordinateHelper(object):
    """
    Helper class to control one of the coordinates in the
    :class:`~astropy.visualization.wcsaxes.WCSAxes`.

    Parameters
    ----------
    parent_axes : :class:`~astropy.visualization.wcsaxes.WCSAxes`
        The axes the coordinate helper belongs to.
    parent_map : :class:`~astropy.visualization.wcsaxes.CoordinatesMap`
        The :class:`~astropy.visualization.wcsaxes.CoordinatesMap` object this
        coordinate belongs to.
    transform : `~matplotlib.transforms.Transform`
        The transform corresponding to this coordinate system.
    coord_index : int
        The index of this coordinate in the
        :class:`~astropy.visualization.wcsaxes.CoordinatesMap`.
    coord_type : {'longitude', 'latitude', 'scalar'}
        The type of this coordinate, which is used to determine the wrapping and
        boundary behavior of coordinates. Longitudes wrap at ``coord_wrap``,
        latitudes have to be in the range -90 to 90, and scalars are unbounded
        and do not wrap.
    coord_unit : `~astropy.units.Unit`
        The unit that this coordinate is in given the output of transform.
    coord_wrap : float
        The angle at which the longitude wraps (defaults to 360)
    frame : `~astropy.visualization.wcsaxes.frame.BaseFrame`
        The frame of the :class:`~astropy.visualization.wcsaxes.WCSAxes`.
    """

    def __init__(self, parent_axes=None, parent_map=None, transform=None, coord_index=None,
                 coord_type='scalar', coord_unit=None, coord_wrap=None, frame=None):

        # Keep a reference to the parent axes and the transform
        self.parent_axes = parent_axes
        self.parent_map = parent_map
        self.transform = transform
        self.coord_index = coord_index
        self.coord_unit = coord_unit
        self.frame = frame

        self.set_coord_type(coord_type, coord_wrap)

        # Initialize ticks
        self.dpi_transform = Affine2D()
        self.offset_transform = ScaledTranslation(0, 0, self.dpi_transform)
        self.ticks = Ticks(transform=parent_axes.transData + self.offset_transform)

        # Initialize tick labels
        self.ticklabels = TickLabels(self.frame,
                                     transform=None,  # display coordinates
                                     figure=parent_axes.get_figure())
        self.ticks.display_minor_ticks(False)
        self.minor_frequency = 5

        # Initialize axis labels
        self.axislabels = AxisLabels(self.frame,
                                     transform=None,  # display coordinates
                                     figure=parent_axes.get_figure())

        # Initialize container for the grid lines
        self.grid_lines = []

        # Initialize grid style. Take defaults from matplotlib.rcParams.
        # Based on matplotlib.axis.YTick._get_gridline.
        #
        # Matplotlib's gridlines use Line2D, but ours use PathPatch.
        # Patches take a slightly different format of linestyle argument.
        lines_to_patches_linestyle = {'-': 'solid',
                                      '--': 'dashed',
                                      '-.': 'dashdot',
                                      ':': 'dotted',
                                      'none': 'none',
                                      'None': 'none',
                                      ' ': 'none',
                                      '': 'none'}
        self.grid_lines_kwargs = {'visible': False,
                                  'facecolor': 'none',
                                  'edgecolor': rcParams['grid.color'],
                                  'linestyle': lines_to_patches_linestyle[rcParams['grid.linestyle']],
                                  'linewidth': rcParams['grid.linewidth'],
                                  'alpha': rcParams.get('grid.alpha', 1.0),
                                  'transform': self.parent_axes.transData}

    def grid(self, draw_grid=True, grid_type='lines', **kwargs):
        """
        Plot grid lines for this coordinate.

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

        if grid_type in ('lines', 'contours'):
            self._grid_type = grid_type
        else:
            raise ValueError("grid_type should be 'lines' or 'contours'")

        if 'color' in kwargs:
            kwargs['edgecolor'] = kwargs.pop('color')

        self.grid_lines_kwargs.update(kwargs)

        if self.grid_lines_kwargs['visible']:
            if not draw_grid:
                self.grid_lines_kwargs['visible'] = False
        else:
            self.grid_lines_kwargs['visible'] = True

    def set_coord_type(self, coord_type, coord_wrap=None):
        """
        Set the coordinate type for the axis.

        Parameters
        ----------
        coord_type : str
            One of 'longitude', 'latitude' or 'scalar'
        coord_wrap : float, optional
            The value to wrap at for angular coordinates
        """

        self.coord_type = coord_type

        if coord_type == 'longitude' and coord_wrap is None:
            self.coord_wrap = 360
        elif coord_type != 'longitude' and coord_wrap is not None:
            raise NotImplementedError('coord_wrap is not yet supported '
                                      'for non-longitude coordinates')
        else:
            self.coord_wrap = coord_wrap

        # Initialize tick formatter/locator
        if coord_type == 'scalar':
            self._coord_unit_scale = None
            self._formatter_locator = ScalarFormatterLocator(unit=self.coord_unit)
        elif coord_type in ['longitude', 'latitude']:
            if self.coord_unit is u.deg:
                self._coord_unit_scale = None
            else:
                self._coord_unit_scale = self.coord_unit.to(u.deg)
            self._formatter_locator = AngleFormatterLocator()
        else:
            raise ValueError("coord_type should be one of 'scalar', 'longitude', or 'latitude'")

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

    def format_coord(self, value):
        """
        Given the value of a coordinate, will format it according to the
        format of the formatter_locator.
        """

        if not hasattr(self, "_fl_spacing"):
            return ""  # _update_ticks has not been called yet

        fl = self._formatter_locator
        if isinstance(fl, AngleFormatterLocator):

            # Convert to degrees if needed
            if self._coord_unit_scale is not None:
                value *= self._coord_unit_scale

            if self.coord_type == 'longitude':
                value = wrap_angle_at(value, self.coord_wrap)
            value = value * u.degree
            value = value.to_value(fl._unit)

        spacing = self._fl_spacing
        string = fl.formatter(values=[value] * fl._unit, spacing=spacing)

        return string[0]

    def set_separator(self, separator):
        """
        Set the separator to use for the angle major tick labels.

        Parameters
        ----------
        separator : The separator between numbers in sexagesimal
        representation. Can be either a string or a tuple.
        """
        if not (self._formatter_locator.__class__ == AngleFormatterLocator):
            raise TypeError("Separator can only be specified for angle coordinates")
        if isinstance(separator, six.string_types) or isinstance(separator, tuple):
            self._formatter_locator.sep = separator
        else:
            raise TypeError("separator should be a string or a tuple")

    def set_format_unit(self, unit):
        """
        Set the unit for the major tick labels.

        Parameters
        ----------
        unit : class:`~astropy.units.Unit`
            The unit to which the tick labels should be converted to.
        """
        if not issubclass(unit.__class__, u.UnitBase):
            raise TypeError("unit should be an astropy UnitBase subclass")
        self._formatter_locator.format_unit = unit

    def set_ticks(self, values=None, spacing=None, number=None, size=None,
                  width=None, color=None, alpha=None, exclude_overlapping=False):
        """
        Set the location and properties of the ticks.

        At most one of the options from ``values``, ``spacing``, or
        ``number`` can be specified.

        Parameters
        ----------
        values : iterable, optional
            The coordinate values at which to show the ticks.
        spacing : float, optional
            The spacing between ticks.
        number : float, optional
            The approximate number of ticks shown.
        size : float, optional
            The length of the ticks in points
        color : str or tuple
            A valid Matplotlib color for the ticks
        exclude_overlapping : bool, optional
            Whether to exclude tick labels that overlap over each other.
        """

        if sum([values is None, spacing is None, number is None]) < 2:
            raise ValueError("At most one of values, spacing, or number should "
                             "be specified")

        if values is not None:
            self._formatter_locator.values = values
        elif spacing is not None:
            self._formatter_locator.spacing = spacing
        elif number is not None:
            self._formatter_locator.number = number

        if size is not None:
            self.ticks.set_ticksize(size)

        if width is not None:
            self.ticks.set_linewidth(width)

        if color is not None:
            self.ticks.set_color(color)

        if alpha is not None:
            self.ticks.set_alpha(alpha)

        self.ticklabels.set_exclude_overlapping(exclude_overlapping)

    def set_ticks_position(self, position):
        """
        Set where ticks should appear

        Parameters
        ----------
        position : str
            The axes on which the ticks for this coordinate should appear.
            Should be a string containing zero or more of ``'b'``, ``'t'``,
            ``'l'``, ``'r'``. For example, ``'lb'`` will lead the ticks to be
            shown on the left and bottom axis.
        """
        self.ticks.set_visible_axes(position)

    def set_ticks_visible(self, visible):
        """
        Set whether ticks are visible or not.

        Parameters
        ----------
        visible : bool
            The visibility of ticks. Setting as ``False`` will hide ticks
            along this coordinate.
        """
        self.ticks.set_visible(visible)

    def set_ticklabel(self, **kwargs):
        """
        Set the visual properties for the tick labels.

        Parameters
        ----------
        kwargs
            Keyword arguments are passed to :class:`matplotlib.text.Text`. These
            can include keywords to set the ``color``, ``size``, ``weight``, and
            other text properties.
        """
        self.ticklabels.set(**kwargs)

    def set_ticklabel_position(self, position):
        """
        Set where tick labels should appear

        Parameters
        ----------
        position : str
            The axes on which the tick labels for this coordinate should
            appear. Should be a string containing zero or more of ``'b'``,
            ``'t'``, ``'l'``, ``'r'``. For example, ``'lb'`` will lead the
            tick labels to be shown on the left and bottom axis.
        """
        self.ticklabels.set_visible_axes(position)

    def set_ticklabel_visible(self, visible):
        """
        Set whether the tick labels are visible or not.

        Parameters
        ----------
        visible : bool
            The visibility of ticks. Setting as ``False`` will hide this
            coordinate's tick labels.
        """
        self.ticklabels.set_visible(visible)

    def set_axislabel(self, text, minpad=1, **kwargs):
        """
        Set the text and optionally visual properties for the axis label.

        Parameters
        ----------
        text : str
            The axis label text.
        minpad : float, optional
            The padding for the label in terms of axis label font size.
        kwargs
            Keywords are passed to :class:`matplotlib.text.Text`. These
            can include keywords to set the ``color``, ``size``, ``weight``, and
            other text properties.
        """

        fontdict = kwargs.pop('fontdict', None)

        # NOTE: When using plt.xlabel/plt.ylabel, minpad can get set explicitly
        # to None so we need to make sure that in that case we change to a
        # default numerical value.
        if minpad is None:
            minpad = 1

        self.axislabels.set_text(text)
        self.axislabels.set_minpad(minpad)
        self.axislabels.set(**kwargs)

        if fontdict is not None:
            self.axislabels.update(fontdict)

    def get_axislabel(self):
        """
        Get the text for the axis label

        Returns
        -------
        label : str
            The axis label
        """
        return self.axislabels.get_text()

    def set_axislabel_position(self, position):
        """
        Set where axis labels should appear

        Parameters
        ----------
        position : str
            The axes on which the axis label for this coordinate should
            appear. Should be a string containing zero or more of ``'b'``,
            ``'t'``, ``'l'``, ``'r'``. For example, ``'lb'`` will lead the
            axis label to be shown on the left and bottom axis.
        """
        self.axislabels.set_visible_axes(position)

    @property
    def locator(self):
        return self._formatter_locator.locator

    @property
    def formatter(self):
        return self._formatter_locator.formatter

    def _draw_grid(self, renderer):

        renderer.open_group('grid lines')

        self._update_ticks()

        if self.grid_lines_kwargs['visible']:

            if self._grid_type == 'lines':
                self._update_grid_lines()
            else:
                self._update_grid_contour()

            if self._grid_type == 'lines':

                frame_patch = self.frame.patch
                for path in self.grid_lines:
                    p = PathPatch(path, **self.grid_lines_kwargs)
                    p.set_clip_path(frame_patch)
                    p.draw(renderer)

            elif self._grid is not None:

                for line in self._grid.collections:
                    line.set(**self.grid_lines_kwargs)
                    line.draw(renderer)

        renderer.close_group('grid lines')

    def _draw_ticks(self, renderer, bboxes, ticklabels_bbox):

        renderer.open_group('ticks')

        self.ticks.draw(renderer)
        self.ticklabels.draw(renderer, bboxes=bboxes,
                             ticklabels_bbox=ticklabels_bbox)

        renderer.close_group('ticks')

    def _draw_axislabels(self, renderer, bboxes, ticklabels_bbox, visible_ticks):

        renderer.open_group('axis labels')

        self.axislabels.draw(renderer, bboxes=bboxes,
                             ticklabels_bbox_list=ticklabels_bbox,
                             visible_ticks=visible_ticks)

        renderer.close_group('axis labels')

    def _update_ticks(self):

        # TODO: this method should be optimized for speed

        # Here we determine the location and rotation of all the ticks. For
        # each axis, we can check the intersections for the specific
        # coordinate and once we have the tick positions, we can use the WCS
        # to determine the rotations.

        # Find the range of coordinates in all directions
        coord_range = self.parent_map.get_coord_range()

        # First find the ticks we want to show
        tick_world_coordinates, self._fl_spacing = self.locator(*coord_range[self.coord_index])
        if self.ticks.get_display_minor_ticks():
            minor_ticks_w_coordinates = self._formatter_locator.minor_locator(self._fl_spacing, self.get_minor_frequency(), *coord_range[self.coord_index])

        # We want to allow non-standard rectangular frames, so we just rely on
        # the parent axes to tell us what the bounding frame is.
        from . import conf
        frame = self.frame.sample(conf.frame_boundary_samples)

        self.ticks.clear()
        self.ticklabels.clear()
        self.lblinfo = []
        self.lbl_world = []
        # Look up parent axes' transform from data to figure coordinates.
        #
        # See:
        # http://matplotlib.org/users/transforms_tutorial.html#the-transformation-pipeline
        transData = self.parent_axes.transData
        invertedTransLimits = transData.inverted()

        for axis, spine in six.iteritems(frame):

            # Determine tick rotation in display coordinates and compare to
            # the normal angle in display coordinates.

            pixel0 = spine.data
            world0 = spine.world[:, self.coord_index]
            world0 = self.transform.transform(pixel0)[:, self.coord_index]
            axes0 = transData.transform(pixel0)

            # Advance 2 pixels in figure coordinates
            pixel1 = axes0.copy()
            pixel1[:, 0] += 2.0
            pixel1 = invertedTransLimits.transform(pixel1)
            world1 = self.transform.transform(pixel1)[:, self.coord_index]

            # Advance 2 pixels in figure coordinates
            pixel2 = axes0.copy()
            pixel2[:, 1] += 2.0 if self.frame.origin == 'lower' else -2.0
            pixel2 = invertedTransLimits.transform(pixel2)
            world2 = self.transform.transform(pixel2)[:, self.coord_index]

            dx = (world1 - world0)
            dy = (world2 - world0)

            # Rotate by 90 degrees
            dx, dy = -dy, dx

            if self._coord_unit_scale is not None:
                dx *= self._coord_unit_scale
                dy *= self._coord_unit_scale

            if self.coord_type == 'longitude':
                # Here we wrap at 180 not self.coord_wrap since we want to
                # always ensure abs(dx) < 180 and abs(dy) < 180
                dx = wrap_angle_at(dx, 180.)
                dy = wrap_angle_at(dy, 180.)

            tick_angle = np.degrees(np.arctan2(dy, dx))

            normal_angle_full = np.hstack([spine.normal_angle, spine.normal_angle[-1]])
            with np.errstate(invalid='ignore'):
                reset = (((normal_angle_full - tick_angle) % 360 > 90.) &
                         ((tick_angle - normal_angle_full) % 360 > 90.))
            tick_angle[reset] -= 180.

            # We find for each interval the starting and ending coordinate,
            # ensuring that we take wrapping into account correctly for
            # longitudes.
            w1 = spine.world[:-1, self.coord_index]
            w2 = spine.world[1:, self.coord_index]

            if self._coord_unit_scale is not None:
                w1 = w1 * self._coord_unit_scale
                w2 = w2 * self._coord_unit_scale

            if self.coord_type == 'longitude':
                w1 = wrap_angle_at(w1, self.coord_wrap)
                w2 = wrap_angle_at(w2, self.coord_wrap)
                with np.errstate(invalid='ignore'):
                    w1[w2 - w1 > 180.] += 360
                    w2[w1 - w2 > 180.] += 360

            # For longitudes, we need to check ticks as well as ticks + 360,
            # since the above can produce pairs such as 359 to 361 or 0.5 to
            # 1.5, both of which would match a tick at 0.75. Otherwise we just
            # check the ticks determined above.
            self._compute_ticks(tick_world_coordinates, spine, axis, w1, w2, tick_angle)

            if self.ticks.get_display_minor_ticks():
                self._compute_ticks(minor_ticks_w_coordinates, spine, axis, w1,
                                    w2, tick_angle, ticks='minor')

        # format tick labels, add to scene
        text = self.formatter(self.lbl_world * tick_world_coordinates.unit, spacing=self._fl_spacing)
        for kwargs, txt in zip(self.lblinfo, text):
            self.ticklabels.add(text=txt, **kwargs)

    def _compute_ticks(self, tick_world_coordinates, spine, axis, w1, w2, tick_angle, ticks='major'):
        tick_world_coordinates_values = tick_world_coordinates.value
        if self.coord_type == 'longitude':
            tick_world_coordinates_values = np.hstack([tick_world_coordinates_values,
                                                       tick_world_coordinates_values + 360])

        for t in tick_world_coordinates_values:

            # Find steps where a tick is present. We have to check
            # separately for the case where the tick falls exactly on the
            # frame points, otherwise we'll get two matches, one for w1 and
            # one for w2.
            with np.errstate(invalid='ignore'):
                intersections = np.hstack([np.nonzero((t - w1) == 0)[0],
                                           np.nonzero(((t - w1) * (t - w2)) < 0)[0]])

            # But we also need to check for intersection with the last w2
            if t - w2[-1] == 0:
                intersections = np.append(intersections, len(w2) - 1)

            # Loop over ticks, and find exact pixel coordinates by linear
            # interpolation
            for imin in intersections:

                imax = imin + 1

                if np.allclose(w1[imin], w2[imin], rtol=1.e-13, atol=1.e-13):
                    continue  # tick is exactly aligned with frame
                else:
                    frac = (t - w1[imin]) / (w2[imin] - w1[imin])
                    x_data_i = spine.data[imin, 0] + frac * (spine.data[imax, 0] - spine.data[imin, 0])
                    y_data_i = spine.data[imin, 1] + frac * (spine.data[imax, 1] - spine.data[imin, 1])
                    x_pix_i = spine.pixel[imin, 0] + frac * (spine.pixel[imax, 0] - spine.pixel[imin, 0])
                    y_pix_i = spine.pixel[imin, 1] + frac * (spine.pixel[imax, 1] - spine.pixel[imin, 1])
                    delta_angle = tick_angle[imax] - tick_angle[imin]
                    if delta_angle > 180.:
                        delta_angle -= 360.
                    elif delta_angle < -180.:
                        delta_angle += 360.
                    angle_i = tick_angle[imin] + frac * delta_angle

                if self.coord_type == 'longitude':
                    world = wrap_angle_at(t, self.coord_wrap)
                else:
                    world = t

                if ticks == 'major':

                    self.ticks.add(axis=axis,
                                   pixel=(x_data_i, y_data_i),
                                   world=world,
                                   angle=angle_i,
                                   axis_displacement=imin + frac)

                    # store information to pass to ticklabels.add
                    # it's faster to format many ticklabels at once outside
                    # of the loop
                    self.lblinfo.append(dict(axis=axis,
                                             pixel=(x_pix_i, y_pix_i),
                                             world=world,
                                             angle=spine.normal_angle[imin],
                                             axis_displacement=imin + frac))
                    self.lbl_world.append(world)

                else:
                    self.ticks.add_minor(minor_axis=axis,
                                         minor_pixel=(x_data_i, y_data_i),
                                         minor_world=world,
                                         minor_angle=angle_i,
                                         minor_axis_displacement=imin + frac)

    def display_minor_ticks(self, display_minor_ticks):
        """
        Display minor ticks for this coordinate.

        Parameters
        ----------
        display_minor_ticks : bool
            Whether or not to display minor ticks.
        """
        self.ticks.display_minor_ticks(display_minor_ticks)

    def get_minor_frequency(self):
        return self.minor_frequency

    def set_minor_frequency(self, frequency):
        """
        Set the frequency of minor ticks per major ticks.

        Parameters
        ----------
        frequency : int
            The number of minor ticks per major ticks.
        """
        self.minor_frequency = frequency

    def _update_grid_lines(self):

        # For 3-d WCS with a correlated third axis, the *proper* way of
        # drawing a grid should be to find the world coordinates of all pixels
        # and drawing contours. What we are doing here assumes that we can
        # define the grid lines with just two of the coordinates (and
        # therefore assumes that the other coordinates are fixed and set to
        # the value in the slice). Here we basically assume that if the WCS
        # had a third axis, it has been abstracted away in the transformation.

        coord_range = self.parent_map.get_coord_range()

        tick_world_coordinates, spacing = self.locator(*coord_range[self.coord_index])
        tick_world_coordinates_values = tick_world_coordinates.value

        n_coord = len(tick_world_coordinates_values)

        from . import conf
        n_samples = conf.grid_samples

        xy_world = np.zeros((n_samples * n_coord, 2))

        self.grid_lines = []
        for iw, w in enumerate(tick_world_coordinates_values):
            subset = slice(iw * n_samples, (iw + 1) * n_samples)
            if self.coord_index == 0:
                xy_world[subset, 0] = np.repeat(w, n_samples)
                xy_world[subset, 1] = np.linspace(coord_range[1][0], coord_range[1][1], n_samples)
            else:
                xy_world[subset, 0] = np.linspace(coord_range[0][0], coord_range[0][1], n_samples)
                xy_world[subset, 1] = np.repeat(w, n_samples)

        # We now convert all the world coordinates to pixel coordinates in a
        # single go rather than doing this in the gridline to path conversion
        # to fully benefit from vectorized coordinate transformations.

        # Currently xy_world is in deg, but transform function needs it in
        # native units
        if self._coord_unit_scale is not None:
            xy_world /= self._coord_unit_scale

        # Transform line to pixel coordinates
        pixel = self.transform.inverted().transform(xy_world)

        # Create round-tripped values for checking
        xy_world_round = self.transform.transform(pixel)

        for iw in range(n_coord):
            subset = slice(iw * n_samples, (iw + 1) * n_samples)
            self.grid_lines.append(self._get_gridline(xy_world[subset], pixel[subset], xy_world_round[subset]))

    def _get_gridline(self, xy_world, pixel, xy_world_round):
        if self.coord_type == 'scalar':
            return get_gridline_path(xy_world, pixel)
        else:
            return get_lon_lat_path(xy_world, pixel, xy_world_round)

    def _update_grid_contour(self):

        if hasattr(self, '_grid') and self._grid:
            for line in self._grid.collections:
                line.remove()

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        from . import conf
        res = conf.contour_grid_samples

        x, y = np.meshgrid(np.linspace(xmin, xmax, res),
                           np.linspace(ymin, ymax, res))
        pixel = np.array([x.ravel(), y.ravel()]).T
        world = self.transform.transform(pixel)
        field = world[:, self.coord_index].reshape(res, res).T

        coord_range = self.parent_map.get_coord_range()

        tick_world_coordinates, spacing = self.locator(*coord_range[self.coord_index])

        # tick_world_coordinates is a Quantities array and we only needs its values
        tick_world_coordinates_values = tick_world_coordinates.value

        if self.coord_type == 'longitude':

            # Find biggest gap in tick_world_coordinates and wrap in middle
            # For now just assume spacing is equal, so any mid-point will do
            mid = 0.5 * (tick_world_coordinates_values[0] + tick_world_coordinates_values[1])
            field = wrap_angle_at(field, mid)
            tick_world_coordinates_values = wrap_angle_at(tick_world_coordinates_values, mid)

            # Replace wraps by NaN
            reset = (np.abs(np.diff(field[:, :-1], axis=0)) > 180) | (np.abs(np.diff(field[:-1, :], axis=1)) > 180)
            field[:-1, :-1][reset] = np.nan
            field[1:, :-1][reset] = np.nan
            field[:-1, 1:][reset] = np.nan
            field[1:, 1:][reset] = np.nan

        if len(tick_world_coordinates_values) > 0:
            self._grid = self.parent_axes.contour(x, y, field.transpose(), levels=np.sort(tick_world_coordinates_values))
        else:
            self._grid = None

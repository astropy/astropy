# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This file defines the classes used to represent a 'coordinate', which includes
axes, ticks, tick labels, and grid lines.
"""

import warnings

import numpy as np
from matplotlib import rcParams
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.transforms import Affine2D, ScaledTranslation

from astropy import units as u
from astropy.utils.exceptions import AstropyDeprecationWarning

from .axislabels import AxisLabels
from .formatter_locator import AngleFormatterLocator, ScalarFormatterLocator
from .frame import EllipticalFrame, RectangularFrame1D
from .grid_paths import get_gridline_path, get_lon_lat_path
from .ticklabels import TickLabels
from .ticks import Ticks
from .utils import MATPLOTLIB_LT_3_8

__all__ = ["CoordinateHelper"]

# Matplotlib's gridlines use Line2D, but ours use PathPatch.
# Patches take a slightly different format of linestyle argument.
LINES_TO_PATCHES_LINESTYLE = {
    "-": "solid",
    "--": "dashed",
    "-.": "dashdot",
    ":": "dotted",
    "none": "none",
    "None": "none",
    " ": "none",
    "": "none",
}


def wrap_angle_at(values, coord_wrap):
    # On ARM processors, np.mod emits warnings if there are NaN values in the
    # array, although this doesn't seem to happen on other processors.
    with np.errstate(invalid="ignore"):
        return np.mod(values - coord_wrap, 360.0) - (360.0 - coord_wrap)


class CoordinateHelper:
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
    format_unit : `~astropy.units.Unit`, optional
        The unit to use to display the coordinates.
    coord_wrap : `astropy.units.Quantity`
        The angle at which the longitude wraps (defaults to 360 degrees).
    frame : `~astropy.visualization.wcsaxes.frame.BaseFrame`
        The frame of the :class:`~astropy.visualization.wcsaxes.WCSAxes`.
    """

    def __init__(
        self,
        parent_axes=None,
        parent_map=None,
        transform=None,
        coord_index=None,
        coord_type="scalar",
        coord_unit=None,
        coord_wrap=None,
        frame=None,
        format_unit=None,
        default_label=None,
    ):
        # Keep a reference to the parent axes and the transform
        self.parent_axes = parent_axes
        self.parent_map = parent_map
        self.transform = transform
        self.coord_index = coord_index
        self.coord_unit = coord_unit
        self._format_unit = format_unit
        self.frame = frame
        self.default_label = default_label or ""
        self._auto_axislabel = True
        self._axislabel_set = False

        # Disable auto label for elliptical frames as it puts labels in
        # annoying places.
        if issubclass(self.parent_axes.frame_class, EllipticalFrame):
            self._auto_axislabel = False

        self.set_coord_type(coord_type, coord_wrap)

        # Initialize ticks
        self.dpi_transform = Affine2D()
        self.offset_transform = ScaledTranslation(0, 0, self.dpi_transform)
        self.ticks = Ticks(transform=parent_axes.transData + self.offset_transform)

        # Initialize tick labels
        self.ticklabels = TickLabels(
            self.frame,
            transform=None,  # display coordinates
            figure=parent_axes.get_figure(),
        )
        self.ticks.display_minor_ticks(rcParams["xtick.minor.visible"])
        self.minor_frequency = 5

        # Initialize axis labels
        self.axislabels = AxisLabels(
            self.frame,
            transform=None,  # display coordinates
            figure=parent_axes.get_figure(),
        )

        # Initialize container for the grid lines
        self.grid_lines = []

        # Initialize grid style. Take defaults from matplotlib.rcParams.
        # Based on matplotlib.axis.YTick._get_gridline.
        self.grid_lines_kwargs = {
            "visible": False,
            "facecolor": "none",
            "edgecolor": rcParams["grid.color"],
            "linestyle": LINES_TO_PATCHES_LINESTYLE[rcParams["grid.linestyle"]],
            "linewidth": rcParams["grid.linewidth"],
            "alpha": rcParams["grid.alpha"],
            "transform": self.parent_axes.transData,
        }

    def grid(self, draw_grid=True, grid_type=None, **kwargs):
        """
        Plot grid lines for this coordinate.

        Standard matplotlib appearance options (color, alpha, etc.) can be
        passed as keyword arguments.

        Parameters
        ----------
        draw_grid : bool
            Whether to show the gridlines
        grid_type : {'lines', 'contours'}
            Whether to plot the contours by determining the grid lines in
            world coordinates and then plotting them in world coordinates
            (``'lines'``) or by determining the world coordinates at many
            positions in the image and then drawing contours
            (``'contours'``). The first is recommended for 2-d images, while
            for 3-d (or higher dimensional) cubes, the ``'contours'`` option
            is recommended. By default, 'lines' is used if the transform has
            an inverse, otherwise 'contours' is used.
        """
        if grid_type == "lines" and not self.transform.has_inverse:
            raise ValueError(
                "The specified transform has no inverse, so the "
                "grid cannot be drawn using grid_type='lines'"
            )

        if grid_type is None:
            grid_type = "lines" if self.transform.has_inverse else "contours"

        if grid_type in ("lines", "contours"):
            self._grid_type = grid_type
        else:
            raise ValueError("grid_type should be 'lines' or 'contours'")

        if "color" in kwargs:
            kwargs["edgecolor"] = kwargs.pop("color")

        self.grid_lines_kwargs.update(kwargs)

        if draw_grid is None:
            draw_grid = True
        self.grid_lines_kwargs["visible"] = draw_grid

    def set_coord_type(self, coord_type, coord_wrap=None):
        """
        Set the coordinate type for the axis.

        Parameters
        ----------
        coord_type : str
            One of 'longitude', 'latitude' or 'scalar'
        coord_wrap : `~astropy.units.Quantity`, optional
            The value to wrap at for angular coordinates.
        """
        self.coord_type = coord_type

        if coord_wrap is not None and not isinstance(coord_wrap, u.Quantity):
            warnings.warn(
                "Passing 'coord_wrap' as a number is deprecated. Use a Quantity with units convertible to angular degrees instead.",
                AstropyDeprecationWarning,
            )
            coord_wrap = coord_wrap * u.deg

        if coord_type == "longitude" and coord_wrap is None:
            self.coord_wrap = 360 * u.deg
        elif coord_type != "longitude" and coord_wrap is not None:
            raise NotImplementedError(
                "coord_wrap is not yet supported for non-longitude coordinates"
            )
        else:
            self.coord_wrap = coord_wrap

        # Initialize tick formatter/locator
        if coord_type == "scalar":
            self._coord_scale_to_deg = None
            self._formatter_locator = ScalarFormatterLocator(unit=self.coord_unit)
        elif coord_type in ["longitude", "latitude"]:
            if self.coord_unit is u.deg:
                self._coord_scale_to_deg = None
            else:
                self._coord_scale_to_deg = self.coord_unit.to(u.deg)
            self._formatter_locator = AngleFormatterLocator(
                unit=self.coord_unit, format_unit=self._format_unit
            )
        else:
            raise ValueError(
                "coord_type should be one of 'scalar', 'longitude', or 'latitude'"
            )

    def set_major_formatter(self, formatter):
        """
        Set the format string to use for the major tick labels.

        Parameters
        ----------
        formatter : str
            The format string to use.
        """
        # TODO: Figure out how to accept a Formatter instance
        if isinstance(formatter, str):
            self._formatter_locator.format = formatter
        else:
            raise TypeError("formatter should be a string")

    def format_coord(self, value, format="auto"):
        """
        Given the value of a coordinate, will format it according to the
        format of the formatter_locator.

        Parameters
        ----------
        value : float
            The value to format
        format : {'auto', 'ascii', 'latex'}, optional
            The format to use - by default the formatting will be adjusted
            depending on whether Matplotlib is using LaTeX or MathTex. To
            get plain ASCII strings, use format='ascii'.
        """
        if not hasattr(self, "_fl_spacing"):
            return ""  # _update_ticks has not been called yet

        fl = self._formatter_locator
        if isinstance(fl, AngleFormatterLocator):
            # Convert to degrees if needed
            if self._coord_scale_to_deg is not None:
                value *= self._coord_scale_to_deg

            if self.coord_type == "longitude":
                value = wrap_angle_at(value, self.coord_wrap.to_value(u.deg))
            value = value * u.degree
            value = value.to_value(fl._unit)

        spacing = self._fl_spacing
        string = fl.formatter(values=[value] * fl._unit, spacing=spacing, format=format)

        return string[0]

    def set_separator(self, separator):
        """
        Set the separator to use for the angle major tick labels.

        Parameters
        ----------
        separator : str or tuple or None
            The separator between numbers in sexagesimal representation. Can be
            either a string or a tuple (or `None` for default).
        """
        if not (self._formatter_locator.__class__ == AngleFormatterLocator):
            raise TypeError("Separator can only be specified for angle coordinates")
        if isinstance(separator, (str, tuple)) or separator is None:
            self._formatter_locator.sep = separator
        else:
            raise TypeError("separator should be a string, a tuple, or None")

    def set_format_unit(self, unit, decimal=None, show_decimal_unit=True):
        """
        Set the unit for the major tick labels.

        Parameters
        ----------
        unit : class:`~astropy.units.Unit`
            The unit to which the tick labels should be converted to.
        decimal : bool, optional
            Whether to use decimal formatting. By default this is `False`
            for degrees or hours (which therefore use sexagesimal formatting)
            and `True` for all other units.
        show_decimal_unit : bool, optional
            Whether to include units when in decimal mode.
        """
        self._formatter_locator.format_unit = u.Unit(unit)
        self._formatter_locator.decimal = decimal
        self._formatter_locator.show_decimal_unit = show_decimal_unit

    def get_format_unit(self):
        """
        Get the unit for the major tick labels.
        """
        return self._formatter_locator.format_unit

    def set_ticks(
        self,
        values=None,
        spacing=None,
        number=None,
        size=None,
        width=None,
        color=None,
        alpha=None,
        direction=None,
        exclude_overlapping=None,
    ):
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
        color : str or tuple, optional
            A valid Matplotlib color for the ticks
        alpha : float, optional
            The alpha value (transparency) for the ticks.
        direction : {'in','out'}, optional
            Whether the ticks should point inwards or outwards.
        """
        if sum([values is None, spacing is None, number is None]) < 2:
            raise ValueError(
                "At most one of values, spacing, or number should be specified"
            )

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

        if direction is not None:
            if direction in ("in", "out"):
                self.ticks.set_tick_out(direction == "out")
            else:
                raise ValueError("direction should be 'in' or 'out'")

        if exclude_overlapping is not None:
            warnings.warn(
                "exclude_overlapping= should be passed to "
                "set_ticklabel instead of set_ticks",
                AstropyDeprecationWarning,
            )
            self.ticklabels.set_exclude_overlapping(exclude_overlapping)

    def set_ticks_position(self, position):
        """
        Set where ticks should appear.

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

    def set_ticklabel(
        self, color=None, size=None, pad=None, exclude_overlapping=None, **kwargs
    ):
        """
        Set the visual properties for the tick labels.

        Parameters
        ----------
        size : float, optional
            The size of the ticks labels in points
        color : str or tuple, optional
            A valid Matplotlib color for the tick labels
        pad : float, optional
            Distance in points between tick and label.
        exclude_overlapping : bool, optional
            Whether to exclude tick labels that overlap over each other.
        **kwargs
            Other keyword arguments are passed to :class:`matplotlib.text.Text`.
        """
        if size is not None:
            self.ticklabels.set_size(size)
        if color is not None:
            self.ticklabels.set_color(color)
        if pad is not None:
            self.ticklabels.set_pad(pad)
        if exclude_overlapping is not None:
            self.ticklabels.set_exclude_overlapping(exclude_overlapping)
        self.ticklabels.set(**kwargs)

    def set_ticklabel_position(self, position):
        """
        Set where tick labels should appear.

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
        **kwargs
            Keywords are passed to :class:`matplotlib.text.Text`. These
            can include keywords to set the ``color``, ``size``, ``weight``, and
            other text properties.
        """
        fontdict = kwargs.pop("fontdict", None)

        # NOTE: When using plt.xlabel/plt.ylabel, minpad can get set explicitly
        # to None so we need to make sure that in that case we change to a
        # default numerical value.
        if minpad is None:
            minpad = 1

        self._axislabel_set = True

        self.axislabels.set_text(text)
        self.axislabels.set_minpad(minpad)
        self.axislabels.set(**kwargs)

        if fontdict is not None:
            self.axislabels.update(fontdict)

    def get_axislabel(self):
        """
        Get the text for the axis label.

        Returns
        -------
        label : str
            The axis label
        """
        if self._auto_axislabel and not self._axislabel_set:
            return self._get_default_axislabel()
        else:
            return self.axislabels.get_text()

    def set_auto_axislabel(self, auto_label):
        """
        Render default axis labels if no explicit label is provided.

        Parameters
        ----------
        auto_label : `bool`
            `True` if default labels will be rendered.
        """
        self._auto_axislabel = bool(auto_label)

    def get_auto_axislabel(self):
        """
        Render default axis labels if no explicit label is provided.

        Returns
        -------
        auto_axislabel : `bool`
            `True` if default labels will be rendered.
        """
        return self._auto_axislabel

    def _get_default_axislabel(self):
        unit = self.get_format_unit() or self.coord_unit

        if not unit or unit is u.one or self.coord_type in ("longitude", "latitude"):
            return f"{self.default_label}"
        else:
            return f"{self.default_label} [{unit:latex}]"

    def set_axislabel_position(self, position):
        """
        Set where axis labels should appear.

        Parameters
        ----------
        position : str
            The axes on which the axis label for this coordinate should
            appear. Should be a string containing zero or more of ``'b'``,
            ``'t'``, ``'l'``, ``'r'``. For example, ``'lb'`` will lead the
            axis label to be shown on the left and bottom axis.
        """
        self.axislabels.set_visible_axes(position)

    def set_axislabel_visibility_rule(self, rule):
        """
        Set the rule used to determine when the axis label is drawn.

        Parameters
        ----------
        rule : str
            If the rule is 'always' axis labels will always be drawn on the
            axis. If the rule is 'ticks' the label will only be drawn if ticks
            were drawn on that axis. If the rule is 'labels' the axis label
            will only be drawn if tick labels were drawn on that axis.
        """
        self.axislabels.set_visibility_rule(rule)

    def get_axislabel_visibility_rule(self, rule):
        """
        Get the rule used to determine when the axis label is drawn.
        """
        return self.axislabels.get_visibility_rule()

    @property
    def locator(self):
        return self._formatter_locator.locator

    @property
    def formatter(self):
        return self._formatter_locator.formatter

    def _draw_grid(self, renderer):
        renderer.open_group("grid lines")

        self._update_ticks()

        if self.grid_lines_kwargs["visible"]:
            if isinstance(self.frame, RectangularFrame1D):
                self._update_grid_lines_1d()
            else:
                if self._grid_type == "lines":
                    self._update_grid_lines()
                else:
                    self._update_grid_contour()

            if self._grid_type == "lines":
                frame_patch = self.frame.patch
                for path in self.grid_lines:
                    p = PathPatch(path, **self.grid_lines_kwargs)
                    p.set_clip_path(frame_patch)
                    p.draw(renderer)

            elif self._grid is not None:
                if MATPLOTLIB_LT_3_8:
                    for line in self._grid.collections:
                        line.set(**self.grid_lines_kwargs)
                        line.draw(renderer)
                else:
                    self._grid.set(**self.grid_lines_kwargs)
                    self._grid.draw(renderer)

        renderer.close_group("grid lines")

    def _draw_ticks(self, renderer, existing_bboxes):
        """
        Draw all ticks and ticklabels.

        Parameters
        ----------
        existing_bboxes : list[Bbox]
            All bboxes for ticks that have already been drawn by other
            coordinates.
        """
        renderer.open_group("ticks")
        self.ticks.draw(renderer)

        self.ticklabels._tick_out_size = self.ticks.out_size

        self.ticklabels._set_existing_bboxes(existing_bboxes)
        self.ticklabels.draw(renderer)

        renderer.close_group("ticks")

    def _draw_axislabels(self, renderer, bboxes, ticklabels_bbox, visible_ticks):
        # Render the default axis label if no axis label is set.
        if self._auto_axislabel and not self._axislabel_set:
            self.set_axislabel(self._get_default_axislabel())

        renderer.open_group("axis labels")

        self.axislabels.draw(
            renderer,
            bboxes=bboxes,
            ticklabels_bbox=ticklabels_bbox,
            coord_ticklabels_bbox=ticklabels_bbox[self],
            ticks_locs=self.ticks.ticks_locs,
            visible_ticks=visible_ticks,
        )

        renderer.close_group("axis labels")

    def _update_ticks(self):
        if self.coord_index is None:
            return

        # TODO: this method should be optimized for speed

        # Here we determine the location and rotation of all the ticks. For
        # each axis, we can check the intersections for the specific
        # coordinate and once we have the tick positions, we can use the WCS
        # to determine the rotations.

        # Find the range of coordinates in all directions
        coord_range = self.parent_map.get_coord_range()

        # First find the ticks we want to show
        tick_world_coordinates, self._fl_spacing = self.locator(
            *coord_range[self.coord_index]
        )

        if self.ticks.get_display_minor_ticks():
            minor_ticks_w_coordinates = self._formatter_locator.minor_locator(
                self._fl_spacing,
                self.get_minor_frequency(),
                *coord_range[self.coord_index],
            )

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
        # https://matplotlib.org/stable/tutorials/advanced/transforms_tutorial.html#the-transformation-pipeline
        transData = self.parent_axes.transData
        invertedTransLimits = transData.inverted()

        for axis, spine in frame.items():
            if spine.data.size == 0:
                continue

            if not isinstance(self.frame, RectangularFrame1D):
                # Determine tick rotation in display coordinates and compare to
                # the normal angle in display coordinates.

                pixel0 = spine.data
                world0 = spine.world[:, self.coord_index]
                if np.isnan(world0).all():
                    continue
                axes0 = transData.transform(pixel0)

                # Advance 2 pixels in figure coordinates
                pixel1 = axes0.copy()
                pixel1[:, 0] += 2.0
                pixel1 = invertedTransLimits.transform(pixel1)
                with np.errstate(invalid="ignore"):
                    world1 = self.transform.transform(pixel1)[:, self.coord_index]

                # Advance 2 pixels in figure coordinates
                pixel2 = axes0.copy()
                pixel2[:, 1] += 2.0 if self.frame.origin == "lower" else -2.0
                pixel2 = invertedTransLimits.transform(pixel2)
                with np.errstate(invalid="ignore"):
                    world2 = self.transform.transform(pixel2)[:, self.coord_index]

                dx = world1 - world0
                dy = world2 - world0

                # Rotate by 90 degrees
                dx, dy = -dy, dx

                if self.coord_type == "longitude":
                    if self._coord_scale_to_deg is not None:
                        dx *= self._coord_scale_to_deg
                        dy *= self._coord_scale_to_deg

                    # Here we wrap at 180 not self.coord_wrap since we want to
                    # always ensure abs(dx) < 180 and abs(dy) < 180
                    dx = wrap_angle_at(dx, 180.0)
                    dy = wrap_angle_at(dy, 180.0)

                tick_angle = np.degrees(np.arctan2(dy, dx))

                normal_angle_full = np.hstack(
                    [spine.normal_angle, spine.normal_angle[-1]]
                )
                with np.errstate(invalid="ignore"):
                    reset = ((normal_angle_full - tick_angle) % 360 > 90.0) & (
                        (tick_angle - normal_angle_full) % 360 > 90.0
                    )
                tick_angle[reset] -= 180.0

            else:
                rotation = 90 if axis == "b" else -90
                tick_angle = np.zeros((conf.frame_boundary_samples,)) + rotation

            # We find for each interval the starting and ending coordinate,
            # ensuring that we take wrapping into account correctly for
            # longitudes.
            w1 = spine.world[:-1, self.coord_index]
            w2 = spine.world[1:, self.coord_index]

            if self.coord_type == "longitude":
                if self._coord_scale_to_deg is not None:
                    w1 = w1 * self._coord_scale_to_deg
                    w2 = w2 * self._coord_scale_to_deg

                w1 = wrap_angle_at(w1, self.coord_wrap.to_value(u.deg))
                w2 = wrap_angle_at(w2, self.coord_wrap.to_value(u.deg))
                with np.errstate(invalid="ignore"):
                    w1[w2 - w1 > 180.0] += 360
                    w2[w1 - w2 > 180.0] += 360

                if self._coord_scale_to_deg is not None:
                    w1 = w1 / self._coord_scale_to_deg
                    w2 = w2 / self._coord_scale_to_deg

            # For longitudes, we need to check ticks as well as ticks + 360,
            # since the above can produce pairs such as 359 to 361 or 0.5 to
            # 1.5, both of which would match a tick at 0.75. Otherwise we just
            # check the ticks determined above.
            self._compute_ticks(tick_world_coordinates, spine, axis, w1, w2, tick_angle)

            if self.ticks.get_display_minor_ticks():
                self._compute_ticks(
                    minor_ticks_w_coordinates,
                    spine,
                    axis,
                    w1,
                    w2,
                    tick_angle,
                    ticks="minor",
                )

        # format tick labels, add to scene
        text = self.formatter(
            self.lbl_world * tick_world_coordinates.unit, spacing=self._fl_spacing
        )
        for kwargs, txt in zip(self.lblinfo, text):
            self.ticklabels.add(text=txt, **kwargs)

    def _compute_ticks(
        self, tick_world_coordinates, spine, axis, w1, w2, tick_angle, ticks="major"
    ):
        if self.coord_type == "longitude":
            tick_world_coordinates_values = tick_world_coordinates.to_value(u.deg)
            tick_world_coordinates_values = np.hstack(
                [tick_world_coordinates_values, tick_world_coordinates_values + 360]
            )
            tick_world_coordinates_values *= u.deg.to(self.coord_unit)
        else:
            tick_world_coordinates_values = tick_world_coordinates.to_value(
                self.coord_unit
            )

        for t in tick_world_coordinates_values:
            # Find steps where a tick is present. We have to check
            # separately for the case where the tick falls exactly on the
            # frame points, otherwise we'll get two matches, one for w1 and
            # one for w2.
            with np.errstate(invalid="ignore"):
                intersections = np.hstack(
                    [
                        np.nonzero((t - w1) == 0)[0],
                        np.nonzero(((t - w1) * (t - w2)) < 0)[0],
                    ]
                )

            # But we also need to check for intersection with the last w2
            if t - w2[-1] == 0:
                intersections = np.append(intersections, len(w2) - 1)

            # Loop over ticks, and find exact pixel coordinates by linear
            # interpolation
            for imin in intersections:
                imax = imin + 1

                if np.allclose(w1[imin], w2[imin], rtol=1.0e-13, atol=1.0e-13):
                    continue  # tick is exactly aligned with frame
                else:
                    frac = (t - w1[imin]) / (w2[imin] - w1[imin])
                    x_data_i = spine.data[imin, 0] + frac * (
                        spine.data[imax, 0] - spine.data[imin, 0]
                    )
                    y_data_i = spine.data[imin, 1] + frac * (
                        spine.data[imax, 1] - spine.data[imin, 1]
                    )
                    delta_angle = tick_angle[imax] - tick_angle[imin]
                    if delta_angle > 180.0:
                        delta_angle -= 360.0
                    elif delta_angle < -180.0:
                        delta_angle += 360.0
                    angle_i = tick_angle[imin] + frac * delta_angle

                if self.coord_type == "longitude":
                    if self._coord_scale_to_deg is not None:
                        t *= self._coord_scale_to_deg

                    world = wrap_angle_at(t, self.coord_wrap.to_value(u.deg))

                    if self._coord_scale_to_deg is not None:
                        world /= self._coord_scale_to_deg

                else:
                    world = t

                if ticks == "major":
                    self.ticks.add(
                        axis=axis,
                        pixel=(x_data_i, y_data_i),
                        world=world,
                        angle=angle_i,
                        axis_displacement=imin + frac,
                    )

                    # store information to pass to ticklabels.add
                    # it's faster to format many ticklabels at once outside
                    # of the loop
                    self.lblinfo.append(
                        dict(
                            axis=axis,
                            data=(x_data_i, y_data_i),
                            world=world,
                            angle=spine.normal_angle[imin],
                            axis_displacement=imin + frac,
                        )
                    )
                    self.lbl_world.append(world)

                else:
                    self.ticks.add_minor(
                        minor_axis=axis,
                        minor_pixel=(x_data_i, y_data_i),
                        minor_world=world,
                        minor_angle=angle_i,
                        minor_axis_displacement=imin + frac,
                    )

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

    def _update_grid_lines_1d(self):
        if self.coord_index is None:
            return

        x_ticks_pos = [a[0] for a in self.ticks.pixel["b"]]

        ymin, ymax = self.parent_axes.get_ylim()

        self.grid_lines = []
        for x_coord in x_ticks_pos:
            pixel = [[x_coord, ymin], [x_coord, ymax]]
            self.grid_lines.append(Path(pixel))

    def _update_grid_lines(self):
        # For 3-d WCS with a correlated third axis, the *proper* way of
        # drawing a grid should be to find the world coordinates of all pixels
        # and drawing contours. What we are doing here assumes that we can
        # define the grid lines with just two of the coordinates (and
        # therefore assumes that the other coordinates are fixed and set to
        # the value in the slice). Here we basically assume that if the WCS
        # had a third axis, it has been abstracted away in the transformation.

        if self.coord_index is None:
            return

        coord_range = self.parent_map.get_coord_range()

        tick_world_coordinates, spacing = self.locator(*coord_range[self.coord_index])
        tick_world_coordinates_values = tick_world_coordinates.to_value(self.coord_unit)

        n_coord = len(tick_world_coordinates_values)
        if n_coord == 0:
            return

        from . import conf

        n_samples = conf.grid_samples

        xy_world = np.zeros((n_samples * n_coord, 2))

        self.grid_lines = []

        for iw, w in enumerate(tick_world_coordinates_values):
            subset = slice(iw * n_samples, (iw + 1) * n_samples)
            if self.coord_index == 0:
                xy_world[subset, 0] = np.repeat(w, n_samples)
                xy_world[subset, 1] = np.linspace(
                    coord_range[1][0], coord_range[1][1], n_samples
                )
            else:
                xy_world[subset, 0] = np.linspace(
                    coord_range[0][0], coord_range[0][1], n_samples
                )
                xy_world[subset, 1] = np.repeat(w, n_samples)

        # We now convert all the world coordinates to pixel coordinates in a
        # single go rather than doing this in the gridline to path conversion
        # to fully benefit from vectorized coordinate transformations.

        # Transform line to pixel coordinates
        pixel = self.transform.inverted().transform(xy_world)

        # Create round-tripped values for checking
        xy_world_round = self.transform.transform(pixel)

        for iw in range(n_coord):
            subset = slice(iw * n_samples, (iw + 1) * n_samples)
            self.grid_lines.append(
                self._get_gridline(
                    xy_world[subset], pixel[subset], xy_world_round[subset]
                )
            )

    def add_tickable_gridline(self, name, constant):
        """
        Define a gridline that can be used for ticks and labels.

        This gridline is not itself drawn, but instead can be specified in calls to
        methods such as
        :meth:`~astropy.visualization.wcsaxes.coordinate_helpers.CoordinateHelper.set_ticklabel_position`
        for drawing ticks and labels.  Since the gridline has a constant value in this
        coordinate, and thus would not have any ticks or labels for the same coordinate,
        the call to
        :meth:`~astropy.visualization.wcsaxes.coordinate_helpers.CoordinateHelper.set_ticklabel_position`
        would typically be made on the complementary coordinate.

        Parameters
        ----------
        name : str
            The name for the gridline, usually a single character, but can be longer
        constant : `~astropy.units.Quantity`
            The constant coordinate value of the gridline

        Notes
        -----
        A limitation is that the tickable part of the gridline must be contiguous.  If
        the gridline consists of more than one disconnected segment within the plot
        extent, only one of those segments will be made tickable.
        """
        if self.coord_index is None:
            return

        if name in self.frame:
            raise ValueError(f"The frame already has a spine with the name '{name}'")

        coord_range = self.parent_map.get_coord_range()
        constant = constant.to_value(self.coord_unit)

        from . import conf

        n_samples = conf.grid_samples

        # See comment in _update_grid_lines() about a WCS with more than 2 axes

        xy_world = np.zeros((n_samples, 2))
        xy_world[:, self.coord_index] = np.repeat(constant, n_samples)

        # If the complementary coordinate is longitude, we attempt to close the gridline
        # If such closure is a discontinuity, it will be filtered out later
        if self.parent_map[1 - self.coord_index].coord_type == "longitude":
            xy_world[:-1, 1 - self.coord_index] = np.linspace(
                coord_range[1 - self.coord_index][0],
                coord_range[1 - self.coord_index][1],
                n_samples - 1,
            )
            xy_world[-1, 1 - self.coord_index] = coord_range[1 - self.coord_index][0]
        else:
            xy_world[:, 1 - self.coord_index] = np.linspace(
                coord_range[1 - self.coord_index][0],
                coord_range[1 - self.coord_index][1],
                n_samples,
            )

        # Transform line to pixel coordinates
        pixel = self.transform.inverted().transform(xy_world)

        # Create round-tripped values for checking
        xy_world_round = self.transform.transform(pixel)

        # Get the path of the gridline, which masks hidden parts
        gridline = self._get_gridline(xy_world, pixel, xy_world_round)

        def data_for_spine(spine):
            vertices = gridline.vertices.copy()
            codes = gridline.codes.copy()

            # Retain the parts of the gridline within the rectangular plot bounds.
            # We ought to use the potentially non-rectangular plot frame, but
            # calculating that patch requires updating all spines first, which is a
            # catch-22.
            xmin, xmax = spine.parent_axes.get_xlim()
            ymin, ymax = spine.parent_axes.get_ylim()
            keep = (
                (vertices[:, 0] >= xmin)
                & (vertices[:, 0] <= xmax)
                & (vertices[:, 1] >= ymin)
                & (vertices[:, 1] <= ymax)
            )
            codes[~keep] = Path.MOVETO
            codes[1:][~keep[:-1]] = Path.MOVETO

            # We isolate the last segment (the last run of LINETOs), which must be preceded
            # by at least one MOVETO and may be succeeded by MOVETOs.
            # We have to account for longitude wrapping as well.

            # Bail out if there is no visible segment
            lineto = np.flatnonzero(codes == Path.LINETO)
            if np.size(lineto) == 0:
                return np.zeros((0, 2))

            # Find the start of the last segment (the last MOVETO before the LINETOs)
            last_segment = np.flatnonzero(codes[: lineto[-1]] == Path.MOVETO)[-1]

            # Double the gridline if it is closed (i.e., spans all longitudes)
            if vertices[0, 0] == vertices[-1, 0] and vertices[0, 1] == vertices[-1, 1]:
                codes = np.concatenate([codes, codes[1:]])
                vertices = np.vstack([vertices, vertices[1:, :]])

            # Stop the last segment before any trailing MOVETOs
            moveto = np.flatnonzero(codes[last_segment + 1 :] == Path.MOVETO)
            if np.size(moveto) > 0:
                return vertices[last_segment : last_segment + moveto[0] + 1, :]
            else:
                return vertices[last_segment:n_samples, :]

        self.frame[name] = self.frame.spine_class(
            self.frame.parent_axes, self.frame.transform, data_func=data_for_spine
        )

    def _get_gridline(self, xy_world, pixel, xy_world_round):
        if self.coord_type == "scalar":
            return get_gridline_path(xy_world, pixel)
        else:
            return get_lon_lat_path(xy_world, pixel, xy_world_round)

    def _clear_grid_contour(self):
        if hasattr(self, "_grid") and self._grid:
            if MATPLOTLIB_LT_3_8:
                for line in self._grid.collections:
                    line.remove()
            else:
                self._grid.remove()

    def _update_grid_contour(self):
        if self.coord_index is None:
            return

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        from . import conf

        res = conf.contour_grid_samples

        x, y = np.meshgrid(np.linspace(xmin, xmax, res), np.linspace(ymin, ymax, res))
        pixel = np.array([x.ravel(), y.ravel()]).T
        world = self.transform.transform(pixel)
        field = world[:, self.coord_index].reshape(res, res).T

        coord_range = self.parent_map.get_coord_range()

        tick_world_coordinates, spacing = self.locator(*coord_range[self.coord_index])

        # tick_world_coordinates is a Quantities array and we only needs its values
        tick_world_coordinates_values = tick_world_coordinates.value

        if self.coord_type == "longitude":
            # Find biggest gap in tick_world_coordinates and wrap in middle
            # For now just assume spacing is equal, so any mid-point will do
            mid = 0.5 * (
                tick_world_coordinates_values[0] + tick_world_coordinates_values[1]
            )
            field = wrap_angle_at(field, mid)
            tick_world_coordinates_values = wrap_angle_at(
                tick_world_coordinates_values, mid
            )

            # Replace wraps by NaN
            with np.errstate(invalid="ignore"):
                reset = (np.abs(np.diff(field[:, :-1], axis=0)) > 180) | (
                    np.abs(np.diff(field[:-1, :], axis=1)) > 180
                )
            field[:-1, :-1][reset] = np.nan
            field[1:, :-1][reset] = np.nan
            field[:-1, 1:][reset] = np.nan
            field[1:, 1:][reset] = np.nan

        if len(tick_world_coordinates_values) > 0:
            with np.errstate(invalid="ignore"):
                self._grid = self.parent_axes.contour(
                    x,
                    y,
                    field.transpose(),
                    levels=np.sort(tick_world_coordinates_values),
                )
        else:
            self._grid = None

    def tick_params(self, which="both", **kwargs):
        """
        Method to set the tick and tick label parameters in the same way as the
        :meth:`~matplotlib.axes.Axes.tick_params` method in Matplotlib.

        This is provided for convenience, but the recommended API is to use
        :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.set_ticks`,
        :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.set_ticklabel`,
        :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.set_ticks_position`,
        :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.set_ticklabel_position`,
        and :meth:`~astropy.visualization.wcsaxes.CoordinateHelper.grid`.

        Parameters
        ----------
        which : {'both', 'major', 'minor'}, optional
            Which ticks to apply the settings to. By default, setting are
            applied to both major and minor ticks. Note that if ``'minor'`` is
            specified, only the length of the ticks can be set currently.
        direction : {'in', 'out'}, optional
            Puts ticks inside the axes, or outside the axes.
        length : float, optional
            Tick length in points.
        width : float, optional
            Tick width in points.
        color : color, optional
            Tick color (accepts any valid Matplotlib color)
        pad : float, optional
            Distance in points between tick and label.
        labelsize : float or str, optional
            Tick label font size in points or as a string (e.g., 'large').
        labelcolor : color, optional
            Tick label color (accepts any valid Matplotlib color)
        colors : color, optional
            Changes the tick color and the label color to the same value
             (accepts any valid Matplotlib color).
        bottom, top, left, right : bool, optional
            Where to draw the ticks. Note that this will not work correctly if
            the frame is not rectangular.
        labelbottom, labeltop, labelleft, labelright : bool, optional
            Where to draw the tick labels. Note that this will not work
            correctly if the frame is not rectangular.
        grid_color : color, optional
            The color of the grid lines (accepts any valid Matplotlib color).
        grid_alpha : float, optional
            Transparency of grid lines: 0 (transparent) to 1 (opaque).
        grid_linewidth : float, optional
            Width of grid lines in points.
        grid_linestyle : str, optional
            The style of the grid lines (accepts any valid Matplotlib line
            style).
        """
        # First do some sanity checking on the keyword arguments

        # colors= is a fallback default for color and labelcolor
        if "colors" in kwargs:
            if "color" not in kwargs:
                kwargs["color"] = kwargs["colors"]
            if "labelcolor" not in kwargs:
                kwargs["labelcolor"] = kwargs["colors"]

        # The only property that can be set *specifically* for minor ticks is
        # the length. In future we could consider having a separate Ticks instance
        # for minor ticks so that e.g. the color can be set separately.
        if which == "minor":
            if len(set(kwargs) - {"length"}) > 0:
                raise ValueError(
                    "When setting which='minor', the only "
                    "property that can be set at the moment is "
                    "'length' (the minor tick length)"
                )
            else:
                if "length" in kwargs:
                    self.ticks.set_minor_ticksize(kwargs["length"])
            return

        # At this point, we can now ignore the 'which' argument.

        # Set the tick arguments
        self.set_ticks(
            size=kwargs.get("length"),
            width=kwargs.get("width"),
            color=kwargs.get("color"),
            direction=kwargs.get("direction"),
        )

        # Set the tick position
        position = None
        for arg in ("bottom", "left", "top", "right"):
            if arg in kwargs and position is None:
                position = ""
            if kwargs.get(arg):
                position += arg[0]
        if position is not None:
            self.set_ticks_position(position)

        # Set the tick label arguments.
        self.set_ticklabel(
            color=kwargs.get("labelcolor"),
            size=kwargs.get("labelsize"),
            pad=kwargs.get("pad"),
        )

        # Set the tick label position
        position = None
        for arg in ("bottom", "left", "top", "right"):
            if "label" + arg in kwargs and position is None:
                position = ""
            if kwargs.get("label" + arg):
                position += arg[0]
        if position is not None:
            self.set_ticklabel_position(position)

        # And the grid settings
        if "grid_color" in kwargs:
            self.grid_lines_kwargs["edgecolor"] = kwargs["grid_color"]
        if "grid_alpha" in kwargs:
            self.grid_lines_kwargs["alpha"] = kwargs["grid_alpha"]
        if "grid_linewidth" in kwargs:
            self.grid_lines_kwargs["linewidth"] = kwargs["grid_linewidth"]
        if "grid_linestyle" in kwargs:
            if kwargs["grid_linestyle"] in LINES_TO_PATCHES_LINESTYLE:
                self.grid_lines_kwargs["linestyle"] = LINES_TO_PATCHES_LINESTYLE[
                    kwargs["grid_linestyle"]
                ]
            else:
                self.grid_lines_kwargs["linestyle"] = kwargs["grid_linestyle"]

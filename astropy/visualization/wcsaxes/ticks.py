# Licensed under a 3-clause BSD style license - see LICENSE.rst

from collections import defaultdict

import numpy as np
from matplotlib import rcParams
from matplotlib.lines import Line2D, Path
from matplotlib.transforms import Affine2D


class Ticks(Line2D):
    """
    Ticks are derived from Line2D, and note that ticks themselves
    are markers. Thus, you should use set_mec, set_mew, etc.

    To change the tick size (length), you need to use
    set_ticksize. To change the direction of the ticks (ticks are
    in opposite direction of ticklabels by default), use
    set_tick_out(False).

    Note that Matplotlib's defaults dictionary :data:`~matplotlib.rcParams`
    contains default settings (color, size, width) of the form `xtick.*` and
    `ytick.*`. In a WCS projection, there may not be a clear relationship
    between axes of the projection and 'x' or 'y' axes. For this reason,
    we read defaults from `xtick.*`. The following settings affect the
    default appearance of ticks:

    * `xtick.direction`
    * `xtick.major.size`
    * `xtick.major.width`
    * `xtick.minor.size`
    * `xtick.color`

    Attributes
    ----------
    ticks_locs : dict
        This is set when the ticks are drawn, and is a mapping from axis to
        the locations of the ticks for that axis.
    """

    def __init__(self, frame=None, ticksize=None, **kwargs):
        self._frame = frame
        if ticksize is None:
            ticksize = rcParams["xtick.major.size"]
        self.set_ticksize(ticksize)
        self.set_minor_ticksize(rcParams["xtick.minor.size"])
        self.set_tick_out(rcParams["xtick.direction"] == "out")
        self.clear()
        line2d_kwargs = {
            "color": rcParams["xtick.color"],
            "linewidth": rcParams["xtick.major.width"],
        }
        line2d_kwargs.update(kwargs)
        Line2D.__init__(self, [0.0], [0.0], **line2d_kwargs)
        self.set_visible_axes("all")
        self._display_minor_ticks = False

    def display_minor_ticks(self, display_minor_ticks):
        self._display_minor_ticks = display_minor_ticks

    def get_display_minor_ticks(self):
        return self._display_minor_ticks

    def set_tick_out(self, tick_out):
        """
        set True if tick need to be rotated by 180 degree.
        """
        self._tick_out = tick_out

    def get_tick_out(self):
        """
        Return True if the tick will be rotated by 180 degree.
        """
        return self._tick_out

    def set_ticksize(self, ticksize):
        """
        set length of the ticks in points.
        """
        self._ticksize = ticksize

    def get_ticksize(self):
        """
        Return length of the ticks in points.
        """
        return self._ticksize

    def set_minor_ticksize(self, ticksize):
        """
        set length of the minor ticks in points.
        """
        self._minor_ticksize = ticksize

    def get_minor_ticksize(self):
        """
        Return length of the minor ticks in points.
        """
        return self._minor_ticksize

    @property
    def out_size(self):
        if self._tick_out:
            return self._ticksize
        else:
            return 0.0

    def set_visible_axes(self, visible_axes):
        self._visible_axes = self._frame._validate_positions(visible_axes)

    def get_visible_axes(self):
        if self._visible_axes == "all":
            return list(self._frame.keys())
        else:
            return [x for x in self._visible_axes if x in self._frame or x == "#"]

    def clear(self):
        self.world = defaultdict(list)
        self.pixel = defaultdict(list)
        self.angle = defaultdict(list)
        self.disp = defaultdict(list)
        self.minor_world = defaultdict(list)
        self.minor_pixel = defaultdict(list)
        self.minor_angle = defaultdict(list)
        self.minor_disp = defaultdict(list)

    def add(self, axis, world, pixel, angle, axis_displacement):
        self.world[axis].append(world)
        self.pixel[axis].append(pixel)
        self.angle[axis].append(angle)
        self.disp[axis].append(axis_displacement)

    def get_minor_world(self):
        return self.minor_world

    def add_minor(
        self, minor_axis, minor_world, minor_pixel, minor_angle, minor_axis_displacement
    ):
        self.minor_world[minor_axis].append(minor_world)
        self.minor_pixel[minor_axis].append(minor_pixel)
        self.minor_angle[minor_axis].append(minor_angle)
        self.minor_disp[minor_axis].append(minor_axis_displacement)

    def __len__(self):
        return len(self.world)

    _tickvert_path = Path([[0.0, 0.0], [1.0, 0.0]])

    def draw(self, renderer):
        """
        Draw the ticks.
        """
        self.ticks_locs = defaultdict(list)

        if not self.get_visible():
            return

        offset = renderer.points_to_pixels(self.get_ticksize())
        self._draw_ticks(renderer, self.pixel, self.angle, offset)
        if self._display_minor_ticks:
            offset = renderer.points_to_pixels(self.get_minor_ticksize())
            self._draw_ticks(renderer, self.minor_pixel, self.minor_angle, offset)

    def _draw_ticks(self, renderer, pixel_array, angle_array, offset):
        """
        Draw the minor ticks.
        """
        path_trans = self.get_transform()

        gc = renderer.new_gc()
        gc.set_foreground(self.get_color())
        gc.set_alpha(self.get_alpha())
        gc.set_linewidth(self.get_linewidth())

        marker_scale = Affine2D().scale(offset, offset)
        marker_rotation = Affine2D()
        marker_transform = marker_scale + marker_rotation

        initial_angle = 180.0 if self.get_tick_out() else 0.0

        for axis in self.get_visible_axes():
            if axis == "#":
                continue

            if axis not in pixel_array:
                continue

            for loc, angle in zip(pixel_array[axis], angle_array[axis]):
                # Set the rotation for this tick
                marker_rotation.rotate_deg(initial_angle + angle)

                # Draw the markers
                locs = path_trans.transform_non_affine(np.array([loc, loc]))
                renderer.draw_markers(
                    gc,
                    self._tickvert_path,
                    marker_transform,
                    Path(locs),
                    path_trans.get_affine(),
                )

                # Reset the tick rotation before moving to the next tick
                marker_rotation.clear()

                self.ticks_locs[axis].append(locs)

        gc.restore()

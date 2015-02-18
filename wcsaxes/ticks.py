# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from matplotlib.lines import Path, Line2D
from matplotlib.transforms import Affine2D
from matplotlib import rcParams


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

    * `xtick.major.size`
    * `xtick.major.width`
    * `xtick.color`
    """

    def __init__(self, ticksize=None, tick_out=False, **kwargs):
        if ticksize is None:
            ticksize = rcParams['xtick.major.size']
        self.set_ticksize(ticksize)
        self.set_tick_out(tick_out)
        # FIXME: tick_out is incompatible with Matplotlib tickdir option
        self.clear()
        line2d_kwargs = {
            'color': rcParams['xtick.color'],
            'linewidth': rcParams['xtick.major.width']
        }
        line2d_kwargs.update(kwargs)
        Line2D.__init__(self, [0.], [0.], **line2d_kwargs)
        self.set_visible_axes('all')
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

    def set_visible_axes(self, visible_axes):
        self._visible_axes = visible_axes

    def get_visible_axes(self):
        if self._visible_axes == 'all':
            return self.world.keys()
        else:
            return [x for x in self._visible_axes if x in self.world]

    def clear(self):
        self.world = {}
        self.pixel = {}
        self.angle = {}
        self.disp = {}
        self.minor_world = {}
        self.minor_pixel = {}
        self.minor_angle = {}
        self.minor_disp = {}

    def add(self, axis, world, pixel, angle, axis_displacement):
        if axis not in self.world:
            self.world[axis] = [world]
            self.pixel[axis] = [pixel]
            self.angle[axis] = [angle]
            self.disp[axis] = [axis_displacement]
        else:
            self.world[axis].append(world)
            self.pixel[axis].append(pixel)
            self.angle[axis].append(angle)
            self.disp[axis].append(axis_displacement)

    def get_minor_world(self):
        return self.minor_world

    def add_minor(self, minor_axis, minor_world, minor_pixel, minor_angle,
                  minor_axis_displacement):
        if minor_axis not in self.minor_world:
            self.minor_world[minor_axis] = [minor_world]
            self.minor_pixel[minor_axis] = [minor_pixel]
            self.minor_angle[minor_axis] = [minor_angle]
            self.minor_disp[minor_axis] = [minor_axis_displacement]
        else:
            self.minor_world[minor_axis].append(minor_world)
            self.minor_pixel[minor_axis].append(minor_pixel)
            self.minor_angle[minor_axis].append(minor_angle)
            self.minor_disp[minor_axis].append(minor_axis_displacement)

    def __len__(self):
        return len(self.world)

    _tickvert_path = Path([[0., 0.], [1., 0.]])

    def draw(self, renderer):
        """
        Draw the ticks.
        """

        if not self.get_visible():
            return

        offset = renderer.points_to_pixels(self.get_ticksize())
        self._draw_ticks(renderer, self.pixel, self.angle, offset)
        if self._display_minor_ticks:
            offset = offset * 0.5  # for minor ticksize
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

        initial_angle = 180. if self.get_tick_out() else 0.

        for axis in self.get_visible_axes():

            if not axis in pixel_array:
                continue

            for loc, angle in zip(pixel_array[axis], angle_array[axis]):

                # Set the rotation for this tick
                marker_rotation.rotate_deg(initial_angle + angle)

                # Draw the markers
                locs = path_trans.transform_non_affine(np.array([loc, loc]))
                renderer.draw_markers(gc, self._tickvert_path, marker_transform,
                                      Path(locs), path_trans.get_affine())

                # Reset the tick rotation before moving to the next tick
                marker_rotation.clear()

        gc.restore()

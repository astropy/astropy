import numpy as np

from matplotlib.lines import Path, Line2D
from matplotlib.transforms import Affine2D


class Ticks(Line2D):
    """
    Ticks are derived from Line2D, and note that ticks themselves
    are markers. Thus, you should use set_mec, set_mew, etc.

    To change the tick size (length), you need to use
    set_ticksize. To change the direction of the ticks (ticks are
    in opposite direction of ticklabels by default), use
    set_tick_out(False).
    """

    def __init__(self, ticksize=5., tick_out=False, **kwargs):
        self.set_ticksize(ticksize)
        self.set_tick_out(tick_out)
        self.clear()
        Line2D.__init__(self, [0.], [0.], **kwargs)
        self.set_color('black')
        self.set_visible_axes('all')

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

    def __len__(self):
        return len(self.world)

    _tickvert_path = Path([[0., 0.], [1., 0.]])

    def draw(self, renderer):
        """
        Draw the ticks.
        """

        if not self.get_visible():
            return

        path_trans = self.get_transform()

        gc = renderer.new_gc()
        gc.set_foreground(self.get_color())
        gc.set_alpha(self.get_alpha())

        offset = renderer.points_to_pixels(self.get_ticksize())
        marker_scale = Affine2D().scale(offset, offset)
        marker_rotation = Affine2D()
        marker_transform = marker_scale + marker_rotation

        initial_angle = 180. if self.get_tick_out() else 0.

        for axis in self.get_visible_axes():

            for loc, angle in zip(self.pixel[axis], self.angle[axis]):

                # Set the rotation for this tick
                marker_rotation.rotate_deg(initial_angle + angle)

                # Draw the markers
                locs = path_trans.transform_non_affine(np.array([loc, loc]))
                renderer.draw_markers(gc, self._tickvert_path, marker_transform,
                                      Path(locs), path_trans.get_affine())

                # Reset the tick rotation before moving to the next tick
                marker_rotation.clear()

        gc.restore()

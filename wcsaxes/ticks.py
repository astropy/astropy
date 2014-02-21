import numpy as np

from matplotlib.lines import Path, Line2D
from matplotlib.transforms import Affine2D

# Code partially adapted from mpl_toolkits.axisartist.
# TODO: include license


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
        self.set_color('white')
        Line2D.__init__(self, [0.], [0.], **kwargs)


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

    def set_color(self, color):
        """
        Set the color of the ticks
        """
        self._color = color

    def get_color(self):
        """
        Return the color of the ticks
        """
        return self._color

    def set_alpha(self, color):
        """
        Set the color of the ticks
        """
        self._alpha = alpha

    def get_alpha(self):
        """
        Return the color of the ticks
        """
        return self._alpha

    def clear(self):
        self.world = []
        self.pixel = []
        self.angle = []

    def add(self, world, pixel, angle):
        self.world.append(world)
        self.pixel.append(pixel)
        self.angle.append(angle)

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

        for loc, angle in zip(self.pixel, self.angle):

            # Set the rotation for this tick
            marker_rotation.rotate_deg(initial_angle + angle)

            # Draw the markers
            locs = path_trans.transform_non_affine(np.array([loc, loc]))
            renderer.draw_markers(gc, self._tickvert_path, marker_transform,
                                  Path(locs), path_trans.get_affine())

            # Reset the tick rotation before moving to the next tick
            marker_rotation.clear()

        gc.restore()

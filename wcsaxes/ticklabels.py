import numpy as np

from matplotlib.text import Text


class TickLabels(Text):

    def __init__(self, *args, **kwargs):
        self.clear()
        self._x_pixel_to_data = 1.
        self._y_pixel_to_data = 1.
        super(TickLabels, self).__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes('all')

    def clear(self):
        self.world = {}
        self.pixel = {}
        self.angle = {}
        self.text = {}

    def add(self, axis, world, pixel, angle, text):
        if axis not in self.world:
            self.world[axis] = [world]
            self.pixel[axis] = [pixel]
            self.angle[axis] = [angle]
            self.text[axis] = [text]
        else:
            self.world[axis].append(world)
            self.pixel[axis].append(pixel)
            self.angle[axis].append(angle)
            self.text[axis].append(text)

    def set_pixel_to_data_scaling(self, xscale, yscale):
        self._x_pixel_to_data = xscale
        self._y_pixel_to_data = yscale

    def set_visible_axes(self, visible_axes):
        self._visible_axes = visible_axes

    def get_visible_axes(self):
        if self._visible_axes == 'all':
            return self.world.keys()
        else:
            return self._visible_axes

    def draw(self, renderer):

        text_size = renderer.points_to_pixels(self.get_size())

        pad = text_size * 0.4

        for axis in self.get_visible_axes():

            for i in range(len(self.world[axis])):

                self.set_text(self.text[axis][i])

                # TODO: do something smarter for arbitrary directions
                if np.abs(self.angle[axis][i]) < 45.:
                    ha = 'right'
                    va = 'center'
                    dx = -pad
                    dy = 0.
                elif np.abs(self.angle[axis][i] - 90.) < 45:
                    ha = 'center'
                    va = 'top'
                    dx = 0
                    dy = -pad
                elif np.abs(self.angle[axis][i] - 180.) < 45:
                    ha = 'left'
                    va = 'center'
                    dx = pad
                    dy = 0
                else:
                    ha = 'center'
                    va = 'bottom'
                    dx = 0
                    dy = pad

                dx *= self._x_pixel_to_data
                dy *= self._y_pixel_to_data

                x, y = self.pixel[axis][i]

                self.set_position((x + dx, y + dy))
                self.set_ha(ha)
                self.set_va(va)

                super(TickLabels, self).draw(renderer)

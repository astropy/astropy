import numpy as np

from matplotlib.text import Text


class TickLabels(Text):

    def __init__(self, *args, **kwargs):
        self.clear()
        self._x_pixel_to_data = 1.
        self._y_pixel_to_data = 1.
        super(TickLabels, self).__init__(*args, **kwargs)
        self.set_clip_on(True)

    def clear(self):
        self.world = []
        self.pixel = []
        self.angle = []
        self.text = []

    def add(self, world, pixel, angle, text):
        self.world.append(world)
        self.pixel.append(pixel)
        self.angle.append(angle)
        self.text.append(text)

    def set_pixel_to_data_scaling(self, xscale, yscale):
        self._x_pixel_to_data = xscale
        self._y_pixel_to_data = yscale

    def draw(self, renderer):

        text_size = renderer.points_to_pixels(self.get_size())

        pad = text_size * 0.4

        for i in range(len(self.world)):

            self.set_text(self.text[i])

            # TODO: do something smarter for arbitrary directions
            if np.abs(self.angle[i]) < 45.:
                ha = 'right'
                va = 'center'
                dx = -pad
                dy = 0.
            elif np.abs(self.angle[i] - 90.) < 45:
                ha = 'center'
                va = 'top'
                dx = 0
                dy = -pad
            elif np.abs(self.angle[i] - 180.) < 45:
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

            x, y = self.pixel[i]

            self.set_position((x+dx, y+dy))
            self.set_ha(ha)
            self.set_va(va)

            super(TickLabels, self).draw(renderer)

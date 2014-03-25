import numpy as np

from matplotlib.text import Text


def sort_using(X, Y):
    return [x for (y,x) in sorted(zip(Y,X))]


class TickLabels(Text):

    def __init__(self, *args, **kwargs):
        self.clear()
        super(TickLabels, self).__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes('all')

    def clear(self):
        self.world = {}
        self.pixel = {}
        self.angle = {}
        self.text = {}
        self.disp = {}

    def add(self, axis, world, pixel, angle, text, axis_displacement):
        if axis not in self.world:
            self.world[axis] = [world]
            self.pixel[axis] = [pixel]
            self.angle[axis] = [angle]
            self.text[axis] = [text]
            self.disp[axis] = [axis_displacement]
        else:
            self.world[axis].append(world)
            self.pixel[axis].append(pixel)
            self.angle[axis].append(angle)
            self.text[axis].append(text)
            self.disp[axis].append(axis_displacement)

    def sort(self):
        """
        Sort by axis displacement, which allows us to figure out which parts
        of labels to not repeat.
        """
        for axis in self.world:
            self.world[axis] = sort_using(self.world[axis], self.disp[axis])
            self.pixel[axis] = sort_using(self.pixel[axis], self.disp[axis])
            self.angle[axis] = sort_using(self.angle[axis], self.disp[axis])
            self.text[axis] = sort_using(self.text[axis], self.disp[axis])
            self.disp[axis] = sort_using(self.disp[axis], self.disp[axis])

    def simplify_labels(self):
        """
        Figure out which parts of labels can be dropped to avoid repetition.
        """
        self.sort()
        for axis in self.world:
            t1 = self.text[axis][0]
            for i in range(1, len(self.world[axis])):
                t2 = self.text[axis][i]
                if len(t1) != len(t2):
                    t1 = self.text[axis][i]
                    continue
                start = 0
                for j in range(len(t1)):
                    if t1[j] != t2[j]:
                        break
                    if t1[j] not in '-0123456789.':
                        start = j + 1
                if start == 0:
                    t1 = self.text[axis][i]
                else:
                    self.text[axis][i] = self.text[axis][i][start:]

    def set_visible_axes(self, visible_axes):
        self._visible_axes = visible_axes

    def get_visible_axes(self):
        if self._visible_axes == 'all':
            return self.world.keys()
        else:
            return [x for x in self._visible_axes if x in self.world]

    def draw(self, renderer, bboxes):

        self.simplify_labels()

        text_size = renderer.points_to_pixels(self.get_size())

        for axis in self.get_visible_axes():

            for i in range(len(self.world[axis])):

                self.set_text(self.text[axis][i])

                # TODO: do something smarter for arbitrary directions
                if np.abs(self.angle[axis][i]) < 45.:
                    ha = 'right'
                    va = 'bottom'
                    dx = - text_size * 0.5
                    dy = - text_size * 0.5
                elif np.abs(self.angle[axis][i] - 90.) < 45:
                    ha = 'center'
                    va = 'bottom'
                    dx = 0
                    dy = - text_size * 1.5
                elif np.abs(self.angle[axis][i] - 180.) < 45:
                    ha = 'left'
                    va = 'bottom'
                    dx = text_size * 0.5
                    dy = - text_size * 0.5
                else:
                    ha = 'center'
                    va = 'bottom'
                    dx = 0
                    dy = text_size * 0.2

                x, y = self.pixel[axis][i]

                self.set_position((x + dx, y + dy))
                self.set_ha(ha)
                self.set_va(va)

                bb = super(TickLabels, self).get_window_extent(renderer)

                # TODO: the problem here is that we might get rid of a label
                # that has a key starting bit such as -0:30 where the -0
                # might be dropped from all other labels.

                if bb.count_overlaps(bboxes) == 0:
                    super(TickLabels, self).draw(renderer)
                    bboxes.append(bb)

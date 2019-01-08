# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np

from matplotlib import rcParams
from matplotlib.text import Text

from .frame import RectangularFrame


def sort_using(X, Y):
    return [x for (y, x) in sorted(zip(Y, X))]


class TickLabels(Text):

    def __init__(self, frame, *args, **kwargs):
        self.clear()
        self._frame = frame
        super().__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes('all')
        self.set_pad(rcParams['xtick.major.pad'])
        self._exclude_overlapping = False

        # Check rcParams

        if 'color' not in kwargs:
            self.set_color(rcParams['xtick.color'])

        if 'size' not in kwargs:
            self.set_size(rcParams['xtick.labelsize'])

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
                # In the following loop, we need to ignore the last character,
                # hence the len(t1) - 1. This is because if we have two strings
                # like 13d14m15s we want to make sure that we keep the last
                # part (15s) even if the two labels are identical.
                for j in range(len(t1) - 1):
                    if t1[j] != t2[j]:
                        break
                    if t1[j] not in '-0123456789.':
                        start = j + 1
                t1 = self.text[axis][i]
                if start != 0:
                    starts_dollar = self.text[axis][i].startswith('$')
                    self.text[axis][i] = self.text[axis][i][start:]
                    if starts_dollar:
                        self.text[axis][i] = '$' + self.text[axis][i]

    def set_pad(self, value):
        self._pad = value

    def get_pad(self):
        return self._pad

    def set_visible_axes(self, visible_axes):
        self._visible_axes = visible_axes

    def get_visible_axes(self):
        if self._visible_axes == 'all':
            return self.world.keys()
        else:
            return [x for x in self._visible_axes if x in self.world]

    def set_exclude_overlapping(self, exclude_overlapping):
        self._exclude_overlapping = exclude_overlapping

    def draw(self, renderer, bboxes, ticklabels_bbox, tick_out_size):

        if not self.get_visible():
            return

        self.simplify_labels()

        text_size = renderer.points_to_pixels(self.get_size())

        for axis in self.get_visible_axes():

            for i in range(len(self.world[axis])):

                # In the event that the label is empty (which is not expected
                # but could happen in unforeseen corner cases), we should just
                # skip to the next label.
                if self.text[axis][i] == '':
                    continue

                self.set_text(self.text[axis][i])

                x, y = self.pixel[axis][i]

                pad = renderer.points_to_pixels(self.get_pad() + tick_out_size)

                if isinstance(self._frame, RectangularFrame):

                    # This is just to preserve the current results, but can be
                    # removed next time the reference images are re-generated.

                    if np.abs(self.angle[axis][i]) < 45.:
                        ha = 'right'
                        va = 'bottom'
                        dx = -pad
                        dy = -text_size * 0.5
                    elif np.abs(self.angle[axis][i] - 90.) < 45:
                        ha = 'center'
                        va = 'bottom'
                        dx = 0
                        dy = -text_size - pad
                    elif np.abs(self.angle[axis][i] - 180.) < 45:
                        ha = 'left'
                        va = 'bottom'
                        dx = pad
                        dy = -text_size * 0.5
                    else:
                        ha = 'center'
                        va = 'bottom'
                        dx = 0
                        dy = pad

                    self.set_position((x + dx, y + dy))
                    self.set_ha(ha)
                    self.set_va(va)

                else:

                    # This is the more general code for arbitrarily oriented
                    # axes

                    # Set initial position and find bounding box
                    self.set_position((x, y))
                    bb = super().get_window_extent(renderer)

                    # Find width and height, as well as angle at which we
                    # transition which side of the label we use to anchor the
                    # label.
                    width = bb.width
                    height = bb.height

                    # Project axis angle onto bounding box
                    ax = np.cos(np.radians(self.angle[axis][i]))
                    ay = np.sin(np.radians(self.angle[axis][i]))

                    # Set anchor point for label
                    if np.abs(self.angle[axis][i]) < 45.:
                        dx = width
                        dy = ay * height
                    elif np.abs(self.angle[axis][i] - 90.) < 45:
                        dx = ax * width
                        dy = height
                    elif np.abs(self.angle[axis][i] - 180.) < 45:
                        dx = -width
                        dy = ay * height
                    else:
                        dx = ax * width
                        dy = -height

                    dx *= 0.5
                    dy *= 0.5

                    # Find normalized vector along axis normal, so as to be
                    # able to nudge the label away by a constant padding factor

                    dist = np.hypot(dx, dy)

                    ddx = dx / dist
                    ddy = dy / dist

                    dx += ddx * pad
                    dy += ddy * pad

                    self.set_position((x - dx, y - dy))
                    self.set_ha('center')
                    self.set_va('center')

                bb = super().get_window_extent(renderer)

                # TODO: the problem here is that we might get rid of a label
                # that has a key starting bit such as -0:30 where the -0
                # might be dropped from all other labels.

                if not self._exclude_overlapping or bb.count_overlaps(bboxes) == 0:
                    super().draw(renderer)
                    bboxes.append(bb)
                    ticklabels_bbox[axis].append(bb)

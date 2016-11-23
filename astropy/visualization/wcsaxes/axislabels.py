# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from matplotlib.text import Text
import matplotlib.transforms as mtransforms

from .frame import RectangularFrame


class AxisLabels(Text):

    def __init__(self, frame, minpad=1, *args, **kwargs):
        self._frame = frame
        super(AxisLabels, self).__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes('all')
        self.set_ha('center')
        self.set_va('center')
        self._minpad = minpad

    def get_minpad(self, axis):
        try:
            return self._minpad[axis]
        except TypeError:
            return self._minpad

    def set_visible_axes(self, visible_axes):
        self._visible_axes = visible_axes

    def get_visible_axes(self):
        if self._visible_axes == 'all':
            return self._frame.keys()
        else:
            return [x for x in self._visible_axes if x in self._frame]

    def set_minpad(self, minpad):
        self._minpad = minpad

    def draw(self, renderer, bboxes, ticklabels_bbox_list, visible_ticks):

        if not self.get_visible():
            return

        text_size = renderer.points_to_pixels(self.get_size())

        for axis in self.get_visible_axes():

            padding = text_size * self.get_minpad(axis)

            # Find position of the axis label. For now we pick the mid-point
            # along the path but in future we could allow this to be a
            # parameter.
            x_disp, y_disp = self._frame[axis].pixel[:, 0], self._frame[axis].pixel[:, 1]
            d = np.hstack([0., np.cumsum(np.sqrt(np.diff(x_disp) ** 2 + np.diff(y_disp) ** 2))])
            xcen = np.interp(d[-1] / 2., d, x_disp)
            ycen = np.interp(d[-1] / 2., d, y_disp)

            # Find segment along which the mid-point lies
            imin = np.searchsorted(d, d[-1] / 2.) - 1

            # Find normal of the axis label facing outwards on that segment
            normal_angle = self._frame[axis].normal_angle[imin] + 180.

            label_angle = (normal_angle - 90.) % 360.
            if label_angle < 225 and label_angle > 135:
                label_angle += 180
            self.set_rotation(label_angle)

            # Find label position by looking at the bounding box of ticks'
            # labels and the image. It sets the default padding at 1 times the
            # axis label font size which can also be changed by setting
            # the minpad parameter.

            if isinstance(self._frame, RectangularFrame):

                if len(ticklabels_bbox_list) > 0:
                    ticklabels_bbox = mtransforms.Bbox.union(ticklabels_bbox_list)
                else:
                    ticklabels_bbox = None

                if axis == 'l':
                    if axis in visible_ticks and ticklabels_bbox is not None:
                        left = ticklabels_bbox.xmin
                    else:
                        left = xcen
                    xpos = left - padding
                    self.set_position((xpos, ycen))

                elif axis == 'r':
                    if axis in visible_ticks and ticklabels_bbox is not None:
                        right = ticklabels_bbox.x1
                    else:
                        right = xcen
                    xpos = right + padding
                    self.set_position((xpos, ycen))

                elif axis == 'b':
                    if axis in visible_ticks and ticklabels_bbox is not None:
                        bottom = ticklabels_bbox.ymin
                    else:
                        bottom = ycen
                    ypos = bottom - padding
                    self.set_position((xcen, ypos))

                elif axis == 't':
                    if axis in visible_ticks and ticklabels_bbox is not None:
                        top = ticklabels_bbox.y1
                    else:
                        top = ycen
                    ypos = top + padding
                    self.set_position((xcen, ypos))

            else:  # arbitrary axis

                dx = np.cos(np.radians(normal_angle)) * (padding + text_size * 1.5)
                dy = np.sin(np.radians(normal_angle)) * (padding + text_size * 1.5)

                self.set_position((xcen + dx, ycen + dy))

            super(AxisLabels, self).draw(renderer)

            bb = super(AxisLabels, self).get_window_extent(renderer)
            bboxes.append(bb)

import numpy as np

from matplotlib.text import Text
import matplotlib.transforms as mtransforms


class AxisLabels(Text):

    def __init__(self, frame, minpad=1, *args, **kwargs):
        self._frame = frame
        super(AxisLabels, self).__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes('all')
        self.set_ha('center')
        self.set_va('center')
        self._minpad = minpad
        # self._pad = None

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

        text_size = renderer.points_to_pixels(self.get_size())
        self.pad = text_size * self._minpad

        for axis in self.get_visible_axes():

            # Find position of the axis label. For now we pick the mid-point
            # along the path but in future we could allow this to be a
            # parameter.
            x_disp, y_disp = self._frame[axis].pixel[:,0], self._frame[axis].pixel[:,1]
            d = np.hstack([0., np.cumsum(np.sqrt(np.diff(x_disp)**2 + np.diff(y_disp)**2))])
            xcen = np.interp(d[-1]/2., d, x_disp)
            ycen = np.interp(d[-1]/2., d, y_disp)

            # Find segment along which the mid-point lies
            imin = np.searchsorted(d, d[-1] / 2.) - 1

            # Find normal of the axis label facing outwards on that segment
            normal_angle = self._frame[axis].normal_angle[imin] + 180.

            label_angle = (normal_angle - 90.) % 360.
            if label_angle < 225 and label_angle > 135:
                label_angle += 180
            self.set_rotation(label_angle)

            # Find label position by looking at the bouding box of ticks'
            # labels and the image. It places the label 1 font size away
            # from the axis, which can also be changed by setting
            # the minpad parameter.

            ticklabels_bbox = mtransforms.Bbox.union(ticklabels_bbox_list)

            if axis == 'l':
                if axis in visible_ticks:
                    left = ticklabels_bbox.xmin
                else:
                    left = xcen
                xpos = left - (self.pad)
                self.set_position((xpos, ycen))

            elif axis == 'r':
                if axis in visible_ticks:
                    right = ticklabels_bbox.x1
                else:
                    right = xcen
                xpos = right + (self.pad)
                self.set_position((xpos, ycen))

            elif axis == 'b':
                if axis in visible_ticks:
                    bottom = ticklabels_bbox.ymin
                else:
                    bottom = ycen
                ypos = bottom - (self.pad)
                self.set_position((xcen, ypos))

            elif axis == 't':
                if axis in visible_ticks:
                    top = ticklabels_bbox.y1
                else:
                    top = ycen
                ypos = top + (self.pad)
                self.set_position((xcen, ypos))

            super(AxisLabels, self).draw(renderer)

            bb = super(AxisLabels, self).get_window_extent(renderer)
            bboxes.append(bb)

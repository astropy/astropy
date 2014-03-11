import numpy as np

from matplotlib.text import Text


class AxisLabels(Text):

    def __init__(self, frame, *args, **kwargs):
        self._frame = frame
        super(AxisLabels, self).__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes('all')
        self.set_ha('center')
        self.set_va('center')

    def set_visible_axes(self, visible_axes):
        self._visible_axes = visible_axes

    def get_visible_axes(self):
        if self._visible_axes == 'all':
            return self._frame.keys()
        else:
            return [x for x in self._visible_axes if x in self._frame]

    def draw(self, renderer, bboxes):

        text_size = renderer.points_to_pixels(self.get_size())

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

            # Find label position, by trying to move it successively further
            # away from the axis and then checking for intersection with tick
            # labels. Obviously this could be optimized if needed.
            dx = np.cos(np.radians(normal_angle)) * text_size
            dy = np.sin(np.radians(normal_angle)) * text_size

            for pad in np.linspace(0.8, 8.8, 21):

                xlabel = xcen + dx * pad
                ylabel = ycen + dy * pad

                self.set_position((xlabel, ylabel))

                bb = super(AxisLabels, self).get_window_extent(renderer)

                if bb.count_overlaps(bboxes) == 0:
                    break

            self.set_position((xlabel + dx * pad * 0.3, ylabel + dy * pad * 0.3))
            super(AxisLabels, self).draw(renderer)

            bb = super(AxisLabels, self).get_window_extent(renderer)
            bboxes.append(bb)

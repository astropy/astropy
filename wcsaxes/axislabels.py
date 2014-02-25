import numpy as np

from matplotlib.text import Text


class AxisLabels(Text):

    def __init__(self, frame, *args, **kwargs):
        self._frame = frame
        self._x_pixel_to_data = 1.
        self._y_pixel_to_data = 1.
        super(AxisLabels, self).__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes('all')
        self.set_ha('center')
        self.set_va('center')

    def set_pixel_to_data_scaling(self, xscale, yscale):
        self._x_pixel_to_data = xscale
        self._y_pixel_to_data = yscale

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
            x, y = self._frame[axis]
            d = np.hstack([0., np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))])
            xcen = np.interp(d[-1]/2., d, x)
            ycen = np.interp(d[-1]/2., d, y)

            # Find segment along which the mid-point lies
            imin = np.searchsorted(d, d[-1] / 2.) - 1

            # Find normal of the axis label facing outwards on that segment
            normal_angle = np.arctan2(x[imin+1] - x[imin], -(y[imin+1] - y[imin])) + np.pi

            label_angle = (np.degrees(normal_angle) - 90.) % 360.
            if label_angle < 225 and label_angle > 135:
                label_angle += 180
            self.set_rotation(label_angle)

            # Find label position, by trying to move it successively further
            # away from the axis and then checking for intersection with tick
            # labels. Obviously this could be optimized if needed.
            for pad in range(1, 20):

                xlabel = xcen + np.cos(normal_angle) * pad * text_size * self._x_pixel_to_data
                ylabel = ycen + np.sin(normal_angle) * pad * text_size * self._y_pixel_to_data

                self.set_position((xlabel, ylabel))

                bb = super(AxisLabels, self).get_window_extent(renderer)

                if bb.count_overlaps(bboxes) == 0:
                    break

            super(AxisLabels, self).draw(renderer)
            bboxes.append(bb)

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import abc
from collections import OrderedDict
import numpy as np

from astropy.extern import six

from matplotlib.lines import Line2D, Path
from matplotlib.patches import PathPatch

__all__ = ['Spine', 'BaseFrame', 'RectangularFrame', 'EllipticalFrame']


class Spine(object):

    def __init__(self, parent_axes, transform):

        self.parent_axes = parent_axes
        self.transform = transform

        self.data = None
        self.pixel = None
        self.world = None

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        if value is None:
            self._data = None
            self._pixel = None
            self._world = None
        else:
            self._data = value
            self._pixel = self.parent_axes.transData.transform(self._data)
            self._world = self.transform.transform(self._data)
            self._update_normal()

    @property
    def pixel(self):
        return self._pixel

    @pixel.setter
    def pixel(self, value):
        if value is None:
            self._data = None
            self._pixel = None
            self._world = None
        else:
            self._data = self.parent_axes.transData.inverted().transform(self._data)
            self._pixel = value
            self._world = self.transform.transform(self._data)
            self._update_normal()

    @property
    def world(self):
        return self._world

    @world.setter
    def world(self, value):
        if value is None:
            self._data = None
            self._pixel = None
            self._world = None
        else:
            self._data = self.transform.transform(value)
            self._pixel = self.parent_axes.transData.transform(self._data)
            self._world = value
            self._update_normal()

    def _update_normal(self):
        # Find angle normal to border and inwards, in display coordinate
        dx = self.pixel[1:, 0] - self.pixel[:-1, 0]
        dy = self.pixel[1:, 1] - self.pixel[:-1, 1]
        self.normal_angle = np.degrees(np.arctan2(dx, -dy))


@six.add_metaclass(abc.ABCMeta)
class BaseFrame(OrderedDict):

    def __init__(self, parent_axes, transform, path=None):

        super(BaseFrame, self).__init__()

        self.parent_axes = parent_axes
        self._transform = transform
        self._linewidth = None
        self._color = 'black'
        self._path = path

        for axis in self.spine_names:
            self[axis] = Spine(parent_axes, transform)

    @property
    def origin(self):
        ymin, ymax = self.parent_axes.get_ylim()
        return 'lower' if ymin < ymax else 'upper'

    @property
    def transform(self):
        return self._transform

    @transform.setter
    def transform(self, value):
        self._transform = value
        for axis in self:
            self[axis].transform = value

    def _update_patch_path(self):

        self.update_spines()
        x, y = [], []
        for axis in self:
            x.append(self[axis].data[:, 0])
            y.append(self[axis].data[:, 1])
        vertices = np.vstack([np.hstack(x), np.hstack(y)]).transpose()

        if self._path is None:
            self._path = Path(vertices)
        else:
            self._path.vertices = vertices

    @property
    def patch(self):
        self._update_patch_path()
        return PathPatch(self._path, transform=self.parent_axes.transData,
                         facecolor='white', edgecolor='white')

    def draw(self, renderer):
        for axis in self:
            x, y = self[axis].pixel[:, 0], self[axis].pixel[:, 1]
            line = Line2D(x, y, linewidth=self._linewidth, color=self._color, zorder=1000)
            line.draw(renderer)

    def sample(self, n_samples):

        self.update_spines()

        spines = OrderedDict()

        for axis in self:

            data = self[axis].data
            p = np.linspace(0., 1., data.shape[0])
            p_new = np.linspace(0., 1., n_samples)
            spines[axis] = Spine(self.parent_axes, self.transform)
            spines[axis].data = np.array([np.interp(p_new, p, data[:, 0]),
                                          np.interp(p_new, p, data[:, 1])]).transpose()

        return spines

    def set_color(self, color):
        """
        Sets the color of the frame.

        Parameters
        ----------
        color : string
            The color of the frame.
        """
        self._color = color

    def get_color(self):
        return self._color

    def set_linewidth(self, linewidth):
        """
        Sets the linewidth of the frame.

        Parameters
        ----------
        linewidth : float
            The linewidth of the frame in points.
        """
        self._linewidth = linewidth

    def get_linewidth(self):
        return self._linewidth

    @abc.abstractmethod
    def update_spines(self):
        raise NotImplementedError("")


class RectangularFrame(BaseFrame):

    spine_names = 'brtl'

    def update_spines(self):

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        self['b'].data = np.array(([xmin, ymin], [xmax, ymin]))
        self['r'].data = np.array(([xmax, ymin], [xmax, ymax]))
        self['t'].data = np.array(([xmax, ymax], [xmin, ymax]))
        self['l'].data = np.array(([xmin, ymax], [xmin, ymin]))


class EllipticalFrame(BaseFrame):

    spine_names = 'c'

    def update_spines(self):

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        xmid = 0.5 * (xmax + xmin)
        ymid = 0.5 * (ymax + ymin)

        dx = xmid - xmin
        dy = ymid - ymin

        theta = np.linspace(0., 2 * np.pi, 1000)
        self['c'].data = np.array([xmid + dx * np.cos(theta),
                                   ymid + dy * np.sin(theta)]).transpose()

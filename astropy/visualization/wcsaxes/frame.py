# Licensed under a 3-clause BSD style license - see LICENSE.rst


import abc
from collections import OrderedDict
import warnings

import numpy as np


from matplotlib import rcParams
from matplotlib.lines import Line2D, Path
from matplotlib.patches import PathPatch

from astropy.utils.exceptions import AstropyDeprecationWarning

__all__ = ['RectangularFrame1D', 'Spine', 'BaseFrame', 'RectangularFrame', 'EllipticalFrame']


class Spine:
    """
    A single side of an axes.

    This does not need to be a straight line, but represents a 'side' when
    determining which part of the frame to put labels and ticks on.
    """

    def __init__(self, parent_axes, transform):

        self.parent_axes = parent_axes
        self.transform = transform

        self._data = None
        self._world = None

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        if value is None:
            self._data = None
            self._world = None
        else:
            self._data = value
            with np.errstate(invalid='ignore'):
                self._world = self.transform.transform(self._data)
            self._update_normal()

    def _get_pixel(self):
        return self.parent_axes.transData.transform(self._data)

    @property
    def pixel(self):
        warnings.warn('Pixel coordinates cannot be accurately calculated unless '
                      'Matplotlib is currently drawing a figure, so the .pixel '
                      'attribute is deprecated and will be removed in a future '
                      'astropy release.', AstropyDeprecationWarning)
        return self._get_pixel()

    @pixel.setter
    def pixel(self, value):
        warnings.warn('Manually setting pixel values of a Spine can lead to incorrect results '
                      'as these can only be accuractley calculated when Matplotlib is drawing '
                      'a figure. As such the .pixel setter now does nothing, is deprecated, '
                      'and will be removed in a future astropy release.', AstropyDeprecationWarning)

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
        pixel = self._get_pixel()
        # Find angle normal to border and inwards, in display coordinate
        dx = pixel[1:, 0] - pixel[:-1, 0]
        dy = pixel[1:, 1] - pixel[:-1, 1]
        self.normal_angle = np.degrees(np.arctan2(dx, -dy))

    def _halfway_x_y_angle(self):
        """
        Return the x, y, normal_angle values halfway along the spine
        """
        pixel = self._get_pixel()
        x_disp, y_disp = pixel[:, 0], pixel[:, 1]
        # Get distance along the path
        d = np.hstack([0., np.cumsum(np.sqrt(np.diff(x_disp) ** 2 + np.diff(y_disp) ** 2))])
        xcen = np.interp(d[-1] / 2., d, x_disp)
        ycen = np.interp(d[-1] / 2., d, y_disp)

        # Find segment along which the mid-point lies
        imin = np.searchsorted(d, d[-1] / 2.) - 1

        # Find normal of the axis label facing outwards on that segment
        normal_angle = self.normal_angle[imin] + 180.
        return xcen, ycen, normal_angle


class SpineXAligned(Spine):
    """
    A single side of an axes, aligned with the X data axis.

    This does not need to be a straight line, but represents a 'side' when
    determining which part of the frame to put labels and ticks on.
    """

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        if value is None:
            self._data = None
            self._world = None
        else:
            self._data = value
            with np.errstate(invalid='ignore'):
                self._world = self.transform.transform(self._data[:, 0:1])
            self._update_normal()


class BaseFrame(OrderedDict, metaclass=abc.ABCMeta):
    """
    Base class for frames, which are collections of
    :class:`~astropy.visualization.wcsaxes.frame.Spine` instances.
    """

    spine_class = Spine

    def __init__(self, parent_axes, transform, path=None):

        super().__init__()

        self.parent_axes = parent_axes
        self._transform = transform
        self._linewidth = rcParams['axes.linewidth']
        self._color = rcParams['axes.edgecolor']
        self._path = path

        for axis in self.spine_names:
            self[axis] = self.spine_class(parent_axes, transform)

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
                         facecolor=rcParams['axes.facecolor'], edgecolor='white')

    def draw(self, renderer):
        for axis in self:
            pixel = self[axis]._get_pixel()
            x, y = pixel[:, 0], pixel[:, 1]
            line = Line2D(x, y, linewidth=self._linewidth, color=self._color, zorder=1000)
            line.draw(renderer)

    def sample(self, n_samples):

        self.update_spines()

        spines = OrderedDict()

        for axis in self:

            data = self[axis].data
            p = np.linspace(0., 1., data.shape[0])
            p_new = np.linspace(0., 1., n_samples)
            spines[axis] = self.spine_class(self.parent_axes, self.transform)
            spines[axis].data = np.array([np.interp(p_new, p, d) for d in data.T]).transpose()

        return spines

    def set_color(self, color):
        """
        Sets the color of the frame.

        Parameters
        ----------
        color : str
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


class RectangularFrame1D(BaseFrame):
    """
    A classic rectangular frame.
    """

    spine_names = 'bt'
    spine_class = SpineXAligned

    def update_spines(self):

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        self['b'].data = np.array(([xmin, ymin], [xmax, ymin]))
        self['t'].data = np.array(([xmax, ymax], [xmin, ymax]))

    def _update_patch_path(self):

        self.update_spines()

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        x = [xmin, xmax, xmax, xmin, xmin]
        y = [ymin, ymin, ymax, ymax, ymin]

        vertices = np.vstack([np.hstack(x), np.hstack(y)]).transpose()

        if self._path is None:
            self._path = Path(vertices)
        else:
            self._path.vertices = vertices

    def draw(self, renderer):
        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        x = [xmin, xmax, xmax, xmin, xmin]
        y = [ymin, ymin, ymax, ymax, ymin]

        line = Line2D(x, y, linewidth=self._linewidth, color=self._color, zorder=1000,
                      transform=self.parent_axes.transData)
        line.draw(renderer)


class RectangularFrame(BaseFrame):
    """
    A classic rectangular frame.
    """

    spine_names = 'brtl'

    def update_spines(self):

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        self['b'].data = np.array(([xmin, ymin], [xmax, ymin]))
        self['r'].data = np.array(([xmax, ymin], [xmax, ymax]))
        self['t'].data = np.array(([xmax, ymax], [xmin, ymax]))
        self['l'].data = np.array(([xmin, ymax], [xmin, ymin]))


class EllipticalFrame(BaseFrame):
    """
    An elliptical frame.
    """

    spine_names = 'chv'

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
        self['h'].data = np.array([np.linspace(xmin, xmax, 1000),
                                   np.repeat(ymid, 1000)]).transpose()
        self['v'].data = np.array([np.repeat(xmid, 1000),
                                   np.linspace(ymin, ymax, 1000)]).transpose()

    def _update_patch_path(self):
        """Override path patch to include only the outer ellipse,
        not the major and minor axes in the middle."""

        self.update_spines()
        vertices = self['c'].data

        if self._path is None:
            self._path = Path(vertices)
        else:
            self._path.vertices = vertices

    def draw(self, renderer):
        """Override to draw only the outer ellipse,
        not the major and minor axes in the middle.

        FIXME: we may want to add a general method to give the user control
        over which spines are drawn."""
        axis = 'c'
        pixel = self[axis]._get_pixel()
        x, y = pixel[:, 0], pixel[:, 1]
        line = Line2D(x, y, linewidth=self._linewidth, color=self._color, zorder=1000)
        line.draw(renderer)

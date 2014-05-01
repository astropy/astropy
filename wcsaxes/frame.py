import numpy as np
from astropy.utils import OrderedDict
from matplotlib.lines import Line2D, Path

# TODO: once we want to start writing more complex frames, use an abstract base
# class.

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
        dx = self.pixel[1:,0] - self.pixel[:-1,0]
        dy = self.pixel[1:,1] - self.pixel[:-1,1]
        self.normal_angle = np.degrees(np.arctan2(dx, -dy))


class RectangularFrame(OrderedDict):

    def __init__(self, parent_axes, transform):

        super(RectangularFrame, self).__init__()

        self.parent_axes = parent_axes
        self._transform = transform

        for axis in 'brtl':
            self[axis] = Spine(parent_axes, transform)

        self._update_cache = None
        self._sample_cache = {}

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

    def update(self):

        xmin, xmax = self.parent_axes.get_xlim()
        ymin, ymax = self.parent_axes.get_ylim()

        if self._update_cache == [xmin, xmax, ymin, ymax]:
            return
        self._update_cache = [xmin, xmax, ymin, ymax]
        self._sample_cache = {}

        self['b'].data = np.array(([xmin, xmax], [ymin, ymin])).transpose()
        self['r'].data = np.array(([xmax, xmax], [ymin, ymax])).transpose()
        self['t'].data = np.array(([xmax, xmin], [ymax, ymax])).transpose()
        self['l'].data = np.array(([xmin, xmin], [ymax, ymin])).transpose()

    def sample(self, n_samples):

        self.update()
        if n_samples in self._sample_cache:
            return self._sample_cache[n_samples]

        spines = OrderedDict()

        for axis in self:

            data = self[axis].data
            p = np.linspace(0., 1., data.shape[0])
            p_new = np.linspace(0., 1., n_samples)
            spines[axis] = Spine(self.parent_axes, self.transform)
            spines[axis].data = np.array([np.interp(p_new, p, data[:,0]),
                                          np.interp(p_new, p, data[:,1])]).transpose()

        return spines

    @property
    def path(self):
        x, y = [], []
        for axis in self:
            x.append(self[axis].pixel[:,0])
            y.append(self[axis].pixel[:,1])
        return Path(np.vstack([np.hstack(x), np.hstack(y)]).transpose())

    def draw(self, renderer):
        for axis in self:
            x, y = self[axis].pixel[:,0], self[axis].pixel[:,1]
            line = Line2D(x, y, color='black', zorder=1000)
            line.draw(renderer)

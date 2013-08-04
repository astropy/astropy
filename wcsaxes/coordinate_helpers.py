import abc

import numpy as np

from matplotlib.ticker import Formatter
from mpl_toolkits.axisartist import angle_helper

from .formatters import AngleFormatter


class BaseCoordinateHelper(object):

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def set_ticks_position(self):
        """
        Set the axes on which the ticks for this coordinate should
        appear. Should be a string containing zero or more of ``b``, ``t``,
        ``l``, ``r``.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def set_ticklabel_position(self, position):
        """
        Set the axes on which the ticklabels for this coordinate should
        appear. Should be a string containing zero or more of ``b``, ``t``,
        ``l``, ``r``.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def set_major_formatter(self, formatter):
        """
        Set the major formatter for the ticks - should be either a
        ``Formatter`` instance for world coordinates, or a string such as
        ``dd:mm:ss.s``.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def set_ticks(self, spacing=None, number=None):
        """
        Set the spacing/value of the ticks. This can take:

        * A list of tick position
        * A spacing with the ``spacing=`` option (can take quantities)
        * An approximate number of tick marks with the ``number=`` option
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def grid(self):
        """
        Draw grid lines for just this coordinate. Should return a
        ``LineCollection`` instance.
        """
        raise NotImplementedError()


class FixedAngleLocator(object):
    """
    Differs from FixedLocator because it is compatible with the grid helper
    from mpl_toolkits which requires 2 positional arguments.
    """
    def __init__(self, values):
        self.values = np.array(values)
    def __call__(self, lon_min, lon_max):
        return self.values, len(self.values), 1.0


class SkyCoordinateHelper(BaseCoordinateHelper):

    def __init__(self, index=None, parent=None):
        self._index = index
        self._parent = parent
        self.locator = angle_helper.LocatorDMS(4)
        self.formatter = AngleFormatter(precision=1)

    def set_ticks_position(self):
        raise NotImplementedError()

    def set_ticklabel_position(self):
        raise NotImplementedError()

    def set_major_formatter(self, formatter):
        if isinstance(formatter, Formatter):
            raise NotImplementedError()  # figure out how to swap out formatter
        elif isinstance(formatter, basestring):
            self.formatter.format = formatter
            # need to also update locator in this case
        else:
            raise TypeError("formatter should be a string for Formatter instance")

    def set_ticks(self, values=None, spacing=None, number=None):
        if values is not None:
            self.locator = FixedAngleLocator(values)
        elif number is not None:
            self.locator = angle_helper.LocatorDMS(number)
        else:
            raise NotImplementedError("spacing")

    def grid(self):
        raise NotImplementedError()

    @property
    def locator(self):
        return getattr(self._parent.grid_helper.grid_finder, 'grid_locator' + str(self._index))

    @locator.setter
    def locator(self, value):
        self._parent.grid_helper.update_grid_finder(**{'grid_locator' + str(self._index): value})

    @property
    def formatter(self):
        return getattr(self._parent.grid_helper.grid_finder, 'tick_formatter' + str(self._index))

    @formatter.setter
    def formatter(self, value):
        self._parent.grid_helper.update_grid_finder(**{'tick_formatter' + str(self._index): value})


class ScalarCoordinateHelper(BaseCoordinateHelper):
    pass

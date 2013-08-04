import abc

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


class SkyCoordinateHelper(BaseCoordinateHelper):

    def __init__(self):
        self.grid_locator = angle_helper.LocatorDMS(4)
        self.tick_formatter = AngleFormatter(precision=1)

    def set_ticks_position(self):
        raise NotImplementedError()

    def set_ticklabel_position(self):
        raise NotImplementedError()

    def set_major_formatter(self, formatter):
        if isinstance(formatter, Formatter):
            raise NotImplementedError()  # figure out how to swap out formatter
        elif isinstance(formatter, basestring):
            self.tick_formatter.format = formatter
        else:
            raise TypeError("formatter should be a string for Formatter instance")

    def set_ticks(self, spacing=None, number=None):
        raise NotImplementedError()

    def grid(self):
        raise NotImplementedError()


class ScalarCoordinateHelper(BaseCoordinateHelper):
    pass

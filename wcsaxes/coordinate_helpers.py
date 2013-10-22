import abc

from matplotlib.ticker import Formatter

from .formatter_locator import AngleFormatterLocator
from . import six


class BaseCoordinateHelper(object):

    __metaclass__ = abc.ABCMeta

    def set_ticks_position(self):
        """
        Set the axes on which the ticks for this coordinate should
        appear. Should be a string containing zero or more of ``b``, ``t``,
        ``l``, ``r``.
        """
        raise NotImplementedError()

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

    def grid(self):
        """
        Draw grid lines for just this coordinate. Should return a
        ``LineCollection`` instance.
        """
        raise NotImplementedError()

    @property
    def locator(self):
        return self._fl_helper.locator

    @property
    def formatter(self):
        return self._fl_helper.formatter


class SkyCoordinateHelper(BaseCoordinateHelper):

    def __init__(self):
        self._fl_helper = AngleFormatterLocator()
        self._grid_helper = None

    def set_major_formatter(self, formatter):
        if isinstance(formatter, Formatter):
            raise NotImplementedError()  # figure out how to swap out formatter
        elif isinstance(formatter, six.string_types):
            self._fl_helper.format = formatter
        else:
            raise TypeError("formatter should be a string for Formatter instance")
        self._grid_helper.invalidate()

    def set_ticks(self, values=None, spacing=None, number=None):
        if values is not None:
            self._fl_helper.values = values
        elif spacing is not None:
            self._fl_helper.spacing = spacing
        elif number is not None:
            self._fl_helper.number = number
        else:
            raise ValueError("one of values, spacing, or number should be specified")
        self._grid_helper.invalidate()


class ScalarCoordinateHelper(BaseCoordinateHelper):
    pass

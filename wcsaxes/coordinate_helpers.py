import abc

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
    def set_ticklabel_position(self):
        """
        Set the axes on which the ticklabels for this coordinate should
        appear. Should be a string containing zero or more of ``b``, ``t``,
        ``l``, ``r``.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def set_major_formatter(self):
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

    def set_formatter_precision(self, precision):
        self.tick_formatter.precision = precision

    def set_ticks_position(self):
        raise NotImplementedError()

    def set_ticklabel_position(self):
        raise NotImplementedError()

    def set_major_formatter(self):
        raise NotImplementedError()

    def set_ticks(self, spacing=None, number=None):
        raise NotImplementedError()

    def grid(self):
        raise NotImplementedError()


class ScalarCoordinateHelper(BaseCoordinateHelper):
    pass

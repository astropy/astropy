import numpy as np

from matplotlib.transforms import Affine2D
from axisartist import Axes

from .transforms import WCSWorld2PixelTransform

import axisartist.angle_helper as angle_helper
from axisartist import GridHelperCurveLinear
from matplotlib.ticker import Formatter

from astropy import units as u
from astropy.coordinates import Angle

# Define a formatter based on the Astropy Angle class. Of course, we can add
# formatting options later, this is just a proof of concept.

class AngleFormatter(Formatter):

    def __init__(self, precision=0):
        self.precision = precision

    def __call__(self, axis, other, value, **kwargs):
        if len(value) > 0:
            angles = Angle(value, unit=u.deg)
            string = angles.to_string(precision=self.precision).tolist()
            return string
        else:
            return []


class SkyCoordinate(object):

    def __init__(self):
        self.grid_locator = angle_helper.LocatorDMS(12)
        self.tick_formatter = AngleFormatter(precision=1)

    def set_major_formatter(self, format):
        raise NotImplementedError()

    def set_ticks():
        raise NotImplementedError()


class ScalarCoordinate(object):

    def __init__(self):
        pass

    def set_major_formatter(self, format):
        raise NotImplementedError()

    def set_ticks():
        raise NotImplementedError()


class SkyCoordinates(object):

    def __init__(self, wcs):

        self._wcs = wcs

        # The following code was taken from the demo_floating_axis example in
        # Matplotlib. We make use of helpers in the mpl_toolkits package.

        self._extreme_finder = angle_helper.ExtremeFinderCycle(20, 20,
                                                               lon_cycle = 360,
                                                               lat_cycle = None,
                                                               lon_minmax = (-180., 180),
                                                               lat_minmax = (-90., 90.)
                                                               )

        # Set up coordinates
        self._coords = {}
        self._coords[0] = SkyCoordinate()
        self._coords[1] = SkyCoordinate()

        # Set up aliases for coordinates
        name_1 = self._wcs.wcs.ctype[0][:4]
        self._coords[name_1] = self._coords[0]
        name_2 = self._wcs.wcs.ctype[1][:4]
        self._coords[name_2] = self._coords[1]

        # Set up transform
        self._transform = WCSWorld2PixelTransform(self._wcs)

    @property
    def grid_helper(self):
        return GridHelperCurveLinear(self._transform,
                                     extreme_finder=self._extreme_finder,
                                     grid_locator1=self[0].grid_locator,
                                     grid_locator2=self[1].grid_locator,
                                     tick_formatter1=self[0].tick_formatter,
                                     tick_formatter2=self[1].tick_formatter
                                     )

    def __getitem__(self, item):
        return self._coords[item]


class ScalarCoordinates(object):

    def __init__(self, ):
        raise NotImplementedError()

    def __getitem__(self, item):
        raise NotImplementedError()


class WCSAxes(Axes):

    def __init__(self, fig, rect, wcs=None, adjustable='box'):

        self.wcs = wcs

        # For now, assume WCS is Sky WCS
        self.coords = SkyCoordinates(self.wcs)

        Axes.__init__(self, fig, rect, adjustable=adjustable, grid_helper=self.coords.grid_helper)

    def get_transform(self, frame, equinox=None, obstime=None):

        if self.wcs is None and frame != 'pixel':
            raise ValueError('No WCS specified, so only pixel coordinates are available')

        if frame == 'pixel':

            return Affine2D() + self.transData

        elif frame == 'world':

            return WCSWorld2PixelTransform(self.wcs) + self.transData

        else:

            raise NotImplemented("frame {0} not implemented".format(frame))

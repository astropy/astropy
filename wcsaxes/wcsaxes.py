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


class WCSAxes(Axes):

    def __init__(self, fig, rect, wcs=None, adjustable='box'):

        self.wcs = wcs

        # The following code was taken from the demo_floating_axis example in
        # Matplotlib. We make use of helpers in the mpl_toolkits package.

        extreme_finder = angle_helper.ExtremeFinderCycle(20, 20,
                                                 lon_cycle = 360,
                                                 lat_cycle = None,
                                                 lon_minmax = (-180., 180),
                                                 lat_minmax = (-90., 90.)
                                                 )

        grid_locator1 = angle_helper.LocatorDMS(12)
        grid_locator2 = angle_helper.LocatorDMS(12)

        tick_formatter1 = AngleFormatter(precision=1)
        tick_formatter2 = AngleFormatter(precision=1)

        transform = WCSWorld2PixelTransform(self.wcs)

        grid_helper = GridHelperCurveLinear(transform,
                                            extreme_finder=extreme_finder,
                                            grid_locator1=grid_locator1,
                                            tick_formatter1=tick_formatter1,
                                            tick_formatter2=tick_formatter2
                                            )

        Axes.__init__(self, fig, rect, adjustable=adjustable, grid_helper=grid_helper)

    def get_transform(self, frame, equinox=None, obstime=None):

        if self.wcs is None and frame != 'pixel':
            raise ValueError('No WCS specified, so only pixel coordinates are available')

        if frame == 'pixel':

            return Affine2D() + self.transData

        elif frame == 'world':

            return WCSWorld2PixelTransform(self.wcs) + self.transData

        else:

            raise NotImplemented("frame {0} not implemented".format(frame))

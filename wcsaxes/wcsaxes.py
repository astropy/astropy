from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import Axes

from .transforms import WCSWorld2PixelTransform
from .grid_helpers import SkyGridHelper


class WCSAxes(Axes):

    def __init__(self, fig, rect, wcs=None, adjustable='box'):

        self.wcs = wcs

        # For now, assume WCS is Sky WCS
        self.coords = SkyGridHelper(self, self.wcs)

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

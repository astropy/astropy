from matplotlib.transforms import Affine2D

from .hostaxes import HostAxes
from .transforms import WCSWorld2PixelTransform


class WCSAxes(HostAxes):

    def __init__(self, fig, rect, wcs=None, adjustable='datalim'):

        HostAxes.__init__(self, fig, rect, adjustable=adjustable)

        self.wcs = wcs

    def get_transform(self, frame, equinox=None, obstime=None):

        if self.wcs is None and frame != 'pixel':
            raise ValueError('No WCS specified, so only pixel coordinates are available')

        if frame == 'pixel':

            return Affine2D() + self.transData

        elif frame == 'world':

            return WCSWorld2PixelTransform(self.wcs) + self.transData

        else:

            raise NotImplemented("frame {0} not implemented".format(frame))

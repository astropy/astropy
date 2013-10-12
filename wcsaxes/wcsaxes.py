from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import Axes

from astropy.wcs import WCS

from .transforms import WCSPixel2WorldTransform, WCSWorld2PixelTransform, CoordinateTransform
from .grid_helpers import SkyGridHelper
from .utils import get_coordinate_system


class WCSAxes(Axes):

    def __init__(self, fig, rect, wcs=None, adjustable='box'):

        self.wcs = wcs

        if self.wcs is None:

            Axes.__init__(self, fig, rect, adjustable=adjustable)

        else:

            # For now, assume WCS is Sky WCS
            self.coords = SkyGridHelper(self, self.wcs)

            Axes.__init__(self, fig, rect, adjustable=adjustable, grid_helper=self.coords.grid_helper)

    def get_transform(self, frame, equinox=None, obstime=None):

        if self.wcs is None and frame != 'pixel':
            raise ValueError('No WCS specified, so only pixel coordinates are available')

        if isinstance(frame, WCS):

            coord_in = get_coordinate_system(frame)
            coord_out = get_coordinate_system(self.wcs)

            if coord_in == coord_out:

                return (WCSPixel2WorldTransform(frame)
                        + WCSWorld2PixelTransform(self.wcs)
                        + self.transData)

            else:

                return (WCSPixel2WorldTransform(frame)
                        + CoordinateTransform(coord_in, coord_out)
                        + WCSWorld2PixelTransform(self.wcs)
                        + self.transData)

        elif frame == 'pixel':

            return Affine2D() + self.transData

        else:

            from astropy.coordinates import FK5Coordinates, GalacticCoordinates

            world2pixel = WCSWorld2PixelTransform(self.wcs) + self.transData

            coord_class = get_coordinate_system(self.wcs)

            if frame == 'world':

                return world2pixel

            elif frame == 'fk5':

                if coord_class is FK5Coordinates:
                    return world2pixel
                else:
                    return CoordinateTransform(FK5Coordinates, coord_class) + world2pixel

            elif frame == 'galactic':

                if coord_class is GalacticCoordinates:
                    return world2pixel
                else:
                    return CoordinateTransform(GalacticCoordinates, coord_class) + world2pixel

            else:

                raise NotImplemented("frame {0} not implemented".format(frame))

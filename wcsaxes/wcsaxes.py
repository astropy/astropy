from matplotlib.axes import Axes
from matplotlib.transforms import Affine2D, Bbox

from astropy.wcs import WCS

from .transforms import (WCSPixel2WorldTransform, WCSWorld2PixelTransform,
                         CoordinateTransform)
from .coordinates_map import CoordinatesMap
from .utils import get_coordinate_system
from .coordinate_range import find_coordinate_range

__all__ = ['WCSAxes']

IDENTITY = WCS(naxis=2)
IDENTITY.wcs.ctype = ["X", "Y"]
IDENTITY.wcs.crval = [1., 1.]
IDENTITY.wcs.crpix = [1., 1.]
IDENTITY.wcs.cdelt = [1., 1.]


class WCSAxes(Axes):

    def __init__(self, fig, rect, wcs=IDENTITY, **kwargs):

        self.wcs = wcs

        super(WCSAxes, self).__init__(fig, rect, **kwargs)

        # Turn off spines and current axes

        for s in self.spines.values():
            s.set_visible(False)

        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)

        # Here determine all the coordinate axes that should be shown.

        self.coords = CoordinatesMap(self, self.wcs)

        self._all_coords = [self.coords]

        # Common default settings
        self.coords[0].set_axislabel_position('b')
        self.coords[1].set_axislabel_position('l')
        self.coords[0].set_ticklabel_position('b')
        self.coords[1].set_ticklabel_position('l')

    def get_coord_range(self, transform):
        xmin, xmax = self.get_xlim()
        ymin, ymax = self.get_ylim()
        return find_coordinate_range(transform.inverted(),
                                     [xmin, xmax, ymin, ymax],
                                     x_type=self.coords[0].coord_type,
                                     y_type=self.coords[1].coord_type)

    def draw(self, renderer, inframe=False):

        super(WCSAxes, self).draw(renderer, inframe)

        # Here need to find out range of all coordinates, and update range for
        # each coordinate axis. For now, just assume it covers the whole sky.

        self._bboxes = []

        for coords in self._all_coords:

            coords.frame.update()
            coords[0]._draw(renderer, bboxes=self._bboxes)
            coords[1]._draw(renderer, bboxes=self._bboxes)

        for coords in self._all_coords:

            coords[0]._draw_axislabels(renderer, bboxes=self._bboxes)
            coords[1]._draw_axislabels(renderer, bboxes=self._bboxes)

        self.coords.frame.draw(renderer)

    def get_coords_overlay(self, frame, equinox=None, obstime=None):

        # Here we can't use get_transform because that deals with
        # pixel-to-pixel transformations when passing a WCS object.
        if isinstance(frame, WCS):
            coords = CoordinatesMap(self, frame)
        else:
            transform = self.get_transform(frame, equinox=equinox, obstime=obstime) - self.transData
            coords = CoordinatesMap(self, self.wcs, transform=transform)

        self._all_coords.append(coords)

        # Common settings for overlay
        coords[0].set_axislabel_position('t')
        coords[1].set_axislabel_position('r')
        coords[0].set_ticklabel_position('t')
        coords[1].set_ticklabel_position('r')

        return coords

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

            from astropy.coordinates import FK5, Galactic

            world2pixel = WCSWorld2PixelTransform(self.wcs) + self.transData

            if frame == 'world':

                return world2pixel

            elif frame == 'fk5':

                coord_class = get_coordinate_system(self.wcs)

                if coord_class is FK5:
                    return world2pixel
                else:
                    return (CoordinateTransform(FK5, coord_class)
                            + world2pixel)

            elif frame == 'galactic':

                coord_class = get_coordinate_system(self.wcs)

                if coord_class is Galactic:
                    return world2pixel
                else:
                    return (CoordinateTransform(Galactic, coord_class)
                            + world2pixel)

            else:

                raise NotImplemented("frame {0} not implemented".format(frame))

    def get_tightbbox(self, renderer):

        if not self.get_visible():
            return

        bb = [b for b in self._bboxes if b and (b.width!=0 or b.height!=0)]

        if bb:
            _bbox = Bbox.union(bb)
            return _bbox
        else:
            return []

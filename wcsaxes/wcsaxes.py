from matplotlib.axes import Axes, subplot_class_factory
from matplotlib.transforms import Affine2D, Bbox, Transform

from astropy.wcs import WCS

from .transforms import (WCSPixel2WorldTransform, WCSWorld2PixelTransform,
                         CoordinateTransform)
from .coordinates_map import CoordinatesMap
from .utils import get_coordinate_system
from .coordinate_range import find_coordinate_range

__all__ = ['WCSAxes', 'WCSAxesSubplot']

IDENTITY = WCS(naxis=2)
IDENTITY.wcs.ctype = ["X", "Y"]
IDENTITY.wcs.crval = [1., 1.]
IDENTITY.wcs.crpix = [1., 1.]
IDENTITY.wcs.cdelt = [1., 1.]


class WCSAxes(Axes):

    def __init__(self, fig, rect, wcs=IDENTITY, transData=None,
                 **kwargs):
        super(WCSAxes, self).__init__(fig, rect, **kwargs)
        self._bboxes = []

        self.reset_wcs(wcs)
        self._hide_parent_artists()

        if not (transData is None):
            # User wants to override the transform for the final
            # data->pixel mapping
            self.transData = transData
        
    def _hide_parent_artists(self):
        # Turn off spines and current axes
        for s in self.spines.values():
            s.set_visible(False)

        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)

    def reset_wcs(self, wcs):
        """
        Reset the current Axes, to use a new WCS object.
        """

        # Here determine all the coordinate axes that should be shown.
        if wcs is None:
            wcs = IDENTITY

        self.wcs = wcs
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

    def set_xlabel(self, label):
        self.coords[0].set_axislabel(label)

    def set_ylabel(self, label):
        self.coords[1].set_axislabel(label)

    def get_coords_overlay(self, frame, equinox=None, obstime=None):

        # Here we can't use get_transform because that deals with
        # pixel-to-pixel transformations when passing a WCS object.
        if isinstance(frame, WCS):
            coords = CoordinatesMap(self, frame)
        else:
            transform = self._get_transform_no_transdata(frame, equinox=equinox, obstime=obstime)
            coords = CoordinatesMap(self, self.wcs, transform=transform)

        self._all_coords.append(coords)

        # Common settings for overlay
        coords[0].set_axislabel_position('t')
        coords[1].set_axislabel_position('r')
        coords[0].set_ticklabel_position('t')
        coords[1].set_ticklabel_position('r')

        return coords

    def get_transform(self, frame, equinox=None, obstime=None):
        """
        Return a transform from the specified frame to display coordinates.

        This does not include the transData transformation

        Parameters
        ----------
        frame : :class:`~astropy.wcs.WCS` or :class:`~matplotlib.transforms.Transform` or str
            The ``frame`` parameter can have several possible types:
                * :class:`~astropy.wcs.WCS` instance: assumed to be a
                  transformation from pixel to world coordinates, where the
                  world coordinates are the same as those in the WCS
                  transformation used for this ``WCSAxes`` instance. This is
                  used for example to show contours, since this involves
                  plotting an array in pixel coordinates that are not the
                  final data coordinate and have to be transformed to the
                  common world coordinate system first.
                * :class:`~matplotlib.transforms.Transform` instance: it is
                  assumed to be a transform to the world coordinates that are
                  part of the WCS used to instantiate this ``WCSAxes``
                  instance.
                * ``'pixel'`` or ``'world'``: return a transformation that
                  allows users to plot in pixel/data coordinates (essentially
                  an identity transform) and ``world`` (the default
                  world-to-pixel transformation used to instantiate the
                  ``WCSAxes`` instance).
                * ``'fk5'`` or ``'galactic'``: return a transformation from
                  the specified frame to the pixel/data coordinates.
        """
        return self._get_transform_no_transdata(frame, equinox=equinox, obstime=obstime) + self.transData

    def _get_transform_no_transdata(self, frame, equinox=None, obstime=None):
        """
        Return a transform from the specified frame to data coordinates.

        As for :meth:`wcsaxes.wcsaxes.WCSAxes.get_transform` but does not
        include the last transformation from data to display coordinates.
        """

        if self.wcs is None and frame != 'pixel':
            raise ValueError('No WCS specified, so only pixel coordinates are available')

        if isinstance(frame, WCS):

            coord_in = get_coordinate_system(frame)
            coord_out = get_coordinate_system(self.wcs)

            if coord_in == coord_out:

                return (WCSPixel2WorldTransform(frame)
                        + WCSWorld2PixelTransform(self.wcs))

            else:

                return (WCSPixel2WorldTransform(frame)
                        + CoordinateTransform(coord_in, coord_out)
                        + WCSWorld2PixelTransform(self.wcs))

        elif frame == 'pixel':

            return Affine2D()

        elif isinstance(frame, Transform):

            world2pixel = WCSWorld2PixelTransform(self.wcs)

            return frame + world2pixel

        else:

            from astropy.coordinates import FK5, Galactic

            world2pixel = WCSWorld2PixelTransform(self.wcs)

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

        bb = [b for b in self._bboxes if b and (b.width != 0 or b.height != 0)]

        if bb:
            _bbox = Bbox.union(bb)
            return _bbox
        else:
            return self.get_window_extent(renderer)

    def grid(self, draw_grid=True, **kwargs):
        """
        Plot gridlines for both coordinates.

        Standard matplotlib appearance options (color, alpha, etc.) can be
        passed as keyword arguments.

        Parameters
        ----------
        draw_grid : bool
            Whether to show the gridlines
        """
        if draw_grid:
            self.coords.grid(draw_grid=draw_grid, **kwargs)


WCSAxesSubplot = subplot_class_factory(WCSAxes)

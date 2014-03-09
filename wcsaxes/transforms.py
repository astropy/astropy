# Note: This file incldues code dervived from pywcsgrid2
#
# This file contains Matplotlib transformation objects (e.g. from pixel to world
# coordinates, but also world-to-world).

import abc

import numpy as np
from matplotlib.path import Path
from matplotlib.transforms import Transform
from astropy import units as u


class CurvedTransform(Transform):
    """
    Abstract base class for non-affine curved transforms
    """

    __metaclass__ = abc.ABCMeta

    input_dims = 2
    output_dims = 2
    is_separable = False

    def transform_path(self, path):
        """
        Transform a Matplotlib Path

        Parameters
        ----------
        path : :class:`~matplotlib.path.Path`
            The path to transform

        Returns
        -------
        path : :class:`~matplotlib.path.Path`
            The resulting path
        """
        return Path(self.transform(path.vertices), path.codes)

    transform_path_non_affine = transform_path

    @abc.abstractmethod
    def transform(self, input):
        raise NotImplemented("")

    @abc.abstractmethod
    def inverted(self):
        raise NotImplemented("")


class WCSWorld2PixelTransform(CurvedTransform):

    """
    WCS transformation from world to pixel coordinates
    """

    def __init__(self, wcs):
        super(WCSWorld2PixelTransform, self).__init__()
        self.wcs = wcs

    def transform(self, world):
        """
        Transform world to pixel coordinates
        """

        xw, yw = world[:, 0], world[:, 1]

        xp, yp = self.wcs.wcs_world2pix(xw, yw, 1)
        xp, yp = xp - 1, yp - 1
        pixel = np.concatenate((xp[:, np.newaxis], yp[:, np.newaxis]), 1)

        return pixel

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return WCSPixel2WorldTransform(self.wcs)


class WCSPixel2WorldTransform(CurvedTransform):

    """
    WCS transformation from pixel to world coordinates
    """

    def __init__(self, wcs):
        super(WCSPixel2WorldTransform, self).__init__()
        self.wcs = wcs

    def transform(self, pixel):
        """
        Transform pixel to world coordinates
        """

        xp, yp = pixel[:, 0], pixel[:, 1]

        xp, yp = xp + 1, yp + 1
        xw, yw = self.wcs.wcs_pix2world(xp, yp, 1)

        # At the moment, one has to manually check that the transformation
        # round-trips, otherwise it should be considered invalid.
        xp_check, yp_check = self.wcs.wcs_world2pix(xw, yw, 1)
        invalid = ((np.abs(xp_check - xp) > 1.) |
                   (np.abs(yp_check - yp) > 1.))
        xw[invalid] = np.nan
        yw[invalid] = np.nan

        world = np.concatenate((xw[:, np.newaxis], yw[:, np.newaxis]), 1)

        return world

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return WCSWorld2PixelTransform(self.wcs)


class CoordinateTransform(CurvedTransform):

    def __init__(self, input_system, output_system):
        super(CoordinateTransform, self).__init__()
        self.input_system = input_system
        self.output_system = output_system

    def transform(self, input_coords):
        """
        Transform one set of coordinates to another
        """

        x_in, y_in = input_coords[:, 0], input_coords[:, 1]

        c_in = self.input_system(x_in, y_in, unit=(u.deg, u.deg))

        c_out = c_in.transform_to(self.output_system)

        return np.concatenate((c_out.lonangle.deg[:, np.newaxis], c_out.latangle.deg[:, np.newaxis]), 1)

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return CoordinateTransform(self.output_system, self.input_system)

# Note: This file incldues code dervived from pywcsgrid2

import abc

import numpy as np
from matplotlib.path import Path
from matplotlib.transforms import Transform


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
        world = np.concatenate((xw[:, np.newaxis], yw[:, np.newaxis]), 1)

        return world

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform
        """
        return WCSWorld2PixelTransform(self.wcs)

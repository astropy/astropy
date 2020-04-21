import numpy as np

from .base import BaseWCSWrapper

__all__ = ['ResampledLowLevelWCS']


class ResampledLowLevelWCS(BaseWCSWrapper):
    """
    A wrapper for a low-level WCS object that has down- or
    up-sampled pixel axes.

    Parameters
    ----------
    wcs : `~astropy.wcs.wcsapi.BaseLowLevelWCS`
        The original WCS for which to reorder axes
    factor : int or float or iterable
        The factor by which to increase the pixel size for each pixel
        axis. If a scalar, the same factor is used for all axes.
    """
    def __init__(self, wcs, factor):
        self._wcs = wcs
        if np.isscalar(factor):
            factor = [factor] * self.pixel_n_dim
        self._factor = factor

    def pixel_to_world_values(self, *pixel_arrays):
        pixel_arrays = [np.asarray(pixel_arrays[i]) * self._factor[i]
                        for i in range(self.pixel_n_dim)]
        return self._wcs.pixel_to_world_values(*pixel_arrays)

    def world_to_pixel_values(self, *world_arrays):
        pixel_arrays = self._wcs.world_to_pixel_values(*world_arrays)
        pixel_arrays = [np.asarray(pixel_arrays[i]) / self._factor[i]
                        for i in range(self.pixel_n_dim)]
        return pixel_arrays

    @property
    def pixel_shape(self):
        return tuple(self._wcs.pixel_shape[i] / self._factor[i]
                     for i in range(self.pixel_n_dim))

    @property
    def pixel_bounds(self):
        return tuple((self._wcs.pixel_bounds[i][0] / self._factor[i],
                      self._wcs.pixel_bounds[i][1] / self._factor[i])
                     for i in range(self.pixel_n_dim))

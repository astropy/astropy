import numpy as np

from astropy.wcs.wcsapi import BaseLowLevelWCS, wcs_info_str

__all__ = ['ResampledLowLevelWCS']


class ResampledLowLevelWCS(BaseLowLevelWCS):
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

    @property
    def pixel_n_dim(self):
        return self._wcs.pixel_n_dim

    @property
    def world_n_dim(self):
        return self._wcs.world_n_dim

    @property
    def world_axis_physical_types(self):
        return self._wcs.world_axis_physical_types

    @property
    def world_axis_units(self):
        return self._wcs.world_axis_units

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
    def world_axis_object_components(self):
        return self._wcs.world_axis_object_components

    @property
    def world_axis_object_classes(self):
        return self._wcs.world_axis_object_classes

    @property
    def pixel_shape(self):
        return tuple(self._wcs.pixel_shape[i] / self._factor[i]
                     for i in range(self.pixel_n_dim))

    @property
    def pixel_bounds(self):
        return tuple((self._wcs.pixel_bounds[i][0] / self._factor[i],
                      self._wcs.pixel_bounds[i][1] / self._factor[i])
                     for i in range(self.pixel_n_dim))

    @property
    def pixel_axis_names(self):
        return self._wcs.pixel_axis_names

    @property
    def world_axis_names(self):
        return self._wcs.world_axis_names

    @property
    def axis_correlation_matrix(self):
        return self._wcs.axis_correlation_matrix

    @property
    def serialized_classes(self):
        return self._wcs.serialized_classes

    def __repr__(self):
        return wcs_info_str(self)

    def __str__(self):
        return wcs_info_str(self)

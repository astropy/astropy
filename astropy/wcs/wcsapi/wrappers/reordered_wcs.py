
import numbers

import numpy as np

from astropy.wcs.wcsapi import BaseLowLevelWCS, wcs_info_str
from astropy.utils import isiterable

__all__ = ['ReorderedLowLevelWCS']


class ReorderedLowLevelWCS(BaseLowLevelWCS):
    """
    A wrapper for a low-level WCS object that has re-ordered
    pixel and/or world axes.

    Parameters
    ----------
    wcs : `~astropy.wcs.wcsapi.BaseLowLevelWCS`
        The original WCS for which to reorder axes
    pixel_order : iterable
        The indices of the original axes in the order of the
        new WCS.
    world_order : iterable
        The indices of the original axes in the order of the
        new WCS.
    """
    def __init__(self, wcs, pixel_order, world_order):
        self._wcs = wcs
        self._pixel_order = pixel_order
        self._world_order = world_order
        self._pixel_order_inv = np.argsort(pixel_order)
        self._world_order_inv = np.argsort(world_order)

    @property
    def pixel_n_dim(self):
        return self._wcs.pixel_n_dim

    @property
    def world_n_dim(self):
        return self._wcs.world_n_dim

    @property
    def world_axis_physical_types(self):
        return [self._wcs.world_axis_physical_types[idx] for idx in self._world_order]

    @property
    def world_axis_units(self):
        return [self._wcs.world_axis_units[idx] for idx in self._world_order]

    @property
    def pixel_axis_names(self):
        return [self._wcs.pixel_axis_names[idx] for idx in self._pixel_order]

    @property
    def world_axis_names(self):
        return [self._wcs.world_axis_names[idx] for idx in self._world_order]

    def pixel_to_world_values(self, *pixel_arrays):
        pixel_arrays = [pixel_arrays[idx] for idx in self._pixel_order_inv]
        world_arrays = self._wcs.pixel_to_world_values(*pixel_arrays)
        world_arrays = [world_arrays[idx] for idx in self._world_order]
        return world_arrays

    def world_to_pixel_values(self, *world_arrays):
        world_arrays = [world_arrays[idx] for idx in self._world_order_inv]
        pixel_arrays = self._wcs.world_to_pixel_values(*world_arrays)
        pixel_arrays = [pixel_arrays[idx] for idx in self._pixel_order]
        return pixel_arrays

    @property
    def world_axis_object_components(self):
        return [self._wcs.world_axis_object_components[idx] for idx in self._world_order]

    @property
    def world_axis_object_classes(self):
        return self._wcs.world_axis_object_classes

    @property
    def array_shape(self):
        if self.pixel_shape:
            return self.pixel_shape[::-1]

    @property
    def pixel_shape(self):
        if self._wcs.pixel_shape:
            return tuple([self._wcs.pixel_shape[idx] for idx in self._pixel_order])

    @property
    def pixel_bounds(self):
        if self._wcs.pixel_bounds:
            return tuple([self._wcs.pixel_bounds[idx] for idx in self._pixel_order])

    @property
    def axis_correlation_matrix(self):
        return self._wcs.axis_correlation_matrix[self._world_order][:, self._pixel_order]

    def __repr__(self):
        return wcs_info_str(self)

    def __str__(self):
        return wcs_info_str(self)

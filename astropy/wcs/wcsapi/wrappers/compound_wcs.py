from functools import reduce

import numpy as np

from astropy.wcs.wcsapi import BaseLowLevelWCS, wcs_info_str

__all__ = ['CompoundLowLevelWCS']


def tuplesum(lists):
    return reduce(tuple.__add__, map(tuple, lists))


class CompoundLowLevelWCS(BaseLowLevelWCS):
    """
    A wrapper that takes multiple low level WCS objects and makes a compound
    WCS that combines them.

    Parameters
    ----------
    *wcs : `~astropy.wcs.wcsapi.BaseLowLevelWCS`
        The WCSes to combine
    """
    def __init__(self, *wcs):
        self._wcs = wcs

    @property
    def pixel_n_dim(self):
        return sum([w.pixel_n_dim for w in self._wcs])

    @property
    def world_n_dim(self):
        return sum([w.pixel_n_dim for w in self._wcs])

    @property
    def world_axis_physical_types(self):
        return tuplesum([w.world_axis_physical_types for w in self._wcs])

    @property
    def world_axis_units(self):
        return tuplesum([w.world_axis_units for w in self._wcs])

    def pixel_to_world_values(self, *pixel_arrays):
        world_arrays = []
        for w in self._wcs:
            pixel_arrays_sub = pixel_arrays[:w.pixel_n_dim]
            pixel_arrays = pixel_arrays[w.pixel_n_dim:]
            world_arrays_sub = w.pixel_to_world_values(*pixel_arrays_sub)
            if w.world_n_dim > 1:
                world_arrays.extend(world_arrays_sub)
            else:
                world_arrays.append(world_arrays_sub)
        return tuple(world_arrays)

    def world_to_pixel_values(self, *world_arrays):
        pixel_arrays = []
        for w in self._wcs:
            world_arrays_sub = world_arrays[:w.world_n_dim]
            world_arrays = world_arrays[w.world_n_dim:]
            pixel_arrays_sub = w.world_to_pixel_values(*world_arrays_sub)
            if w.pixel_n_dim > 1:
                pixel_arrays.extend(pixel_arrays_sub)
            else:
                pixel_arrays.append(pixel_arrays_sub)
        return tuple(pixel_arrays)

    @property
    def world_axis_object_components(self):
        all_components = []
        for iw, w in enumerate(self._wcs):
            for component in w.world_axis_object_components:
                all_components.append((f'{component[0]}_{iw}',) + component[1:])
        return all_components

    @property
    def world_axis_object_classes(self):
        # TODO: deal with name conflicts
        all_classes = {}
        for iw, w in enumerate(self._wcs):
            for key, value in w.world_axis_object_classes.items():
                all_classes[f'{key}_{iw}'] = value
        return all_classes

    @property
    def pixel_shape(self):
        if any(w.array_shape is None for w in self._wcs):
            return None
        else:
            return tuplesum(w.pixel_shape for w in self._wcs)

    @property
    def pixel_bounds(self):
        if any(w.pixel_bounds is None for w in self._wcs):
            return None
        else:
            return tuplesum(w.pixel_bounds for w in self._wcs)

    @property
    def pixel_axis_names(self):
        return tuplesum(w.pixel_axis_names for w in self._wcs)

    @property
    def world_axis_names(self):
        return tuplesum(w.world_axis_names for w in self._wcs)

    @property
    def axis_correlation_matrix(self):
        matrix = np.zeros((self.world_n_dim, self.pixel_n_dim), dtype=bool)
        iw = ip = 0
        for w in self._wcs:
            matrix[iw:iw+w.world_n_dim,ip:ip+w.pixel_n_dim] = w.axis_correlation_matrix
            iw += w.world_n_dim
            ip += w.pixel_n_dim
        return matrix

    @property
    def serialized_classes(self):
        return any([w.serialized_classes for w in self._wcs])

    def __repr__(self):
        return wcs_info_str(self)

    def __str__(self):
        return wcs_info_str(self)

import numpy as np
from .low_level_api import BaseLowLevelWCS

__all__ = ['CompoundLowLevelWCS']


from functools import reduce

def listsum(lists):
    return reduce(list.__add__, map(list, lists))


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
        return listsum([w.world_axis_physical_types for w in self._wcs])

    @property
    def world_axis_units(self):
        return listsum([w.world_axis_units for w in self._wcs])

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

    def array_index_to_world_values(self, *index_arrays):
        return self.pixel_to_world_values(*index_arrays[-1])

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
        return tuple(world_arrays)

    def world_to_array_index_values(self, *world_arrays):
        pixel_arrays = self.world_to_pixel_values(*world_arrays)
        array_indices = tuple(np.asarray(np.floor(pixel + 0.5), dtype=np.int_) for pixel in pixel_arrays)
        return array_indices

    @property
    def world_axis_object_components(self):
        return listsum([w.world_axis_object_components for w in self._wcs])

    @property
    def world_axis_object_classes(self):
        # TODO: deal with name conflicts
        all_classes = {}
        for w in self._wcs:
            all_classes.update(w.world_axis_object_components)
        return all_classes

    @property
    def array_shape(self):
        if any(w.array_shape is None for w in self._wcs):
            return None
        else:
            return listsum(w.array_shape for w in self._wcs)

    @property
    def pixel_shape(self):
        if self.array_shape is None:
            return None
        else:
            return self.array_shape[::-1]

    @property
    def pixel_bounds(self):
        if any(w.pixel_bounds is None for w in self._wcs):
            return None
        else:
            return listsum(w.pixel_bounds for w in self._wcs)

    @property
    def pixel_axis_names(self):
        return listsum(w.pixel_axis_names for w in self._wcs)

    @property
    def world_axis_names(self):
        return listsum(w.world_axis_names for w in self._wcs)

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

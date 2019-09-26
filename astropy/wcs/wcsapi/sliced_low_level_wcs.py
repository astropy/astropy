import numbers
import numpy as np
from astropy.wcs.wcsapi import BaseLowLevelWCS

__all__ = ['sanitize_slices', 'SlicedLowLevelWCS']


def sanitize_slices(slices, ndim):
    """
    Given a set of input
    """

    if not isinstance(slices, (tuple, list)):  # We just have a single int
        slices = (slices,)

    slices = list(slices)

    if Ellipsis in slices:
        if slices.count(Ellipsis) > 1:
            raise IndexError("an index can only have a single ellipsis ('...')")

        # Replace the Ellipsis with the correct number of slice(None)s
        e_ind = slices.index(Ellipsis)
        slices.remove(Ellipsis)
        n_e = ndim - len(slices)
        for i in range(n_e):
            ind = e_ind + i
            slices.insert(ind, slice(None))

    for i in range(ndim):
        if i < len(slices):
            slc = slices[i]
            if isinstance(slc, slice):
                if slc.step and slc.step != 1:
                    raise ValueError("Slicing WCS with a step is not supported.")
            elif not isinstance(slc, numbers.Integral):
                raise ValueError("Only integer or range slices are accepted.")
        else:
            slices.append(slice(None))

    return slices


class SlicedLowLevelWCS(BaseLowLevelWCS):

    def __init__(self, wcs, slices):

        self._wcs = wcs
        self._slices_array = sanitize_slices(slices, self._wcs.pixel_n_dim)
        self._slices_pixel = self._slices_array[::-1]

        # figure out which pixel dimensions have been kept, then use axis correlation
        # matrix to figure out which world dims are kept
        self._pixel_keep = np.nonzero([not isinstance(self._slices_pixel[ip], numbers.Integral)
                                       for ip in range(self._wcs.pixel_n_dim)])[0]

        # axis_correlation_matrix[world, pixel]
        self._world_keep = np.nonzero(
            self._wcs.axis_correlation_matrix[:, self._pixel_keep].any(axis=1))[0]

    @property
    def pixel_n_dim(self):
        return len(self._pixel_keep)

    @property
    def world_n_dim(self):
        return len(self._world_keep)

    @property
    def world_axis_physical_types(self):
        return [self._wcs.world_axis_physical_types[i] for i in self._world_keep]

    @property
    def world_axis_units(self):
        return [self._wcs.world_axis_units[i] for i in self._world_keep]

    def pixel_to_world_values(self, *pixel_arrays):
        pixel_arrays_new = []
        ipix_curr = -1
        for ipix in range(self._wcs.pixel_n_dim):
            if isinstance(self._slices_pixel[ipix], int):
                pixel_arrays_new.append(self._slices_pixel[ipix])
            else:
                ipix_curr += 1
                if self._slices_pixel[ipix].start is not None:
                    pixel_arrays_new.append(pixel_arrays[ipix_curr] + self._slices_pixel[ipix].start)
                else:
                    pixel_arrays_new.append(pixel_arrays[ipix_curr])

        pixel_arrays_new = np.broadcast_arrays(*pixel_arrays_new)
        world_arrays = self._wcs.pixel_to_world_values(*pixel_arrays_new)
        return [world_arrays[iw] for iw in self._world_keep]

    def array_index_to_world_values(self, *index_arrays):
        return self.pixel_to_world_values(*index_arrays[::-1])

    def world_to_pixel_values(self, *world_arrays):
        world_arrays_new = []
        iworld_curr = -1
        for iworld in range(self._wcs.world_n_dim):
            if iworld in self._world_keep:
                iworld_curr += 1
                world_arrays_new.append(world_arrays[iworld_curr])
            else:
                world_arrays_new.append(1.)

        world_arrays_new = np.broadcast_arrays(*world_arrays_new)
        pixel_arrays = list(self._wcs.world_to_pixel_values(*world_arrays_new))

        for ipixel in range(self._wcs.pixel_n_dim):
            if isinstance(self._slices_pixel[ipixel], slice) and self._slices_pixel[ipixel].start is not None:
                pixel_arrays[ipixel] -= self._slices_pixel[ipixel].start

        return [pixel_arrays[ip] for ip in self._pixel_keep]

    def world_to_array_index_values(self, *world_arrays):
        pixel_arrays = self.world_to_pixel_values(*world_arrays, 0)[::-1]
        array_indices = tuple(np.asarray(np.floor(pixel + 0.5), dtype=np.int) for pixel in pixel_arrays)
        return array_indices

    @property
    def world_axis_object_components(self):
        return [self._wcs.world_axis_object_components[idx] for idx in self._world_keep]

    @property
    def world_axis_object_classes(self):
        keys_keep = [item[0] for item in self.world_axis_object_components]
        return dict([item for item in self._wcs.world_axis_object_classes.items() if item[0] in keys_keep])

    @property
    def array_shape(self):
        if self._wcs.array_shape:
            return np.broadcast_to(0, self._wcs.array_shape)[tuple(self._slices_array)].shape

    @property
    def pixel_shape(self):
        if self.array_shape:
            return self.array_shape[::-1]

    @property
    def pixel_bounds(self):

        if self._wcs.pixel_bounds is None:
            return None

        bounds = []
        for idx in self._pixel_keep:
            if self._slices_pixel[idx].start is None:
                bounds.append(self._wcs.pixel_bounds[idx])
            else:
                imin, imax = self._wcs.pixel_bounds[idx]
                start = self._slices_pixel[idx].start
                bounds.append((imin - start, imax - start))

        return bounds

    @property
    def axis_correlation_matrix(self):
        return self._wcs.axis_correlation_matrix[self._world_keep][:,self._pixel_keep]

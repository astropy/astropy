from functools import reduce

import numpy as np

from .base import BaseWCSWrapper

__all__ = ['CompoundLowLevelWCS']


def tuplesum(lists):
    return reduce(tuple.__add__, map(tuple, lists))


class Mapping:
    """
    Allows inputs to be reordered, duplicated or dropped.

    This is a very stripped down version of `astropy.modeling.models.Mapping`
    to be able to handle input of arbitrary type.

    Parameters
    ----------
    mapping : tuple
        A tuple of integers representing indices of the inputs to this model
        to return and in what order to return them.  See
        :ref:`compound-model-mappings` for more details.

    """
    def __init__(self, mapping):
        self.mapping = mapping
        self.n_inputs = max(mapping) + 1
        self.n_outputs = len(mapping)

    def __call__(self, *values):
        return tuple(values[idx] for idx in self.mapping)

    @property
    def inverse(self):
        mapping = tuple(self.mapping.index(idx)
                        for idx in range(self.n_inputs))
        return type(self)(mapping)

    def __repr__(self):
        return f'<Mapping({self.mapping})>'


class CompoundLowLevelWCS(BaseWCSWrapper):
    """
    A wrapper that takes multiple low level WCS objects and makes a compound
    WCS that combines them.

    Parameters
    ----------
    *wcs : `~astropy.wcs.wcsapi.BaseLowLevelWCS`
        The WCSes to combine
    mapping : `tuple`
        The pixel dimension mapping between the input pixel dimensions and the
        input pixel dimensions to the underlying WCSes. This should have length
        equal to the total number of pixel dimensions in all input WCSes and
        have a maximum of the number of input pixel dimensions to the resulting
        compound WCS -1 (counts from 0). For example ``(0, 1, 2, 1)`` would end
        up with the second and fourth pixel dimensions in the input WCSes being
        shared, so the compound WCS would have 3 pixel dimensions ``(2 + 1)``.
        See :ref:`compound-model-mappings` for more examples of this input
        format.
    pixel_atol : `float`
        A tolerance used to check that the resulting pixel coordinates from
        ``world_to_pixel`` are the same from all WCSes.
    """
    def __init__(self, *wcs, mapping=None, pixel_atol=1e-8):
        self._wcs = wcs

        if not mapping:
            mapping = tuple(range(self._all_pixel_n_dim))

        if not len(mapping) == self._all_pixel_n_dim:
            raise ValueError(
                "The length of the mapping must equal the total number of pixel dimensions in all input WCSes.")

        self.mapping = Mapping(mapping)
        self.atol = pixel_atol

        # Validate the pixel bounds and shape are consistent
        self.pixel_bounds
        self.pixel_shape

    @property
    def _all_pixel_n_dim(self):
        return sum([w.pixel_n_dim for w in self._wcs])

    @property
    def pixel_n_dim(self):
        return self.mapping.n_inputs

    @property
    def world_n_dim(self):
        return sum([w.world_n_dim for w in self._wcs])

    @property
    def world_axis_physical_types(self):
        return tuplesum([w.world_axis_physical_types for w in self._wcs])

    @property
    def world_axis_units(self):
        return tuplesum([w.world_axis_units for w in self._wcs])

    def pixel_to_world_values(self, *pixel_arrays):
        pixel_arrays = self.mapping(*pixel_arrays)
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

        pixel_arrays = tuple(pixel_arrays)
        for i, ix in enumerate(self.mapping.mapping):
            if not np.allclose(pixel_arrays[ix], pixel_arrays[i], atol=self.atol):
                raise ValueError(
                    f"The world inputs for shared pixel axes did not result in a pixel coordinate to within {self.atol} relative accuracy."
                )
        return self.mapping.inverse(*pixel_arrays)

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
        if not any(w.array_shape is None for w in self._wcs):
            pixel_shape = tuplesum(w.pixel_shape for w in self._wcs)
            out_shape = self.mapping.inverse(*pixel_shape)
            for i, ix in enumerate(self.mapping.mapping):
                if out_shape[ix] != pixel_shape[i]:
                    raise ValueError(
                        "The pixel shapes of the supplied WCSes do not match for the dimensions shared by the supplied mapping.")
            return out_shape

    @property
    def pixel_bounds(self):
        if not any(w.pixel_bounds is None for w in self._wcs):
            pixel_bounds = tuplesum(w.pixel_bounds for w in self._wcs)
            out_bounds = self.mapping.inverse(*pixel_bounds)
            for i, ix in enumerate(self.mapping.mapping):
                if out_bounds[ix] != pixel_bounds[i]:
                    raise ValueError(
                        "The pixel bounds of the supplied WCSes do not match for the dimensions shared by the supplied mapping.")
            return out_bounds

    @property
    def pixel_axis_names(self):
        pixel_names = tuplesum(w.pixel_axis_names for w in self._wcs)
        out_names = self.mapping.inverse(*pixel_names)

        for i, ix in enumerate(self.mapping.mapping):
            if out_names[ix] != pixel_names[i]:
                out_names[ix] = ' / '.join([out_names[ix], pixel_names[i]])

        return out_names

    @property
    def world_axis_names(self):
        return tuplesum(w.world_axis_names for w in self._wcs)

    @property
    def axis_correlation_matrix(self):
        full_matrix = np.zeros((self.world_n_dim, self._all_pixel_n_dim), dtype=bool)
        iw = ip = 0
        for w in self._wcs:
            full_matrix[iw:iw + w.world_n_dim, ip:ip + w.pixel_n_dim] = w.axis_correlation_matrix
            iw += w.world_n_dim
            ip += w.pixel_n_dim

        if self._all_pixel_n_dim == self.pixel_n_dim:
            return full_matrix

        matrix = np.zeros((self.world_n_dim, self.pixel_n_dim), dtype=bool)
        for i, ix in enumerate(self.mapping.mapping):
            matrix[:, ix] = np.logical_or(matrix[:, ix], full_matrix[:, i])

        return matrix

    @property
    def serialized_classes(self):
        return any([w.serialized_classes for w in self._wcs])

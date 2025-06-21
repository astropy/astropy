import numpy as np
from astropy.wcs import WCS
from astropy.wcs.wcsapi import HighLevelWCSMixin, BaseLowLevelWCS
from astropy import units as u


class ExactHeaderWCS(BaseLowLevelWCS, HighLevelWCSMixin):
    """
    This is an APE-14-compliant WCS class which takes a Header
    and will respect the units of the header.

    The only public interface to the projection parameters is
    via the header property of the WCS.
    """

    def __init__(self, header):
        self.header = header
        self._wcs = WCS(header)
        self._determine_unit_scales()

    def _determine_unit_scales(self):
        # For now, assume all unit changes are simply scales
        scales = []
        for i, unit in enumerate(self._wcs.wcs.cunit):
            original_unit = self.header.get(f'CUNIT{i+1}', '')
            scales.append(u.Unit(original_unit).to(unit))
        self._scales = scales

    # Properties/methods that need to take into account different units

    @property
    def world_axis_units(self):
        return self._wcs.world_axis_units

    def pixel_to_world_values(self, *pixel_arrays):
        world_arrays = self._wcs.pixel_to_world_values(*pixel_arrays)

        # Detect the case of a length 0 array
        if isinstance(world_arrays, np.ndarray) and not world_arrays.shape:
            return world_arrays / self._scales[0]

        world_arrays = [w / s for (w, s) in zip(world_arrays, self._scales)]

        return world_arrays

    def world_to_pixel_values(self, *world_arrays):
        world_arrays = [w * s for (w, s) in zip(world_arrays, self._scales)]
        pixel_arrays = self._wcs.world_to_pixel_values(*world_arrays)
        return pixel_arrays

    @property
    def world_axis_object_classes(self):
        # TODO!
        return self._wcs.world_axis_object_classes

    # Properties that are plain wrappers of underlying WCS

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
    def pixel_axis_names(self):
        return self._wcs.pixel_axis_names

    @property
    def world_axis_names(self):
        return self._wcs.world_axis_names

    @property
    def array_shape(self):
        return self._wcs.array_shape

    @property
    def pixel_shape(self):
        return self._wcs.pixel_shape

    @property
    def pixel_bounds(self):
        return self._wcs.pixel_bounds

    @property
    def axis_correlation_matrix(self):
        return self._wcs.axis_correlation_matrix

    @property
    def world_axis_object_components(self):
        # TODO
        return self._wcs.world_axis_object_components

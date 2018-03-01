from .base_low_level_api import BaseLowLevelWCS

__all__ = ['FITSWCSLowLevelWCS']


class FITSWCSLowLevelWCS(BaseLowLevelWCS):
    """
    A wrapper around the :class:`astropy.wcs.WCS` class that provides the
    low-level WCS API from APE 14.
    """

    def __init__(self, wcs):
        self._wcs = wcs

    @property
    def pixel_n_dim(self):
        return self._wcs.naxis

    @property
    def world_n_dim(self):
        return len(self._wcs.wcs.ctype)

    @property
    def pixel_shape(self):
        return None

    @property
    def pixel_bounds(self):
        return None

    @property
    def world_axis_physical_types(self):
        raise NotImplementedError()

    @property
    def world_axis_units(self):
        return self._wcs.wcs.cunit

    @property
    def axis_correlation_matrix(self):

        # If there are any distortions present, we assume that there may be
        # correlations between all axes. Maybe if some distortions only apply
        # to the image plane we can improve this
        for distortion_attribute in ('sip', 'det2im1', 'det2im2'):
            if getattr(self._wcs, distortion_attribute):
                return np.ones((self.n_world, self.n_pixel), dtype=bool)

        # Assuming linear world coordinates along each axis, the correlation
        # matrix would be given by whether or not the PC matrix is zero
        matrix = self._wcs.wcs.get_pc() != 0

        # We now need to check specifically for celestial coordinates since
        # these can assume correlations because of spherical distortions. For
        # each celestial coordinate we copy over the pixel dependencies from
        # the other celestial coordinates.
        celestial = (self._wcs.wcs.axis_types // 1000) % 10 == 2
        celestial_indices = np.nonzero(celestial)[0]
        for world1 in celestial_indices:
            for world2 in celestial_indices:
                if world1 != world2:
                    matrix[world1] |= matrix[world2]
                    matrix[world2] |= matrix[world1]

        return matrix

    def pixel_to_world_values(self, *pixel_arrays):
        return self.all_pixel_to_world(*pixel_arrays, 0)

    def world_to_pixel_values(self, *world_arrays):
        return self.all_world_to_pixel(*world_arrays, 0)

    @property
    def world_axis_object_components(self):
        raise NotImplementedError()

    @property
    def world_axis_object_classes(self):
        raise NotImplementedError()

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

        # We flip the PC matrix with [::-1] because WCS and numpy index conventions
        # are reversed.
        pc = np.array(self._wcs.wcs.get_pc()[::-1, ::-1])
        axes = self._wcs.get_axis_types()[::-1]

        # There might be smarter ways to do this with matrix arithmetic
        for world in range(self.n_world):
            for pix in range(self.n_pixel):
                matrix = pc != 0

        # We now need to check specifically for celestial coordinates since
        # these can assume correlations because of spherical distortions.
        for world1 in range(self.n_world):
            if axes[world1]['coordinate_type'] == 'celestial':
                for world2 in range(self.n_world):
                    if world1 != world2:
                        if axes[world2]['coordinate_type'] == 'celestial':
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

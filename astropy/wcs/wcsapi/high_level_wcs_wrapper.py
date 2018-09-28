from .high_level_api import HighLevelWCSMixin

__all__ = ['HighLevelWCSWrapper']


class HighLevelWCSWrapper(HighLevelWCSMixin):
    """
    Wrapper class that can take any `BaseLowLevelWCS` object and expose the
    high-level WCS API.
    """

    def __init__(self, low_level_wcs):
        self._low_level_wcs = low_level_wcs

    @property
    def low_level_wcs(self):
        return self._low_level_wcs

    @property
    def pixel_n_dim(self):
        """
        See `~BaseLowLevelWCS.world_n_dim`
        """
        return self.low_level_wcs.pixel_n_dim

    @property
    def world_n_dim(self):
        """
        See `~BaseLowLevelWCS.world_n_dim`
        """
        return self.low_level_wcs.world_n_dim

    @property
    def world_axis_physical_types(self):
        """
        See `~BaseLowLevelWCS.world_axis_physical_types`
        """
        return self.low_level_wcs.world_axis_physical_types

    @property
    def world_axis_units(self):
        """
        See `~BaseLowLevelWCS.world_axis_units`
        """
        return self.low_level_wcs.world_axis_units

    @property
    def array_shape(self):
        """
        See `~BaseLowLevelWCS.array_shape`
        """
        return self.low_level_wcs.array_shape

    @property
    def pixel_bounds(self):
        """
        See `~BaseLowLevelWCS.pixel_bounds`
        """
        return self.low_level_wcs.pixel_bounds

    @property
    def axis_correlation_matrix(self):
        """
        See `~BaseLowLevelWCS.axis_correlation_matrix`
        """
        return self.low_level_wcs.pixel_bounds

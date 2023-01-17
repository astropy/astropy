from .high_level_api import HighLevelWCSMixin
from .low_level_api import BaseLowLevelWCS
from .utils import wcs_info_str

__all__ = ["HighLevelWCSWrapper"]


class HighLevelWCSWrapper(HighLevelWCSMixin):
    """
    Wrapper class that can take any :class:`~astropy.wcs.wcsapi.BaseLowLevelWCS`
    object and expose the high-level WCS API.
    """

    def __init__(self, low_level_wcs):
        if not isinstance(low_level_wcs, BaseLowLevelWCS):
            raise TypeError(
                "Input to a HighLevelWCSWrapper must be a low level WCS object"
            )

        self._low_level_wcs = low_level_wcs

    @property
    def low_level_wcs(self):
        return self._low_level_wcs

    @property
    def pixel_n_dim(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim`.
        """
        return self.low_level_wcs.pixel_n_dim

    @property
    def world_n_dim(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim`.
        """
        return self.low_level_wcs.world_n_dim

    @property
    def world_axis_physical_types(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_physical_types`.
        """
        return self.low_level_wcs.world_axis_physical_types

    @property
    def world_axis_units(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_units`.
        """
        return self.low_level_wcs.world_axis_units

    @property
    def array_shape(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.array_shape`.
        """
        return self.low_level_wcs.array_shape

    @property
    def pixel_bounds(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_bounds`.
        """
        return self.low_level_wcs.pixel_bounds

    @property
    def axis_correlation_matrix(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.axis_correlation_matrix`.
        """
        return self.low_level_wcs.axis_correlation_matrix

    def _as_mpl_axes(self):
        """
        See `~astropy.wcs.wcsapi.BaseLowLevelWCS._as_mpl_axes`.
        """
        return self.low_level_wcs._as_mpl_axes()

    def __str__(self):
        return wcs_info_str(self.low_level_wcs)

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

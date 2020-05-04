import abc

from astropy.wcs.wcsapi import BaseLowLevelWCS, wcs_info_str


class BaseWCSWrapper(BaseLowLevelWCS, metaclass=abc.ABCMeta):
    """
    A base wrapper class for things that modify Low Level WCSes.

    This wrapper implements a transparent wrapper to many of the properties,
    with the idea that not all of them would need to be overridden in your
    wrapper, but some probably will.

    Parameters
    ----------
    wcs : `astropy.wcs.wcsapi.BaseLowLevelWCS`
        The WCS object to wrap
    """
    def __init__(self, wcs, *args, **kwargs):
        self._wcs = wcs

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
    def world_axis_units(self):
        return self._wcs.world_axis_units

    @property
    def world_axis_object_components(self):
        return self._wcs.world_axis_object_components

    @property
    def world_axis_object_classes(self):
        return self._wcs.world_axis_object_classes

    @property
    def pixel_shape(self):
        return self._wcs.pixel_shape

    @property
    def pixel_bounds(self):
        return self._wcs.pixel_bounds

    @property
    def pixel_axis_names(self):
        return self._wcs.pixel_axis_names

    @property
    def world_axis_names(self):
        return self._wcs.world_axis_names

    @property
    def axis_correlation_matrix(self):
        return self._wcs.axis_correlation_matrix

    @property
    def serialized_classes(self):
        return self._wcs.serialized_classes

    @abc.abstractmethod
    def pixel_to_world_values(self, *pixel_arrays):
        pass

    @abc.abstractmethod
    def world_to_pixel_values(self, *world_arrays):
        pass

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def __str__(self):
        return wcs_info_str(self)

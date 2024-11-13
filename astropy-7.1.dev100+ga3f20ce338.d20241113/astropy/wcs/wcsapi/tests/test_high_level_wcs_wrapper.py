import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.coordinates import SkyCoord
from astropy.wcs.wcsapi.high_level_wcs_wrapper import HighLevelWCSWrapper
from astropy.wcs.wcsapi.low_level_api import BaseLowLevelWCS


class CustomLowLevelWCS(BaseLowLevelWCS):
    @property
    def pixel_n_dim(self):
        return 2

    @property
    def world_n_dim(self):
        return 2

    @property
    def world_axis_physical_types(self):
        return ["pos.eq.ra", "pos.eq.dec"]

    @property
    def world_axis_units(self):
        return ["deg", "deg"]

    def pixel_to_world_values(self, *pixel_arrays):
        return [np.asarray(pix) * 2 for pix in pixel_arrays]

    def world_to_pixel_values(self, *world_arrays):
        return [np.asarray(world) / 2 for world in world_arrays]

    @property
    def world_axis_object_components(self):
        return [
            ("test", 0, "spherical.lon.degree"),
            ("test", 1, "spherical.lat.degree"),
        ]

    @property
    def world_axis_object_classes(self):
        return {"test": (SkyCoord, (), {"unit": "deg"})}


def test_wrapper():
    wcs = CustomLowLevelWCS()

    wrapper = HighLevelWCSWrapper(wcs)

    coord = wrapper.pixel_to_world(1, 2)

    assert isinstance(coord, SkyCoord)
    assert coord.isscalar

    x, y = wrapper.world_to_pixel(coord)

    assert_allclose(x, 1)
    assert_allclose(y, 2)

    assert wrapper.low_level_wcs is wcs
    assert wrapper.pixel_n_dim == 2
    assert wrapper.world_n_dim == 2
    assert wrapper.world_axis_physical_types == ["pos.eq.ra", "pos.eq.dec"]
    assert wrapper.world_axis_units == ["deg", "deg"]
    assert wrapper.array_shape is None
    assert wrapper.pixel_bounds is None
    assert np.all(wrapper.axis_correlation_matrix)


def test_wrapper_invalid():
    class InvalidCustomLowLevelWCS(CustomLowLevelWCS):
        @property
        def world_axis_object_classes(self):
            return {}

    wcs = InvalidCustomLowLevelWCS()

    wrapper = HighLevelWCSWrapper(wcs)

    with pytest.raises(KeyError):
        wrapper.pixel_to_world(1, 2)

import re

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy.wcs import WCS
from astropy.wcs.wcsapi.high_level_api import (
    HighLevelWCSMixin,
    high_level_objects_to_values,
    values_to_high_level_objects,
)
from astropy.wcs.wcsapi.low_level_api import BaseLowLevelWCS


class DoubleLowLevelWCS(BaseLowLevelWCS):
    """
    Basic dummy transformation that doubles values.
    """

    def pixel_to_world_values(self, *pixel_arrays):
        return [np.asarray(pix) * 2 for pix in pixel_arrays]

    def world_to_pixel_values(self, *world_arrays):
        return [np.asarray(world) / 2 for world in world_arrays]


class SimpleDuplicateWCS(DoubleLowLevelWCS, HighLevelWCSMixin):
    """
    This example WCS has two of the world coordinates that use the same class,
    which triggers a different path in the high level WCS code.
    """

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

    @property
    def world_axis_object_components(self):
        return [("test1", 0, "value"), ("test2", 0, "value")]

    @property
    def world_axis_object_classes(self):
        return {
            "test1": (Quantity, (), {"unit": "deg"}),
            "test2": (Quantity, (), {"unit": "deg"}),
        }


def test_simple_duplicate():
    # Make sure that things work properly when the low-level WCS uses the same
    # class for two of the coordinates.

    wcs = SimpleDuplicateWCS()
    q1, q2 = wcs.pixel_to_world(1, 2)

    assert isinstance(q1, Quantity)
    assert isinstance(q2, Quantity)

    x, y = wcs.world_to_pixel(q1, q2)

    assert_allclose(x, 1)
    assert_allclose(y, 2)


class SkyCoordDuplicateWCS(DoubleLowLevelWCS, HighLevelWCSMixin):
    """
    This example WCS returns two SkyCoord objects which, which triggers a
    different path in the high level WCS code.
    """

    @property
    def pixel_n_dim(self):
        return 4

    @property
    def world_n_dim(self):
        return 4

    @property
    def world_axis_physical_types(self):
        return ["pos.eq.ra", "pos.eq.dec", "pos.galactic.lon", "pos.galactic.lat"]

    @property
    def world_axis_units(self):
        return ["deg", "deg", "deg", "deg"]

    @property
    def world_axis_object_components(self):
        # Deliberately use 'ra'/'dec' here to make sure that string argument
        # names work properly.
        return [
            ("test1", "ra", "spherical.lon.degree"),
            ("test1", "dec", "spherical.lat.degree"),
            ("test2", 0, "spherical.lon.degree"),
            ("test2", 1, "spherical.lat.degree"),
        ]

    @property
    def world_axis_object_classes(self):
        return {
            "test1": (SkyCoord, (), {"unit": "deg"}),
            "test2": (SkyCoord, (), {"unit": "deg", "frame": "galactic"}),
        }


def test_skycoord_duplicate():
    # Make sure that things work properly when the low-level WCS uses the same
    # class, and specifically a SkyCoord for two of the coordinates.

    wcs = SkyCoordDuplicateWCS()
    c1, c2 = wcs.pixel_to_world(1, 2, 3, 4)

    assert isinstance(c1, SkyCoord)
    assert isinstance(c2, SkyCoord)

    x, y, z, a = wcs.world_to_pixel(c1, c2)

    assert_allclose(x, 1)
    assert_allclose(y, 2)
    assert_allclose(z, 3)
    assert_allclose(a, 4)


class SerializedWCS(DoubleLowLevelWCS, HighLevelWCSMixin):
    """
    WCS with serialized classes
    """

    @property
    def serialized_classes(self):
        return True

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

    @property
    def world_axis_object_components(self):
        return [("test", 0, "value")]

    @property
    def world_axis_object_classes(self):
        return {
            "test": (
                "astropy.units.Quantity",
                (),
                {"unit": ("astropy.units.Unit", ("deg",), {})},
            )
        }


def test_serialized_classes():
    wcs = SerializedWCS()
    q = wcs.pixel_to_world(1)

    assert isinstance(q, Quantity)

    x = wcs.world_to_pixel(q)

    assert_allclose(x, 1)


def test_objects_to_values():
    wcs = SkyCoordDuplicateWCS()
    c1, c2 = wcs.pixel_to_world(1, 2, 3, 4)

    values = high_level_objects_to_values(c1, c2, low_level_wcs=wcs)

    assert np.allclose(values, [2, 4, 6, 8])


def test_values_to_objects():
    wcs = SkyCoordDuplicateWCS()
    c1, c2 = wcs.pixel_to_world(1, 2, 3, 4)

    c1_out, c2_out = values_to_high_level_objects(*[2, 4, 6, 8], low_level_wcs=wcs)
    assert c1.ra == c1_out.ra
    assert c2.l == c2_out.l

    assert c1.dec == c1_out.dec
    assert c2.b == c2_out.b


class InvalidWCSQuantity(SkyCoordDuplicateWCS):
    """
    WCS which defines ``world_axis_object_components`` which returns Quantity
    instead of bare Numpy arrays, which can cause issues. This is for a
    regression test to make sure that we don't return Quantities from
    ``world_axis_object_components``.
    """

    @property
    def world_axis_object_components(self):
        return [
            ("test1", "ra", "spherical.lon"),
            ("test1", "dec", "spherical.lat"),
            ("test2", 0, "spherical.lon"),
            ("test2", 1, "spherical.lat"),
        ]


def test_objects_to_values_invalid_type():
    wcs = InvalidWCSQuantity()
    c1, c2 = wcs.pixel_to_world(1, 2, 3, 4)
    with pytest.raises(
        TypeError,
        match=(
            re.escape(
                "WCS world_axis_object_components results in values which are not "
                "scalars or plain Numpy arrays (got <class "
                "'astropy.coordinates.angles.core.Longitude'>)"
            )
        ),
    ):
        high_level_objects_to_values(c1, c2, low_level_wcs=wcs)


def test_values_to_objects_invalid_type():
    wcs = SkyCoordDuplicateWCS()
    c1, c2 = wcs.pixel_to_world(1, 2, 3, 4)
    with pytest.raises(
        TypeError,
        match=(
            re.escape(
                "Expected world coordinates as scalars or plain Numpy arrays (got "
                "<class 'astropy.units.quantity.Quantity'>)"
            )
        ),
    ):
        values_to_high_level_objects(2 * u.m, 4, 6, 8, low_level_wcs=wcs)


class MinimalHighLevelWCS(HighLevelWCSMixin):
    def __init__(self, low_level_wcs):
        self._low_level_wcs = low_level_wcs

    @property
    def low_level_wcs(self):
        return self._low_level_wcs


def test_minimal_mixin_subclass():
    # Regression test for a bug that caused coordinate conversions to fail
    # unless the WCS dimensions were defined on the high level WCS (which they
    # are not required to be)

    fits_wcs = WCS(naxis=2)
    high_level_wcs = MinimalHighLevelWCS(fits_wcs)

    coord = high_level_wcs.pixel_to_world(1, 2)
    pixel = high_level_wcs.world_to_pixel(*coord)

    coord = high_level_wcs.array_index_to_world(1, 2)
    pixel = high_level_wcs.world_to_array_index(*coord)

    assert_allclose(pixel, (1, 2))


def test_world_to_array_index_nan():
    # see https://github.com/astropy/astropy/issues/17227
    wcs1 = WCS(naxis=1)
    wcs1.wcs.crpix = (1,)
    wcs1.wcs.set()
    wcs1.pixel_bounds = [None]

    res1 = wcs1.world_to_array_index(*wcs1.pixel_to_world((5,)))
    assert not np.any(np.isnan(res1))
    assert res1.ndim == 0
    assert res1.item() == 5

    wcs2 = WCS(naxis=2)
    wcs2.wcs.crpix = (1, 1)
    wcs2.wcs.set()
    wcs2.pixel_bounds = [None, (-0.5, 3.5)]

    res2 = wcs2.world_to_array_index(*wcs2.pixel_to_world(5, 5))
    assert not np.any(np.isnan(res2))
    assert isinstance(res2, tuple)
    assert len(res2) == 2
    assert res2 == (np.iinfo(int).min, 5)

import numpy as np
from numpy.testing import assert_allclose

from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from astropy.wcs.wcsapi.low_level_api import BaseLowLevelWCS
from astropy.wcs.wcsapi.high_level_api import HighLevelWCSMixin


class DoubleLowLevelWCS(BaseLowLevelWCS):
    """
    Basic dummy transformation that doubles values.
    """

    def pixel_to_world_values(self, *pixel_arrays):
        return [np.asarray(pix) * 2 for pix in pixel_arrays]

    def array_index_to_world_values(self, *index_arrays):
        return [np.asarray(pix) * 2 for pix in index_arrays]

    def world_to_pixel_values(self, *world_arrays):
        return [np.asarray(world) / 2 for world in world_arrays]

    def world_to_array_index_values(self, *world_arrays):
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
        return ['pos.eq.ra', 'pos.eq.dec']

    @property
    def world_axis_units(self):
        return ['deg', 'deg']

    @property
    def world_axis_object_components(self):
        return [('test1', 0, 'value'),
                ('test2', 0, 'value')]

    @property
    def world_axis_object_classes(self):
        return {'test1': (Quantity, (), {'unit': 'deg'}),
                'test2': (Quantity, (), {'unit': 'deg'})}


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
        return ['pos.eq.ra', 'pos.eq.dec', 'pos.galactic.lon', 'pos.galactic.lat']

    @property
    def world_axis_units(self):
        return ['deg', 'deg', 'deg', 'deg']

    @property
    def world_axis_object_components(self):
        # Deliberately use 'ra'/'dec' here to make sure that string argument
        # names work properly.
        return [('test1', 'ra', 'spherical.lon.degree'),
                ('test1', 'dec', 'spherical.lat.degree'),
                ('test2', 0, 'spherical.lon.degree'),
                ('test2', 1, 'spherical.lat.degree')]

    @property
    def world_axis_object_classes(self):
        return {'test1': (SkyCoord, (), {'unit': 'deg'}),
                'test2': (SkyCoord, (), {'unit': 'deg', 'frame': 'galactic'})}


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
        return ['pos.eq.ra', 'pos.eq.dec']

    @property
    def world_axis_units(self):
        return ['deg', 'deg']

    @property
    def world_axis_object_components(self):
        return [('test', 0, 'value')]

    @property
    def world_axis_object_classes(self):
        return {'test': ('astropy.units.Quantity', (),
                         {'unit': ('astropy.units.Unit', ('deg',), {})})}


def test_serialized_classes():

    wcs = SerializedWCS()
    q = wcs.pixel_to_world(1)

    assert isinstance(q, Quantity)

    x = wcs.world_to_pixel(q)

    assert_allclose(x, 1)

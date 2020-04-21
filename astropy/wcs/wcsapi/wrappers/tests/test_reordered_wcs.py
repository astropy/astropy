import re

import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy.units import Quantity
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.wcsapi import ReorderedLowLevelWCS, HighLevelWCSWrapper
from astropy.tests.helper import assert_quantity_allclose


@pytest.fixture
def spectral_wcs(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def celestial_wcs(request):
    return request.getfixturevalue(request.param)


EXPECTED_SPECTRAL_CUBE_REPR = """
ReorderedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (6, 3, 7)

Pixel Dim  Axis Name  Data size  Bounds
        0  None               7  (1, 7)
        1  None               3  (1, 2.5)
        2  None               6  (-1, 5)

World Dim  Axis Name        Physical Type  Units
        0  Frequency        em.freq        Hz
        1  Right Ascension  pos.eq.ra      deg
        2  Declination      pos.eq.dec     deg

Correlation between pixel and world axes:

             Pixel Dim
World Dim    0    1    2
        0   no  yes   no
        1  yes   no  yes
        2  yes   no  yes
""".strip()


def test_spectral_cube(spectral_cube_3d_fitswcs):

    wcs = ReorderedLowLevelWCS(spectral_cube_3d_fitswcs,
                               pixel_order=[1, 2, 0],
                               world_order=[2, 0, 1])

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert tuple(wcs.world_axis_physical_types) == ('em.freq', 'pos.eq.ra', 'pos.eq.dec')
    assert tuple(wcs.world_axis_units) == ('Hz', 'deg', 'deg')
    assert tuple(wcs.pixel_axis_names) == ('', '', '')
    assert tuple(wcs.world_axis_names) == ('Frequency',
                                           'Right Ascension',
                                           'Declination')
    assert_equal(wcs.axis_correlation_matrix, np.array([[0, 1, 0],
                                                        [1, 0, 1],
                                                        [1, 0, 1]]))

    assert wcs.pixel_shape == (7, 3, 6)
    assert wcs.array_shape == (6, 3, 7)
    assert wcs.pixel_bounds == ((1, 7), (1, 2.5), (-1, 5))

    pixel_scalar = (1.3, 2.3, 4.3)
    world_scalar = (-1.91e10, 5.4, -9.4)
    assert_allclose(wcs.pixel_to_world_values(*pixel_scalar), world_scalar)
    assert_allclose(wcs.array_index_to_world_values(*pixel_scalar[::-1]), world_scalar)
    assert_allclose(wcs.world_to_pixel_values(*world_scalar), pixel_scalar)
    assert_allclose(wcs.world_to_array_index_values(*world_scalar), [4, 2, 1])

    pixel_array = (np.array([1.3, 1.4]),
                   np.array([2.3, 2.4]),
                   np.array([4.3, 4.4]))
    world_array = (np.array([-1.91e10, -1.88e10]),
                   np.array([5.4, 5.2]),
                   np.array([-9.4, -9.2]))
    assert_allclose(wcs.pixel_to_world_values(*pixel_array), world_array)
    assert_allclose(wcs.array_index_to_world_values(*pixel_array[::-1]), world_array)
    assert_allclose(wcs.world_to_pixel_values(*world_array), pixel_array)
    assert_allclose(wcs.world_to_array_index_values(*world_array),
                    [[4, 4], [2, 2], [1, 1]])

    wcs_hl = HighLevelWCSWrapper(wcs)

    spectral, celestial = wcs_hl.pixel_to_world(*pixel_scalar)
    assert isinstance(spectral, Quantity)
    assert_quantity_allclose(spectral, world_scalar[0] * u.Hz)
    assert isinstance(celestial, SkyCoord)
    assert_quantity_allclose(celestial.ra, world_scalar[1] * u.deg)
    assert_quantity_allclose(celestial.dec, world_scalar[2] * u.deg)

    spectral, celestial = wcs_hl.pixel_to_world(*pixel_array)
    assert isinstance(spectral, Quantity)
    assert_quantity_allclose(spectral, world_array[0] * u.Hz)
    assert isinstance(celestial, SkyCoord)
    assert_quantity_allclose(celestial.ra, world_array[1] * u.deg)
    assert_quantity_allclose(celestial.dec, world_array[2] * u.deg)

    assert str(wcs) == EXPECTED_SPECTRAL_CUBE_REPR
    assert EXPECTED_SPECTRAL_CUBE_REPR in repr(wcs)


@pytest.mark.parametrize('order', [(1,), (1, 2, 2), (0, 1, 2, 3)])
def test_invalid(spectral_cube_3d_fitswcs, order):

    with pytest.raises(ValueError, match=re.escape('pixel_order should be a permutation of [0, 1, 2]')):
        ReorderedLowLevelWCS(spectral_cube_3d_fitswcs,
                             pixel_order=order,
                             world_order=[2, 0, 1])

    with pytest.raises(ValueError, match=re.escape('world_order should be a permutation of [0, 1, 2]')):
        ReorderedLowLevelWCS(spectral_cube_3d_fitswcs,
                             pixel_order=[1, 2, 0],
                             world_order=order)

from itertools import product

import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy.units import Quantity
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.wcsapi import CompoundLowLevelWCS, HighLevelWCSWrapper
from astropy.tests.helper import assert_quantity_allclose


@pytest.fixture
def spectral_wcs(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def celestial_wcs(request):
    return request.getfixturevalue(request.param)


EXPECTED_CELESTIAL_SPECTRAL_APE14_REPR = """
CompoundLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (7, 6, 3)

Pixel Dim  Axis Name  Data size  Bounds
        0  None               3  (1, 2)
        1  None               6  (-1, 5)
        2  None               7  (1, 7)

World Dim  Axis Name        Physical Type  Units
        0  Frequency        em.freq        Hz
        1  Right Ascension  pos.eq.ra      deg
        2  Declination      pos.eq.dec     deg

Correlation between pixel and world axes:

             Pixel Dim
World Dim    0    1    2
        0  yes   no   no
        1   no  yes  yes
        2   no  yes  yes
""".strip()


@pytest.mark.parametrize(('spectral_wcs', 'celestial_wcs'),
                         product(['spectral_1d_ape14_wcs', 'spectral_1d_fitswcs'],
                                 ['celestial_2d_ape14_wcs', 'celestial_2d_fitswcs']),
                         indirect=True)
def test_celestial_spectral_ape14(spectral_wcs, celestial_wcs):

    wcs = CompoundLowLevelWCS(spectral_wcs, celestial_wcs)

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert tuple(wcs.world_axis_physical_types) == ('em.freq', 'pos.eq.ra', 'pos.eq.dec')
    assert tuple(wcs.world_axis_units) == ('Hz', 'deg', 'deg')
    assert tuple(wcs.pixel_axis_names) == ('', '', '')
    assert tuple(wcs.world_axis_names) == ('Frequency',
                                           'Right Ascension',
                                           'Declination')
    assert_equal(wcs.axis_correlation_matrix, np.array([[1, 0, 0],
                                                        [0, 1, 1],
                                                        [0, 1, 1]]))

    # If any of the individual shapes are None, return None overall
    assert wcs.pixel_shape is None
    assert wcs.array_shape is None
    assert wcs.pixel_bounds is None

    # Set the shape and bounds on the spectrum and test again
    spectral_wcs.pixel_shape = (3,)
    spectral_wcs.pixel_bounds = [(1, 2)]
    assert wcs.pixel_shape == (3, 6, 7)
    assert wcs.array_shape == (7, 6, 3)
    assert wcs.pixel_bounds == ((1, 2), (-1, 5), (1, 7))

    pixel_scalar = (2.3, 4.3, 1.3)
    world_scalar = (-1.91e10, 5.4, -9.4)
    assert_allclose(wcs.pixel_to_world_values(*pixel_scalar), world_scalar)
    assert_allclose(wcs.array_index_to_world_values(*pixel_scalar[::-1]), world_scalar)
    assert_allclose(wcs.world_to_pixel_values(*world_scalar), pixel_scalar)
    assert_allclose(wcs.world_to_array_index_values(*world_scalar), [1, 4, 2])

    pixel_array = (np.array([2.3, 2.4]),
                   np.array([4.3, 4.4]),
                   np.array([1.3, 1.4]))
    world_array = (np.array([-1.91e10, -1.88e10]),
                   np.array([5.4, 5.2]),
                   np.array([-9.4, -9.2]))
    assert_allclose(wcs.pixel_to_world_values(*pixel_array), world_array)
    assert_allclose(wcs.array_index_to_world_values(*pixel_array[::-1]), world_array)
    assert_allclose(wcs.world_to_pixel_values(*world_array), pixel_array)
    assert_allclose(wcs.world_to_array_index_values(*world_array),
                    [[1, 1], [4, 4], [2, 2]])

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

    assert str(wcs) == EXPECTED_CELESTIAL_SPECTRAL_APE14_REPR
    assert EXPECTED_CELESTIAL_SPECTRAL_APE14_REPR in repr(wcs)

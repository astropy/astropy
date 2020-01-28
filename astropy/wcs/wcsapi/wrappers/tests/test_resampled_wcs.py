import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.wcsapi import ResampledLowLevelWCS, HighLevelWCSWrapper
from astropy.tests.helper import assert_quantity_allclose


@pytest.fixture
def celestial_wcs(request):
    return request.getfixturevalue(request.param)


EXPECTED_2D_REPR = """
ResampledLowLevelWCS Transformation

This transformation has 2 pixel and 2 world dimensions

Array shape (Numpy order): (2.3333333333333335, 15.0)

Pixel Dim  Axis Name  Data size  Bounds
        0  None              15  (-2.5, 12.5)
        1  None         2.33333  (0.3333333333333333, 2.3333333333333335)

World Dim  Axis Name        Physical Type  Units
        0  Right Ascension  pos.eq.ra      deg
        1  Declination      pos.eq.dec     deg

Correlation between pixel and world axes:

           Pixel Dim
World Dim    0    1
        0  yes  yes
        1  yes  yes
""".strip()


@pytest.mark.parametrize('celestial_wcs',
                         ['celestial_2d_ape14_wcs', 'celestial_2d_fitswcs'],
                         indirect=True)
def test_2d(celestial_wcs):

    # Upsample along the first pixel dimension and downsample along the second
    # pixel dimension.
    wcs = ResampledLowLevelWCS(celestial_wcs, [0.4, 3])

    # The following shouldn't change compared to the original WCS
    assert wcs.pixel_n_dim == 2
    assert wcs.world_n_dim == 2
    assert tuple(wcs.world_axis_physical_types) == ('pos.eq.ra', 'pos.eq.dec')
    assert tuple(wcs.world_axis_units) == ('deg', 'deg')
    assert tuple(wcs.pixel_axis_names) == ('', '')
    assert tuple(wcs.world_axis_names) == ('Right Ascension',
                                           'Declination')
    assert_equal(wcs.axis_correlation_matrix, np.ones((2,2)))

    # Shapes and bounds should be floating-point if needed
    assert_allclose(wcs.pixel_shape, (15, 7/3))
    assert_allclose(wcs.array_shape, (7/3, 15))
    assert_allclose(wcs.pixel_bounds, ((-2.5, 12.5), (1/3, 7/3)))

    pixel_scalar = (2.3, 4.3)
    world_scalar = (12.16, 13.8)
    assert_allclose(wcs.pixel_to_world_values(*pixel_scalar), world_scalar)
    assert_allclose(wcs.array_index_to_world_values(*pixel_scalar[::-1]), world_scalar)
    assert_allclose(wcs.world_to_pixel_values(*world_scalar), pixel_scalar)
    assert_allclose(wcs.world_to_array_index_values(*world_scalar), [4, 2])

    pixel_array = (np.array([2.3, 2.4]),
                   np.array([4.3, 4.4]))
    world_array = (np.array([12.16, 12.08]),
                   np.array([13.8, 14.4]))
    assert_allclose(wcs.pixel_to_world_values(*pixel_array), world_array)
    assert_allclose(wcs.array_index_to_world_values(*pixel_array[::-1]), world_array)
    assert_allclose(wcs.world_to_pixel_values(*world_array), pixel_array)
    assert_allclose(wcs.world_to_array_index_values(*world_array),
                    [[4, 4], [2, 2]])

    wcs_hl = HighLevelWCSWrapper(wcs)

    celestial = wcs_hl.pixel_to_world(*pixel_scalar)
    assert isinstance(celestial, SkyCoord)
    assert_quantity_allclose(celestial.ra, world_scalar[0] * u.deg)
    assert_quantity_allclose(celestial.dec, world_scalar[1] * u.deg)

    celestial = wcs_hl.pixel_to_world(*pixel_array)
    assert isinstance(celestial, SkyCoord)
    assert_quantity_allclose(celestial.ra, world_array[0] * u.deg)
    assert_quantity_allclose(celestial.dec, world_array[1] * u.deg)

    assert str(wcs) == EXPECTED_2D_REPR
    assert repr(wcs) == EXPECTED_2D_REPR


@pytest.mark.parametrize('celestial_wcs',
                         ['celestial_2d_ape14_wcs', 'celestial_2d_fitswcs'],
                         indirect=True)
def test_scalar_factor(celestial_wcs):

    wcs = ResampledLowLevelWCS(celestial_wcs, 2)

    pixel_scalar = (2.3, 4.3)
    world_scalar = (4.8, 5.2)
    assert_allclose(wcs.pixel_to_world_values(*pixel_scalar), world_scalar)
    assert_allclose(wcs.array_index_to_world_values(*pixel_scalar[::-1]), world_scalar)
    assert_allclose(wcs.world_to_pixel_values(*world_scalar), pixel_scalar)
    assert_allclose(wcs.world_to_array_index_values(*world_scalar), [4, 2])

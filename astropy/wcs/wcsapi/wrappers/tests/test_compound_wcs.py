import pytest

import numpy as np
from numpy.testing import assert_allclose

from astropy.units import Quantity
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.wcsapi import CompoundLowLevelWCS, HighLevelWCSWrapper
from astropy.tests.helper import assert_quantity_allclose


def test_celestial_spectral_ape14(spectral_1d_ape14_wcs,
                                  celestial_2d_ape14_wcs):

    wcs = CompoundLowLevelWCS(spectral_1d_ape14_wcs, celestial_2d_ape14_wcs)

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert tuple(wcs.world_axis_physical_types) == ('em.freq', 'pos.eq.ra', 'pos.eq.dec')
    assert tuple(wcs.world_axis_units) == ('GHz', 'deg', 'deg')

    pix_scalar = (2.3, 4.3, 1.3)
    wrl_scalar = (6.9, 8.6, 2.6)
    assert_allclose(wcs.pixel_to_world_values(*pix_scalar), wrl_scalar)
    assert_allclose(wcs.array_index_to_world_values(*pix_scalar[::-1]), wrl_scalar)
    assert_allclose(wcs.world_to_pixel_values(*wrl_scalar), pix_scalar)
    assert_allclose(wcs.world_to_array_index_values(*wrl_scalar), [1, 4, 2])

    pix_array = (np.array([2.3, 2.4]),
                 np.array([4.3, 4.4]),
                 np.array([1.3, 1.4]))
    wrl_array = (np.array([6.9, 7.2]),
                 np.array([8.6, 8.8]),
                 np.array([2.6, 2.8]))
    assert_allclose(wcs.pixel_to_world_values(*pix_array), wrl_array)
    assert_allclose(wcs.array_index_to_world_values(*pix_array[::-1]), wrl_array)
    assert_allclose(wcs.world_to_pixel_values(*wrl_array), pix_array)
    assert_allclose(wcs.world_to_array_index_values(*wrl_array),
                    [[1, 1], [4, 4], [2, 2]])

    wcs_hl = HighLevelWCSWrapper(wcs)
    spectral, celestial = wcs_hl.pixel_to_world(*pix_scalar)
    assert isinstance(spectral, Quantity)
    assert_quantity_allclose(spectral, 6.9 * u.GHz)
    assert isinstance(celestial, SkyCoord)
    assert_quantity_allclose(celestial.ra, 8.6 * u.deg)
    assert_quantity_allclose(celestial.dec, 2.6 * u.deg)

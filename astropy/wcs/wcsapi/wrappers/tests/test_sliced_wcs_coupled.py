"""Test for coupled dimensions in SlicedLowLevelWCS world_to_pixel."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

from astropy.wcs import WCS
from astropy.wcs.wcsapi.wrappers.sliced_wcs import SlicedLowLevelWCS
from astropy.wcs.wcsapi.high_level_wcs_wrapper import HighLevelWCSWrapper
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import wcs_to_celestial_frame
import astropy.units as u


def test_sliced_wcs_world_to_pixel_coupled_dimensions():
    """
    Test world_to_pixel on SlicedLowLevelWCS with coupled spectral and spatial dimensions.

    This is a regression test for a bug where world_to_pixel would return incorrect
    results when a WCS with coupled dimensions (via PC matrix) was sliced.

    See: https://github.com/astropy/astropy/issues/XXXXX
    """
    # Create a 3D WCS with coupled spatial and spectral dimensions
    nx = 100
    ny = 25
    nz = 2
    wcs_header = {
        'WCSAXES': 3,
        'CRPIX1': (nx + 1)/2,
        'CRPIX2': (ny + 1)/2,
        'CRPIX3': 1.0,
        'PC1_1': 0.0,
        'PC1_2': -1.0,
        'PC1_3': 0.0,
        'PC2_1': 1.0,
        'PC2_2': 0.0,
        'PC2_3': -1.0,  # This couples spatial and spectral dimensions
        'CDELT1': 5,
        'CDELT2': 5,
        'CDELT3': 0.055,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'CUNIT3': 'Angstrom',
        'CTYPE1': 'HPLN-TAN',
        'CTYPE2': 'HPLT-TAN',
        'CTYPE3': 'WAVE',
        'CRVAL1': 0.0,
        'CRVAL2': 0.0,
        'CRVAL3': 1.05,
    }
    fits_wcs = WCS(header=wcs_header)

    # Test point at the reference position
    pt = SkyCoord(Tx=0*u.arcsec, Ty=0*u.arcsec,
                  frame=wcs_to_celestial_frame(fits_wcs))

    # world_to_pixel on full WCS
    px_full, py_full, pz_full = fits_wcs.world_to_pixel(pt, 1.05*u.angstrom)

    # Create sliced WCS at first wavelength slice
    ll_sliced_wcs = SlicedLowLevelWCS(fits_wcs, 0)
    hl_sliced_wcs = HighLevelWCSWrapper(ll_sliced_wcs)

    # world_to_pixel on sliced WCS should give same spatial coordinates
    px_sliced, py_sliced = hl_sliced_wcs.world_to_pixel(pt)

    # The spatial pixel coordinates should match
    assert_allclose(px_sliced, px_full, rtol=1e-5, atol=1e-10,
                    err_msg="X pixel coordinate does not match after slicing")
    assert_allclose(py_sliced, py_full, rtol=1e-5, atol=1e-10,
                    err_msg="Y pixel coordinate does not match after slicing")

    # Verify pixel_to_world gives consistent results
    world_sliced = hl_sliced_wcs.pixel_to_world(px_sliced, py_sliced)
    world_full = fits_wcs.pixel_to_world(px_full, py_full, pz_full)

    assert_allclose(world_sliced.Tx.value, world_full[0].Tx.value, rtol=1e-5, atol=1e-10)
    assert_allclose(world_sliced.Ty.value, world_full[0].Ty.value, rtol=1e-5, atol=1e-10)


def test_sliced_wcs_world_to_pixel_coupled_multiple_points():
    """Test with multiple points to ensure broadcasting works correctly."""
    nx = 100
    ny = 25
    nz = 2
    wcs_header = {
        'WCSAXES': 3,
        'CRPIX1': (nx + 1)/2,
        'CRPIX2': (ny + 1)/2,
        'CRPIX3': 1.0,
        'PC1_1': 0.0,
        'PC1_2': -1.0,
        'PC1_3': 0.0,
        'PC2_1': 1.0,
        'PC2_2': 0.0,
        'PC2_3': -1.0,
        'CDELT1': 5,
        'CDELT2': 5,
        'CDELT3': 0.055,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'CUNIT3': 'Angstrom',
        'CTYPE1': 'HPLN-TAN',
        'CTYPE2': 'HPLT-TAN',
        'CTYPE3': 'WAVE',
        'CRVAL1': 0.0,
        'CRVAL2': 0.0,
        'CRVAL3': 1.05,
    }
    fits_wcs = WCS(header=wcs_header)

    # Test with multiple points
    points = SkyCoord(Tx=[0, 10, -10]*u.arcsec, Ty=[0, 5, -5]*u.arcsec,
                      frame=wcs_to_celestial_frame(fits_wcs))

    # world_to_pixel on full WCS
    px_full, py_full, pz_full = fits_wcs.world_to_pixel(points, 1.05*u.angstrom)

    # Create sliced WCS
    ll_sliced_wcs = SlicedLowLevelWCS(fits_wcs, 0)
    hl_sliced_wcs = HighLevelWCSWrapper(ll_sliced_wcs)

    # world_to_pixel on sliced WCS
    px_sliced, py_sliced = hl_sliced_wcs.world_to_pixel(points)

    # The spatial pixel coordinates should match for all points
    assert_allclose(px_sliced, px_full, rtol=1e-5, atol=1e-10)
    assert_allclose(py_sliced, py_full, rtol=1e-5, atol=1e-10)

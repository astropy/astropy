import os

from numpy.testing import assert_equal, assert_allclose

from ... import units as u
from ...units import Quantity
from ...coordinates import ICRS, Galactic, SkyCoord
from ...io.fits import Header
from ..wcs import WCS
from ..fitswcs_low_level_api import FITSLowLevelWCS

HEADER_SIMPLE_CELESTIAL = """
WCSAXES = 2
CTYPE1  = RA---TAN
CTYPE2  = DEC--TAN
CRVAL1  = 10
CRVAL2  = 20
CRPIX1  = 30
CRPIX2  = 40
CDELT1  = -0.1
CDELT2  =  0.1
CROTA2  = 0.
CUNIT1  = deg
CUNIT2  = deg
"""

WCS_SIMPLE_CELESTIAL = WCS(Header.fromstring(HEADER_SIMPLE_CELESTIAL, sep=os.linesep))


def test_simple_celestial():

    # 2D image with lon,lat axes ordering

    llwcs = FITSLowLevelWCS(WCS_SIMPLE_CELESTIAL)

    assert llwcs.pixel_n_dim == 2
    assert llwcs.world_n_dim == 2
    assert llwcs.world_axis_physical_types == ['pos.eq.ra', 'pos.eq.dec']
    assert llwcs.world_axis_units == ['deg', 'deg']

    assert_equal(llwcs.axis_correlation_matrix, True)

    assert llwcs.world_axis_object_components == [('celestial', 0, 'spherical.lon.degree'),
                                                  ('celestial', 1, 'spherical.lat.degree')]

    assert llwcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert llwcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(llwcs.world_axis_object_classes['celestial'][2]['frame'], ICRS)
    assert llwcs.world_axis_object_classes['celestial'][2]['unit'], 'deg'

    assert_allclose(llwcs.pixel_to_world_values(29, 39), (10, 20))
    assert_allclose(llwcs.numpy_index_to_world_values(39, 29), (10, 20))

    assert_allclose(llwcs.world_to_pixel_values(10, 20), (29., 39.))
    assert_equal(llwcs.world_to_numpy_index_values(10, 20), (39, 29))


HEADER_SPECTRAL_CUBE = """
WCSAXES = 3
CTYPE1  = GLAT-CAR
CTYPE2  = FREQ
CTYPE3  = GLON-CAR
CRVAL1  = 10
CRVAL2  = 20
CRVAL3  = 25
CRPIX1  = 30
CRPIX2  = 40
CRPIX3  = 45
CDELT1  = -0.1
CDELT2  =  0.5
CDELT3  =  0.1
CUNIT1  = deg
CUNIT2  = Hz
CUNIT3  = deg
"""

WCS_SPECTRAL_CUBE = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE, sep=os.linesep))


def test_spectral_cube():

    # Spectral cube with a weird axis ordering

    llwcs = FITSLowLevelWCS(WCS_SPECTRAL_CUBE)

    assert llwcs.pixel_n_dim == 3
    assert llwcs.world_n_dim == 3
    assert llwcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert llwcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(llwcs.axis_correlation_matrix, [[True, False, True], [False, True, False], [True, False, True]])

    assert llwcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                  ('freq', 0, 'value'),
                                                  ('celestial', 0, 'spherical.lon.degree')]

    assert llwcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert llwcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(llwcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert llwcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert llwcs.world_axis_object_classes['freq'][0] is Quantity
    assert llwcs.world_axis_object_classes['freq'][1] == ()
    assert llwcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(llwcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(llwcs.numpy_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(llwcs.world_to_pixel_values(10, 20, 25), (29., 39., 44.))
    assert_equal(llwcs.world_to_numpy_index_values(10, 20, 25), (44, 39, 29))


HEADER_SPECTRAL_CUBE_NONALIGNED = HEADER_SPECTRAL_CUBE.strip() + os.linesep + """
PC2_3 = -0.5
PC3_2 = +0.5
"""

WCS_SPECTRAL_CUBE_NONALIGNED = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE_NONALIGNED, sep=os.linesep))


def test_spectral_cube_nonaligned():

    # Make sure that correlation matrix gets adjusted if there are non-identity
    # CD matrix terms.

    llwcs = FITSLowLevelWCS(WCS_SPECTRAL_CUBE_NONALIGNED)

    assert_equal(llwcs.axis_correlation_matrix, [[True, True, True], [False, True, True], [True, True, True]])

import numpy as np
from numpy.testing import assert_equal, assert_allclose

import pytest
from pytest import raises

from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy import units as u

from astropy.wcs import WCS
from astropy.wcs.wcsapi.utils import deserialize_class, wcs_info_str
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS
from astropy.wcs.wcsapi.utils import (deserialize_class, pixel_to_pixel,
                                      split_matrix, pixel_to_pixel_correlation_matrix,
                                      pixel_to_world_correlation_matrix)
from astropy.utils.misc import unbroadcast


def test_construct():

    result = deserialize_class(('astropy.units.Quantity', (10,), {'unit': 'deg'}))
    assert_quantity_allclose(result, 10 * u.deg)


def test_noconstruct():

    result = deserialize_class(('astropy.units.Quantity', (), {'unit': 'deg'}), construct=False)
    assert result == (u.Quantity, (), {'unit': 'deg'})


def test_invalid():

    with raises(ValueError) as exc:
        deserialize_class(('astropy.units.Quantity', (), {'unit': 'deg'}, ()))
    assert exc.value.args[0] == 'Expected a tuple of three values'


DEFAULT_1D_STR = """
WCS Transformation

This transformation has 1 pixel and 1 world dimensions

Array shape (Numpy order): None

Pixel Dim  Axis Name  Data size  Bounds
        0  None            None  None

World Dim  Axis Name  Physical Type  Units
        0  None       None           unknown

Correlation between pixel and world axes:

           Pixel Dim
World Dim    0
        0  yes
"""


def test_wcs_info_str():

    # The tests in test_sliced_low_level_wcs.py excercise wcs_info_str
    # extensively. This test is to ensure that the function exists and the
    # API of the function works as expected.

    wcs_empty = WCS(naxis=1)

    assert wcs_info_str(wcs_empty).strip() == DEFAULT_1D_STR.strip()


def test_pixel_to_world_correlation_matrix_celestial():

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = 'RA---TAN', 'DEC--TAN'
    wcs.wcs.set()

    assert_equal(wcs.axis_correlation_matrix, [[1, 1], [1, 1]])
    matrix, classes = pixel_to_world_correlation_matrix(wcs)
    assert_equal(matrix, [[1, 1]])
    assert classes == [SkyCoord]


def test_pixel_to_world_correlation_matrix_spectral_cube_uncorrelated():

    wcs = WCS(naxis=3)
    wcs.wcs.ctype = 'RA---TAN', 'FREQ', 'DEC--TAN'
    wcs.wcs.set()

    assert_equal(wcs.axis_correlation_matrix, [[1, 0, 1], [0, 1, 0], [1, 0, 1]])
    matrix, classes = pixel_to_world_correlation_matrix(wcs)
    assert_equal(matrix, [[1, 0, 1], [0, 1, 0]])
    assert classes == [SkyCoord, Quantity]


def test_pixel_to_world_correlation_matrix_spectral_cube_correlated():

    wcs = WCS(naxis=3)
    wcs.wcs.ctype = 'RA---TAN', 'FREQ', 'DEC--TAN'
    wcs.wcs.cd = np.ones((3, 3))
    wcs.wcs.set()

    assert_equal(wcs.axis_correlation_matrix, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])
    matrix, classes = pixel_to_world_correlation_matrix(wcs)
    assert_equal(matrix, [[1, 1, 1], [1, 1, 1]])
    assert classes == [SkyCoord, Quantity]


def test_pixel_to_pixel_correlation_matrix_celestial():

    wcs_in = WCS(naxis=2)
    wcs_in.wcs.ctype = 'RA---TAN', 'DEC--TAN'
    wcs_in.wcs.set()

    wcs_out = WCS(naxis=2)
    wcs_out.wcs.ctype = 'DEC--TAN', 'RA---TAN'
    wcs_out.wcs.set()

    matrix = pixel_to_pixel_correlation_matrix(wcs_in, wcs_out)
    assert_equal(matrix, [[1, 1], [1, 1]])


def test_pixel_to_pixel_correlation_matrix_spectral_cube_uncorrelated():

    wcs_in = WCS(naxis=3)
    wcs_in.wcs.ctype = 'RA---TAN', 'DEC--TAN', 'FREQ'
    wcs_in.wcs.set()

    wcs_out = WCS(naxis=3)
    wcs_out.wcs.ctype = 'DEC--TAN', 'FREQ', 'RA---TAN'
    wcs_out.wcs.set()

    matrix = pixel_to_pixel_correlation_matrix(wcs_in, wcs_out)
    assert_equal(matrix, [[1, 1, 0], [0, 0, 1], [1, 1, 0]])


def test_pixel_to_pixel_correlation_matrix_spectral_cube_correlated():

    # NOTE: only make one of the WCSes have correlated axes to really test this

    wcs_in = WCS(naxis=3)
    wcs_in.wcs.ctype = 'RA---TAN', 'DEC--TAN', 'FREQ'
    wcs_in.wcs.set()

    wcs_out = WCS(naxis=3)
    wcs_out.wcs.ctype = 'DEC--TAN', 'FREQ', 'RA---TAN'
    wcs_out.wcs.cd = np.ones((3, 3))
    wcs_out.wcs.set()

    matrix = pixel_to_pixel_correlation_matrix(wcs_in, wcs_out)
    assert_equal(matrix, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])


def test_pixel_to_pixel_correlation_matrix_mismatch():

    wcs_in = WCS(naxis=2)
    wcs_in.wcs.ctype = 'RA---TAN', 'DEC--TAN'
    wcs_in.wcs.set()

    wcs_out = WCS(naxis=3)
    wcs_out.wcs.ctype = 'DEC--TAN', 'FREQ', 'RA---TAN'
    wcs_out.wcs.set()

    with pytest.raises(ValueError) as exc:
        pixel_to_pixel_correlation_matrix(wcs_in, wcs_out)
    assert exc.value.args[0] == "The two WCS return a different number of world coordinates"

    wcs3 = WCS(naxis=2)
    wcs3.wcs.ctype = 'FREQ', 'PIXEL'
    wcs3.wcs.set()

    with pytest.raises(ValueError) as exc:
        pixel_to_pixel_correlation_matrix(wcs_out, wcs3)
    assert exc.value.args[0] == "The world coordinate types of the two WCS do not match"

    wcs4 = WCS(naxis=4)
    wcs4.wcs.ctype = 'RA---TAN', 'DEC--TAN', 'Q1', 'Q2'
    wcs4.wcs.cunit = ['deg', 'deg', 'm/s', 'm/s']
    wcs4.wcs.set()

    wcs5 = WCS(naxis=4)
    wcs5.wcs.ctype = 'Q1', 'RA---TAN', 'DEC--TAN', 'Q2'
    wcs5.wcs.cunit = ['m/s', 'deg', 'deg', 'm/s']
    wcs5.wcs.set()

    with pytest.raises(ValueError) as exc:
        pixel_to_pixel_correlation_matrix(wcs4, wcs5)
    assert exc.value.args[0] == "World coordinate order doesn't match and automatic matching is ambiguous"


def test_pixel_to_pixel_correlation_matrix_nonsquare():

    # Here we set up an input WCS that maps 3 pixel coordinates to 4 world
    # coordinates - the idea is to make sure that things work fine in cases
    # where the number of input and output pixel coordinates do not match.

    class FakeWCS(object):
        pass

    wcs_in = FakeWCS()
    wcs_in.low_level_wcs = wcs_in
    wcs_in.pixel_n_dim = 3
    wcs_in.world_n_dim = 4
    wcs_in.axis_correlation_matrix = [[True, True, False],
                                      [True, True, False],
                                      [True, True, False],
                                      [False, False, True]]
    wcs_in.world_axis_object_components = [('spat', 'ra', 'ra.degree'),
                                           ('spat', 'dec', 'dec.degree'),
                                           ('spec', 0, 'value'),
                                           ('time', 0, 'utc.value')]
    wcs_in.world_axis_object_classes = {'spat': ('astropy.coordinates.SkyCoord', (),
                                                 {'frame': 'icrs'}),
                                        'spec': ('astropy.units.Wavelength', (None,), {}),
                                        'time': ('astropy.time.Time', (None,),
                                                 {'format': 'mjd', 'scale': 'utc'})}

    wcs_out = FakeWCS()
    wcs_out.low_level_wcs = wcs_out
    wcs_out.pixel_n_dim = 4
    wcs_out.world_n_dim = 4
    wcs_out.axis_correlation_matrix = [[True, False, False, False],
                                       [False, True, True, False],
                                       [False, True, True, False],
                                       [False, False, False, True]]
    wcs_out.world_axis_object_components = [('spec', 0, 'value'),
                                            ('spat', 'ra', 'ra.degree'),
                                            ('spat', 'dec', 'dec.degree'),
                                            ('time', 0, 'utc.value')]
    wcs_out.world_axis_object_classes = wcs_in.world_axis_object_classes

    matrix = pixel_to_pixel_correlation_matrix(wcs_in, wcs_out)

    matrix = matrix.astype(int)

    # The shape should be (n_pixel_out, n_pixel_in)
    assert matrix.shape == (4, 3)

    expected = np.array([[1, 1, 0], [1, 1, 0], [1, 1, 0], [0, 0, 1]])
    assert_equal(matrix, expected)


def test_split_matrix():

    assert split_matrix(np.array([[1]])) == [([0], [0])]

    assert split_matrix(np.array([[1, 1],
                                  [1, 1]])) == [([0, 1], [0, 1])]

    assert split_matrix(np.array([[1, 1, 0],
                                  [1, 1, 0],
                                  [0, 0, 1]])) == [([0, 1], [0, 1]), ([2], [2])]

    assert split_matrix(np.array([[0, 1, 0],
                                  [1, 0, 0],
                                  [0, 0, 1]])) == [([0], [1]), ([1], [0]), ([2], [2])]

    assert split_matrix(np.array([[0, 1, 1],
                                  [1, 0, 0],
                                  [1, 0, 1]])) == [([0, 1, 2], [0, 1, 2])]


def test_efficient_pixel_to_pixel():

    wcs_in = WCS(naxis=3)
    wcs_in.wcs.ctype = 'DEC--TAN', 'FREQ', 'RA---TAN'
    wcs_in.wcs.set()

    wcs_out = WCS(naxis=3)
    wcs_out.wcs.ctype = 'GLON-CAR', 'GLAT-CAR', 'FREQ'
    wcs_out.wcs.set()

    # First try with scalars
    x, y, z = pixel_to_pixel(wcs_in, wcs_out, 1, 2, 3)
    assert x.shape == ()
    assert y.shape == ()
    assert z.shape == ()

    # Now try with broadcasted arrays
    x = np.linspace(10, 20, 10)
    y = np.linspace(10, 20, 20)
    z = np.linspace(10, 20, 30)
    Z1, Y1, X1 = np.meshgrid(z, y, x, indexing='ij', copy=False)
    X2, Y2, Z2 = pixel_to_pixel(wcs_in, wcs_out, X1, Y1, Z1)

    # The final arrays should have the correct shape
    assert X2.shape == (30, 20, 10)
    assert Y2.shape == (30, 20, 10)
    assert Z2.shape == (30, 20, 10)

    # But behind the scenes should also be broadcasted
    assert unbroadcast(X2).shape == (30, 1, 10)
    assert unbroadcast(Y2).shape == (30, 1, 10)
    assert unbroadcast(Z2).shape == (20, 1)

    # We can put the values back through the function to ensure round-tripping
    X3, Y3, Z3 = pixel_to_pixel(wcs_out, wcs_in, X2, Y2, Z2)

    # The final arrays should have the correct shape
    assert X2.shape == (30, 20, 10)
    assert Y2.shape == (30, 20, 10)
    assert Z2.shape == (30, 20, 10)

    # But behind the scenes should also be broadcasted
    assert unbroadcast(X3).shape == (30, 1, 10)
    assert unbroadcast(Y3).shape == (20, 1)
    assert unbroadcast(Z3).shape == (30, 1, 10)

    # And these arrays should match the input
    assert_allclose(X1, X3)
    assert_allclose(Y1, Y3)
    assert_allclose(Z1, Z3)


def test_efficient_pixel_to_pixel_simple():

    # Simple test to make sure that when WCS only returns one world coordinate
    # this still works correctly (since this requires special treatment behind
    # the scenes).

    wcs_in = WCS(naxis=2)
    wcs_in.wcs.ctype = 'DEC--TAN', 'RA---TAN'
    wcs_in.wcs.set()

    wcs_out = WCS(naxis=2)
    wcs_out.wcs.ctype = 'GLON-CAR', 'GLAT-CAR'
    wcs_out.wcs.set()

    # First try with scalars
    x, y = pixel_to_pixel(wcs_in, wcs_out, 1, 2)
    assert x.shape == ()
    assert y.shape == ()

    # Now try with broadcasted arrays
    x = np.linspace(10, 20, 10)
    y = np.linspace(10, 20, 20)
    Y1, X1 = np.meshgrid(y, x, indexing='ij', copy=False)
    Y2, X2 = pixel_to_pixel(wcs_in, wcs_out, X1, Y1)

    # The final arrays should have the correct shape
    assert X2.shape == (20, 10)
    assert Y2.shape == (20, 10)

    # and there are no efficiency gains here since the celestial axes are correlated
    assert unbroadcast(X2).shape == (20, 10)

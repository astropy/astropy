import pytest

from numpy.testing import assert_equal, assert_allclose
from astropy.wcs import WCS
from astropy.io.fits import Header
from astropy.coordinates import SkyCoord, Galactic
from astropy.units import Quantity
from astropy.wcs.wcsapi.sliced_low_level_wcs import SlicedLowLevelWCS, sanitize_slices
import astropy.units as u

# To test the slicing we start off from standard FITS WCS
# objects since those implement the low-level API. We create
# a WCS for a spectral cube with axes in non-standard order
# and with correlated celestial axes and an uncorrelated
# spectral axis.

HEADER_SPECTRAL_CUBE = """
NAXIS   = 3
NAXIS1  = 10
NAXIS2  = 20
NAXIS3  = 30
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

WCS_SPECTRAL_CUBE = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE, sep='\n'))
WCS_SPECTRAL_CUBE.pixel_bounds = [(-1, 11), (-2, 18), (5, 15)]


@pytest.mark.parametrize("item, ndim, expected", (
    ([Ellipsis, 10], 4, [slice(None)] * 3 + [10]),
    ([10, slice(20, 30)], 5, [10, slice(20, 30)] + [slice(None)] * 3),
    ([10, Ellipsis, 8], 10, [10] + [slice(None)] * 8 + [8])
))
def test_sanitize_slice(item, ndim, expected):
    new_item = sanitize_slices(item, ndim)
    # FIXME: do we still need the first two since the third assert
    # should cover it all?
    assert len(new_item) == ndim
    assert all(isinstance(i, (slice, int)) for i in new_item)
    assert new_item == expected


EXPECTED_ELLIPSIS_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20, 10)

Pixel Dim  Data size  Bounds
        0         10  (-1, 11)
        1         20  (-2, 18)
        2         30  (5, 15)

World Dim  Physical Type     Units
        0  pos.galactic.lat  deg
        1  em.freq           Hz
        2  pos.galactic.lon  deg

Correlation between pixel and world axes:

             Pixel Dim
World Dim    0    1    2
        0  yes   no  yes
        1   no  yes   no
        2  yes   no  yes
"""


def test_ellipsis():

    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, Ellipsis)

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape == (30, 20, 10)
    assert wcs.pixel_shape == (10, 20, 30)
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[True, False, True], [False, True, False], [True, False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                  ('freq', 0, 'value'),
                                                  ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert wcs.world_axis_object_classes['freq'][0] is Quantity
    assert wcs.world_axis_object_classes['freq'][1] == ()
    assert wcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(wcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29., 39., 44.))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))

    assert_equal(wcs.pixel_bounds, [(-1, 11), (-2, 18), (5, 15)])

    assert str(wcs) == repr(wcs) == EXPECTED_ELLIPSIS_REPR.strip()


EXPECTED_SPECTRAL_SLICE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 2 pixel and 2 world dimensions

Array shape (Numpy order): (30, 10)

Pixel Dim  Data size  Bounds
        0         10  (-1, 11)
        1         30  (5, 15)

World Dim  Physical Type     Units
        0  pos.galactic.lat  deg
        1  pos.galactic.lon  deg

Correlation between pixel and world axes:

           Pixel Dim
World Dim    0    1
        0  yes  yes
        1  yes  yes
"""

def test_spectral_slice():

    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, [slice(None), 10])

    assert wcs.pixel_n_dim == 2
    assert wcs.world_n_dim == 2
    assert wcs.array_shape == (30, 10)
    assert wcs.pixel_shape == (10, 30)
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[True, True], [True, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert_allclose(wcs.pixel_to_world_values(29, 44), (10, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 29), (10, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 25), (29., 44.))
    assert_equal(wcs.world_to_array_index_values(10, 25), (44, 29))

    assert_equal(wcs.pixel_bounds, [(-1, 11), (5, 15)])

    assert str(wcs) == repr(wcs) == EXPECTED_SPECTRAL_SLICE_REPR.strip()


EXPECTED_SPECTRAL_RANGE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 6, 10)

Pixel Dim  Data size  Bounds
        0         10  (-1, 11)
        1          6  (-6, 14)
        2         30  (5, 15)

World Dim  Physical Type     Units
        0  pos.galactic.lat  deg
        1  em.freq           Hz
        2  pos.galactic.lon  deg

Correlation between pixel and world axes:

             Pixel Dim
World Dim    0    1    2
        0  yes   no  yes
        1   no  yes   no
        2  yes   no  yes
"""


def test_spectral_range():

    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, [slice(None), slice(4, 10)])

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape == (30, 6, 10)
    assert wcs.pixel_shape == (10, 6, 30)
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[True, False, True], [False, True, False], [True, False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                  ('freq', 0, 'value'),
                                                  ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert wcs.world_axis_object_classes['freq'][0] is Quantity
    assert wcs.world_axis_object_classes['freq'][1] == ()
    assert wcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(wcs.pixel_to_world_values(29, 35, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 35, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29., 35., 44.))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 35, 29))

    assert_equal(wcs.pixel_bounds, [(-1, 11), (-6, 14), (5, 15)])

    assert str(wcs) == repr(wcs) == EXPECTED_SPECTRAL_RANGE_REPR.strip()


EXPECTED_CELESTIAL_SLICE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 2 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20)

Pixel Dim  Data size  Bounds
        0         20  (-2, 18)
        1         30  (5, 15)

World Dim  Physical Type     Units
        0  pos.galactic.lat  deg
        1  em.freq           Hz
        2  pos.galactic.lon  deg

Correlation between pixel and world axes:

           Pixel Dim
World Dim    0    1
        0   no  yes
        1  yes   no
        2   no  yes
"""


def test_celestial_slice():

    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, [Ellipsis, 5])

    assert wcs.pixel_n_dim == 2
    assert wcs.world_n_dim == 3
    assert wcs.array_shape == (30, 20)
    assert wcs.pixel_shape == (20, 30)
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[False, True], [True, False], [False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                ('freq', 0, 'value'),
                                                ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert wcs.world_axis_object_classes['freq'][0] is Quantity
    assert wcs.world_axis_object_classes['freq'][1] == ()
    assert wcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(wcs.pixel_to_world_values(39, 44), (12.4, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39), (12.4, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(12.4, 20, 25), (39., 44.))
    assert_equal(wcs.world_to_array_index_values(12.4, 20, 25), (44, 39))

    assert_equal(wcs.pixel_bounds, [(-2, 18), (5, 15)])

    assert str(wcs) == repr(wcs) == EXPECTED_CELESTIAL_SLICE_REPR.strip()


EXPECTED_CELESTIAL_RANGE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20, 5)

Pixel Dim  Data size  Bounds
        0          5  (-6, 6)
        1         20  (-2, 18)
        2         30  (5, 15)

World Dim  Physical Type     Units
        0  pos.galactic.lat  deg
        1  em.freq           Hz
        2  pos.galactic.lon  deg

Correlation between pixel and world axes:

             Pixel Dim
World Dim    0    1    2
        0  yes   no  yes
        1   no  yes   no
        2  yes   no  yes
"""


def test_celestial_range():

    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, [Ellipsis, slice(5, 10)])

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape == (30, 20, 5)
    assert wcs.pixel_shape == (5, 20, 30)
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[True, False, True], [False, True, False], [True, False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                ('freq', 0, 'value'),
                                                ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert wcs.world_axis_object_classes['freq'][0] is Quantity
    assert wcs.world_axis_object_classes['freq'][1] == ()
    assert wcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(wcs.pixel_to_world_values(24, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 24), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (24., 39., 44.))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 24))

    assert_equal(wcs.pixel_bounds, [(-6, 6), (-2, 18), (5, 15)])

    assert str(wcs) == repr(wcs) == EXPECTED_CELESTIAL_RANGE_REPR.strip()


# Now try with a 90 degree rotation

WCS_SPECTRAL_CUBE_ROT = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE, sep='\n'))
WCS_SPECTRAL_CUBE_ROT.wcs.pc = [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
WCS_SPECTRAL_CUBE_ROT.wcs.crval[0] = 0
WCS_SPECTRAL_CUBE_ROT.pixel_bounds = [(-1, 11), (-2, 18), (5, 15)]

EXPECTED_CELESTIAL_RANGE_ROT_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20, 5)

Pixel Dim  Data size  Bounds
        0          5  (-6, 6)
        1         20  (-2, 18)
        2         30  (5, 15)

World Dim  Physical Type     Units
        0  pos.galactic.lat  deg
        1  em.freq           Hz
        2  pos.galactic.lon  deg

Correlation between pixel and world axes:

             Pixel Dim
World Dim    0    1    2
        0  yes   no  yes
        1   no  yes   no
        2  yes   no  yes
"""

def test_celestial_range_rot():

    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE_ROT, [Ellipsis, slice(5, 10)])

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape == (30, 20, 5)
    assert wcs.pixel_shape == (5, 20, 30)
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[True, False, True], [False, True, False], [True, False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                ('freq', 0, 'value'),
                                                ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert wcs.world_axis_object_classes['freq'][0] is Quantity
    assert wcs.world_axis_object_classes['freq'][1] == ()
    assert wcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(wcs.pixel_to_world_values(14, 29, 34), (1, 15, 24))
    assert_allclose(wcs.array_index_to_world_values(34, 29, 14), (1, 15, 24))

    assert_allclose(wcs.world_to_pixel_values(1, 15, 24), (14., 29., 34.))
    assert_equal(wcs.world_to_array_index_values(1, 15, 24), (34, 29, 14))

    assert_equal(wcs.pixel_bounds, [(-6, 6), (-2, 18), (5, 15)])

    assert str(wcs) == repr(wcs) == EXPECTED_CELESTIAL_RANGE_ROT_REPR.strip()


HEADER_NO_SHAPE_CUBE = """
NAXIS   = 3
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

WCS_NO_SHAPE_CUBE = WCS(Header.fromstring(HEADER_NO_SHAPE_CUBE, sep='\n'))

EXPECTED_NO_SHAPE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): None

Pixel Dim  Data size  Bounds
        0       None  None
        1       None  None
        2       None  None

World Dim  Physical Type     Units
        0  pos.galactic.lat  deg
        1  em.freq           Hz
        2  pos.galactic.lon  deg

Correlation between pixel and world axes:

             Pixel Dim
World Dim    0    1    2
        0  yes   no  yes
        1   no  yes   no
        2  yes   no  yes
"""


def test_no_array_shape():

    wcs = SlicedLowLevelWCS(WCS_NO_SHAPE_CUBE, Ellipsis)

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape is None
    assert wcs.pixel_shape is None
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[True, False, True], [False, True, False], [True, False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                ('freq', 0, 'value'),
                                                ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert wcs.world_axis_object_classes['freq'][0] is Quantity
    assert wcs.world_axis_object_classes['freq'][1] == ()
    assert wcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(wcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29., 39., 44.))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))

    assert str(wcs) == repr(wcs) == EXPECTED_NO_SHAPE_REPR.strip()


# Testing the WCS object having some physical types as None/Unknown
HEADER_SPECTRAL_CUBE_NONE_TYPES = {
        'CTYPE1': 'GLAT-CAR',
        'CUNIT1': 'deg',
        'CDELT1': -0.1,
        'CRPIX1': 30,
        'CRVAL1': 10,
        'NAXIS1': 10,
        'CTYPE2': '',
        'CUNIT2': 'Hz',
        'CDELT2':  0.5,
        'CRPIX2': 40,
        'CRVAL2': 20,
        'NAXIS2': 20,
        'CTYPE3': 'GLON-CAR',
        'CUNIT3': 'deg',
        'CDELT3':  0.1,
        'CRPIX3': 45,
        'CRVAL3': 25,
        'NAXIS3': 30
}

WCS_SPECTRAL_CUBE_NONE_TYPES = WCS(header=HEADER_SPECTRAL_CUBE_NONE_TYPES)
WCS_SPECTRAL_CUBE_NONE_TYPES.pixel_bounds = [(-1, 11), (-2, 18), (5, 15)]


EXPECTED_ELLIPSIS_REPR_NONE_TYPES = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20, 10)

Pixel Dim  Data size  Bounds
        0         10  (-1, 11)
        1         20  (-2, 18)
        2         30  (5, 15)

World Dim  Physical Type     Units
        0  pos.galactic.lat  deg
        1  None              unknown
        2  pos.galactic.lon  deg

Correlation between pixel and world axes:

             Pixel Dim
World Dim    0    1    2
        0  yes   no  yes
        1   no  yes   no
        2  yes   no  yes
"""


def test_ellipsis_none_types():

    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE_NONE_TYPES, Ellipsis)

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape == (30, 20, 10)
    assert wcs.pixel_shape == (10, 20, 30)
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', None, 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[True, False, True],
                                               [False, True, False], [True, False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                  ('world', 0, 'value'),
                                                  ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert_allclose(wcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29., 39., 44.))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))

    assert_equal(wcs.pixel_bounds, [(-1, 11), (-2, 18), (5, 15)])

    assert str(wcs) == repr(wcs) == EXPECTED_ELLIPSIS_REPR_NONE_TYPES.strip()

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

import astropy.units as u
from astropy.coordinates import ICRS, Galactic, SkyCoord
from astropy.io.fits import Header
from astropy.io.fits.verify import VerifyWarning
from astropy.time import Time
from astropy.units import Quantity
from astropy.wcs.wcs import WCS, FITSFixedWarning
from astropy.wcs.wcsapi.wrappers.sliced_wcs import (
    SlicedLowLevelWCS,
    combine_slices,
    sanitize_slices,
)

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
CNAME1  = Latitude
CNAME2  = Frequency
CNAME3  = Longitude
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

with warnings.catch_warnings():
    warnings.simplefilter("ignore", VerifyWarning)
    WCS_SPECTRAL_CUBE = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE, sep="\n"))
WCS_SPECTRAL_CUBE.pixel_bounds = [(-1, 11), (-2, 18), (5, 15)]


def test_invalid_slices():
    with pytest.raises(IndexError):
        SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, [None, None, [False, False, False]])

    with pytest.raises(IndexError):
        SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, [None, None, slice(None, None, 2)])

    with pytest.raises(IndexError):
        SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, [None, None, 1000.100])


@pytest.mark.parametrize(
    "item, ndim, expected",
    (
        ([Ellipsis, 10], 4, [slice(None)] * 3 + [10]),
        ([10, slice(20, 30)], 5, [10, slice(20, 30)] + [slice(None)] * 3),
        ([10, Ellipsis, 8], 10, [10] + [slice(None)] * 8 + [8]),
    ),
)
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

Pixel Dim  Axis Name  Data size  Bounds
        0  None              10  (-1, 11)
        1  None              20  (-2, 18)
        2  None              30  (5, 15)

World Dim  Axis Name  Physical Type     Units
        0  Latitude   pos.galactic.lat  deg
        1  Frequency  em.freq           Hz
        2  Longitude  pos.galactic.lon  deg

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
    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        "em.freq",
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]
    assert wcs.pixel_axis_names == ["", "", ""]
    assert wcs.world_axis_names == ["Latitude", "Frequency", "Longitude"]

    assert_equal(
        wcs.axis_correlation_matrix,
        [[True, False, True], [False, True, False], [True, False, True]],
    )

    assert len(wcs.world_axis_object_components) == 3
    assert wcs.world_axis_object_components[0] == (
        "celestial",
        1,
        "spherical.lat.degree",
    )
    assert wcs.world_axis_object_components[1][:2] == ("spectral", 0)
    assert wcs.world_axis_object_components[2] == (
        "celestial",
        0,
        "spherical.lon.degree",
    )

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] is u.deg

    assert wcs.world_axis_object_classes["spectral"][0] is Quantity
    assert wcs.world_axis_object_classes["spectral"][1] == ()
    assert wcs.world_axis_object_classes["spectral"][2] == {}

    assert_allclose(wcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29.0, 39.0, 44.0))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))

    assert_equal(wcs.pixel_bounds, [(-1, 11), (-2, 18), (5, 15)])

    assert str(wcs) == EXPECTED_ELLIPSIS_REPR.strip()
    assert EXPECTED_ELLIPSIS_REPR.strip() in repr(wcs)


def test_pixel_to_world_broadcasting():
    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, Ellipsis)

    assert_allclose(
        wcs.pixel_to_world_values((29, 29), 39, 44), ((10, 10), (20, 20), (25, 25))
    )


def test_world_to_pixel_broadcasting():
    wcs = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, Ellipsis)

    assert_allclose(
        wcs.world_to_pixel_values((10, 10), 20, 25),
        ((29.0, 29.0), (39.0, 39.0), (44.0, 44.0)),
    )


EXPECTED_SPECTRAL_SLICE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 2 pixel and 2 world dimensions

Array shape (Numpy order): (30, 10)

Pixel Dim  Axis Name  Data size  Bounds
        0  None              10  (-1, 11)
        1  None              30  (5, 15)

World Dim  Axis Name  Physical Type     Units
        0  Latitude   pos.galactic.lat  deg
        1  Longitude  pos.galactic.lon  deg

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
    assert wcs.world_axis_physical_types == ["pos.galactic.lat", "pos.galactic.lon"]
    assert wcs.world_axis_units == ["deg", "deg"]
    assert wcs.pixel_axis_names == ["", ""]
    assert wcs.world_axis_names == ["Latitude", "Longitude"]

    assert_equal(wcs.axis_correlation_matrix, [[True, True], [True, True]])

    assert wcs.world_axis_object_components == [
        ("celestial", 1, "spherical.lat.degree"),
        ("celestial", 0, "spherical.lon.degree"),
    ]

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] is u.deg

    assert_allclose(wcs.pixel_to_world_values(29, 44), (10, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 29), (10, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 25), (29.0, 44.0))
    assert_equal(wcs.world_to_array_index_values(10, 25), (44, 29))

    assert_equal(wcs.pixel_bounds, [(-1, 11), (5, 15)])

    assert str(wcs) == EXPECTED_SPECTRAL_SLICE_REPR.strip()
    assert EXPECTED_SPECTRAL_SLICE_REPR.strip() in repr(wcs)


EXPECTED_SPECTRAL_RANGE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 6, 10)

Pixel Dim  Axis Name  Data size  Bounds
        0  None              10  (-1, 11)
        1  None               6  (-6, 14)
        2  None              30  (5, 15)

World Dim  Axis Name  Physical Type     Units
        0  Latitude   pos.galactic.lat  deg
        1  Frequency  em.freq           Hz
        2  Longitude  pos.galactic.lon  deg

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
    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        "em.freq",
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]
    assert wcs.pixel_axis_names == ["", "", ""]
    assert wcs.world_axis_names == ["Latitude", "Frequency", "Longitude"]

    assert_equal(
        wcs.axis_correlation_matrix,
        [[True, False, True], [False, True, False], [True, False, True]],
    )

    assert len(wcs.world_axis_object_components) == 3
    assert wcs.world_axis_object_components[0] == (
        "celestial",
        1,
        "spherical.lat.degree",
    )
    assert wcs.world_axis_object_components[1][:2] == ("spectral", 0)
    assert wcs.world_axis_object_components[2] == (
        "celestial",
        0,
        "spherical.lon.degree",
    )

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] is u.deg

    assert wcs.world_axis_object_classes["spectral"][0] is Quantity
    assert wcs.world_axis_object_classes["spectral"][1] == ()
    assert wcs.world_axis_object_classes["spectral"][2] == {}

    assert_allclose(wcs.pixel_to_world_values(29, 35, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 35, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29.0, 35.0, 44.0))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 35, 29))

    assert_equal(wcs.pixel_bounds, [(-1, 11), (-6, 14), (5, 15)])

    assert str(wcs) == EXPECTED_SPECTRAL_RANGE_REPR.strip()
    assert EXPECTED_SPECTRAL_RANGE_REPR.strip() in repr(wcs)


EXPECTED_CELESTIAL_SLICE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 2 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20)

Pixel Dim  Axis Name  Data size  Bounds
        0  None              20  (-2, 18)
        1  None              30  (5, 15)

World Dim  Axis Name  Physical Type     Units
        0  Latitude   pos.galactic.lat  deg
        1  Frequency  em.freq           Hz
        2  Longitude  pos.galactic.lon  deg

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
    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        "em.freq",
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]
    assert wcs.pixel_axis_names == ["", ""]
    assert wcs.world_axis_names == ["Latitude", "Frequency", "Longitude"]

    assert_equal(
        wcs.axis_correlation_matrix, [[False, True], [True, False], [False, True]]
    )

    assert len(wcs.world_axis_object_components) == 3
    assert wcs.world_axis_object_components[0] == (
        "celestial",
        1,
        "spherical.lat.degree",
    )
    assert wcs.world_axis_object_components[1][:2] == ("spectral", 0)
    assert wcs.world_axis_object_components[2] == (
        "celestial",
        0,
        "spherical.lon.degree",
    )

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] is u.deg

    assert wcs.world_axis_object_classes["spectral"][0] is Quantity
    assert wcs.world_axis_object_classes["spectral"][1] == ()
    assert wcs.world_axis_object_classes["spectral"][2] == {}

    assert_allclose(wcs.pixel_to_world_values(39, 44), (12.4, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39), (12.4, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(12.4, 20, 25), (39.0, 44.0))
    assert_equal(wcs.world_to_array_index_values(12.4, 20, 25), (44, 39))

    assert_equal(wcs.pixel_bounds, [(-2, 18), (5, 15)])

    assert str(wcs) == EXPECTED_CELESTIAL_SLICE_REPR.strip()
    assert EXPECTED_CELESTIAL_SLICE_REPR.strip() in repr(wcs)


EXPECTED_CELESTIAL_RANGE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20, 5)

Pixel Dim  Axis Name  Data size  Bounds
        0  None               5  (-6, 6)
        1  None              20  (-2, 18)
        2  None              30  (5, 15)

World Dim  Axis Name  Physical Type     Units
        0  Latitude   pos.galactic.lat  deg
        1  Frequency  em.freq           Hz
        2  Longitude  pos.galactic.lon  deg

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
    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        "em.freq",
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]
    assert wcs.pixel_axis_names == ["", "", ""]
    assert wcs.world_axis_names == ["Latitude", "Frequency", "Longitude"]

    assert_equal(
        wcs.axis_correlation_matrix,
        [[True, False, True], [False, True, False], [True, False, True]],
    )

    assert len(wcs.world_axis_object_components) == 3
    assert wcs.world_axis_object_components[0] == (
        "celestial",
        1,
        "spherical.lat.degree",
    )
    assert wcs.world_axis_object_components[1][:2] == ("spectral", 0)
    assert wcs.world_axis_object_components[2] == (
        "celestial",
        0,
        "spherical.lon.degree",
    )

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] is u.deg

    assert wcs.world_axis_object_classes["spectral"][0] is Quantity
    assert wcs.world_axis_object_classes["spectral"][1] == ()
    assert wcs.world_axis_object_classes["spectral"][2] == {}

    assert_allclose(wcs.pixel_to_world_values(24, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 24), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (24.0, 39.0, 44.0))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 24))

    assert_equal(wcs.pixel_bounds, [(-6, 6), (-2, 18), (5, 15)])

    assert str(wcs) == EXPECTED_CELESTIAL_RANGE_REPR.strip()
    assert EXPECTED_CELESTIAL_RANGE_REPR.strip() in repr(wcs)


# Now try with a 90 degree rotation

with warnings.catch_warnings():
    warnings.simplefilter("ignore", VerifyWarning)
    WCS_SPECTRAL_CUBE_ROT = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE, sep="\n"))
WCS_SPECTRAL_CUBE_ROT.wcs.pc = [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
WCS_SPECTRAL_CUBE_ROT.wcs.crval[0] = 0
WCS_SPECTRAL_CUBE_ROT.pixel_bounds = [(-1, 11), (-2, 18), (5, 15)]

EXPECTED_CELESTIAL_RANGE_ROT_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20, 5)

Pixel Dim  Axis Name  Data size  Bounds
        0  None               5  (-6, 6)
        1  None              20  (-2, 18)
        2  None              30  (5, 15)

World Dim  Axis Name  Physical Type     Units
        0  Latitude   pos.galactic.lat  deg
        1  Frequency  em.freq           Hz
        2  Longitude  pos.galactic.lon  deg

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
    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        "em.freq",
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]
    assert wcs.pixel_axis_names == ["", "", ""]
    assert wcs.world_axis_names == ["Latitude", "Frequency", "Longitude"]

    assert_equal(
        wcs.axis_correlation_matrix,
        [[True, False, True], [False, True, False], [True, False, True]],
    )

    assert len(wcs.world_axis_object_components) == 3
    assert wcs.world_axis_object_components[0] == (
        "celestial",
        1,
        "spherical.lat.degree",
    )
    assert wcs.world_axis_object_components[1][:2] == ("spectral", 0)
    assert wcs.world_axis_object_components[2] == (
        "celestial",
        0,
        "spherical.lon.degree",
    )

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] is u.deg

    assert wcs.world_axis_object_classes["spectral"][0] is Quantity
    assert wcs.world_axis_object_classes["spectral"][1] == ()
    assert wcs.world_axis_object_classes["spectral"][2] == {}

    assert_allclose(wcs.pixel_to_world_values(14, 29, 34), (1, 15, 24))
    assert_allclose(wcs.array_index_to_world_values(34, 29, 14), (1, 15, 24))

    assert_allclose(wcs.world_to_pixel_values(1, 15, 24), (14.0, 29.0, 34.0))
    assert_equal(wcs.world_to_array_index_values(1, 15, 24), (34, 29, 14))

    assert_equal(wcs.pixel_bounds, [(-6, 6), (-2, 18), (5, 15)])

    assert str(wcs) == EXPECTED_CELESTIAL_RANGE_ROT_REPR.strip()
    assert EXPECTED_CELESTIAL_RANGE_ROT_REPR.strip() in repr(wcs)


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

with warnings.catch_warnings():
    warnings.simplefilter("ignore", VerifyWarning)
    WCS_NO_SHAPE_CUBE = WCS(Header.fromstring(HEADER_NO_SHAPE_CUBE, sep="\n"))

EXPECTED_NO_SHAPE_REPR = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): None

Pixel Dim  Axis Name  Data size  Bounds
        0  None            None  None
        1  None            None  None
        2  None            None  None

World Dim  Axis Name  Physical Type     Units
        0  None       pos.galactic.lat  deg
        1  None       em.freq           Hz
        2  None       pos.galactic.lon  deg

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
    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        "em.freq",
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]

    assert_equal(
        wcs.axis_correlation_matrix,
        [[True, False, True], [False, True, False], [True, False, True]],
    )

    assert len(wcs.world_axis_object_components) == 3
    assert wcs.world_axis_object_components[0] == (
        "celestial",
        1,
        "spherical.lat.degree",
    )
    assert wcs.world_axis_object_components[1][:2] == ("spectral", 0)
    assert wcs.world_axis_object_components[2] == (
        "celestial",
        0,
        "spherical.lon.degree",
    )

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] is u.deg

    assert wcs.world_axis_object_classes["spectral"][0] is Quantity
    assert wcs.world_axis_object_classes["spectral"][1] == ()
    assert wcs.world_axis_object_classes["spectral"][2] == {}

    assert_allclose(wcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29.0, 39.0, 44.0))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))

    assert str(wcs) == EXPECTED_NO_SHAPE_REPR.strip()
    assert EXPECTED_NO_SHAPE_REPR.strip() in repr(wcs)


# Testing the WCS object having some physical types as None/Unknown
HEADER_SPECTRAL_CUBE_NONE_TYPES = {
    "CTYPE1": "GLAT-CAR",
    "CUNIT1": "deg",
    "CDELT1": -0.1,
    "CRPIX1": 30,
    "CRVAL1": 10,
    "NAXIS1": 10,
    "CTYPE2": "",
    "CUNIT2": "Hz",
    "CDELT2": 0.5,
    "CRPIX2": 40,
    "CRVAL2": 20,
    "NAXIS2": 20,
    "CTYPE3": "GLON-CAR",
    "CUNIT3": "deg",
    "CDELT3": 0.1,
    "CRPIX3": 45,
    "CRVAL3": 25,
    "NAXIS3": 30,
}

WCS_SPECTRAL_CUBE_NONE_TYPES = WCS(header=HEADER_SPECTRAL_CUBE_NONE_TYPES)
WCS_SPECTRAL_CUBE_NONE_TYPES.pixel_bounds = [(-1, 11), (-2, 18), (5, 15)]


EXPECTED_ELLIPSIS_REPR_NONE_TYPES = """
SlicedLowLevelWCS Transformation

This transformation has 3 pixel and 3 world dimensions

Array shape (Numpy order): (30, 20, 10)

Pixel Dim  Axis Name  Data size  Bounds
        0  None              10  (-1, 11)
        1  None              20  (-2, 18)
        2  None              30  (5, 15)

World Dim  Axis Name  Physical Type     Units
        0  None       pos.galactic.lat  deg
        1  None       None              Hz
        2  None       pos.galactic.lon  deg

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
    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        None,
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]

    assert_equal(
        wcs.axis_correlation_matrix,
        [[True, False, True], [False, True, False], [True, False, True]],
    )

    assert wcs.world_axis_object_components == [
        ("celestial", 1, "spherical.lat.degree"),
        ("world", 0, "value"),
        ("celestial", 0, "spherical.lon.degree"),
    ]

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] is u.deg

    assert_allclose(wcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29.0, 39.0, 44.0))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))

    assert_equal(wcs.pixel_bounds, [(-1, 11), (-2, 18), (5, 15)])

    assert str(wcs) == EXPECTED_ELLIPSIS_REPR_NONE_TYPES.strip()
    assert EXPECTED_ELLIPSIS_REPR_NONE_TYPES.strip() in repr(wcs)


CASES = [
    (slice(None), slice(None), slice(None)),
    (slice(None), slice(3, None), slice(3, None)),
    (slice(None), slice(None, 16), slice(None, 16)),
    (slice(None), slice(3, 16), slice(3, 16)),
    (slice(2, None), slice(None), slice(2, None)),
    (slice(2, None), slice(3, None), slice(5, None)),
    (slice(2, None), slice(None, 16), slice(2, 18)),
    (slice(2, None), slice(3, 16), slice(5, 18)),
    (slice(None, 10), slice(None), slice(None, 10)),
    (slice(None, 10), slice(3, None), slice(3, 10)),
    (slice(None, 10), slice(None, 16), slice(None, 10)),
    (slice(None, 10), slice(3, 16), slice(3, 10)),
    (slice(2, 10), slice(None), slice(2, 10)),
    (slice(2, 10), slice(3, None), slice(5, 10)),
    (slice(2, 10), slice(None, 16), slice(2, 10)),
    (slice(2, 10), slice(3, 16), slice(5, 10)),
    (slice(None), 3, 3),
    (slice(2, None), 3, 5),
    (slice(None, 10), 3, 3),
    (slice(2, 10), 3, 5),
]


@pytest.mark.parametrize(("slice1", "slice2", "expected"), CASES)
def test_combine_slices(slice1, slice2, expected):
    assert combine_slices(slice1, slice2) == expected


def test_nested_slicing():
    # Make sure that if we call slicing several times, the result is the same
    # as calling the slicing once with the final slice settings.
    wcs = WCS_SPECTRAL_CUBE
    sub1 = SlicedLowLevelWCS(
        SlicedLowLevelWCS(
            SlicedLowLevelWCS(wcs, [slice(None), slice(1, 10), slice(None)]),
            [3, slice(2, None)],
        ),
        [slice(None), slice(2, 8)],
    )
    sub2 = wcs[3, 3:10, 2:8]
    assert_allclose(sub1.pixel_to_world_values(3, 5), sub2.pixel_to_world_values(3, 5))
    assert not isinstance(sub1._wcs, SlicedLowLevelWCS)


def test_too_much_slicing():
    wcs = WCS_SPECTRAL_CUBE
    with pytest.raises(
        ValueError,
        match=(
            "Cannot slice WCS: the resulting WCS "
            "should have at least one pixel and "
            "one world dimension"
        ),
    ):
        wcs[0, 1, 2]


HEADER_TIME_1D = """
SIMPLE  = T
BITPIX  = -32
NAXIS   = 1
NAXIS1  = 2048
TIMESYS = 'UTC'
TREFPOS = 'TOPOCENT'
MJDREF  = 50002.6
CTYPE1  = 'UTC'
CRVAL1  = 5
CUNIT1  = 's'
CRPIX1  = 1.0
CDELT1  = 2
OBSGEO-L= -20
OBSGEO-B= -70
OBSGEO-H= 2530
"""


@pytest.fixture
def header_time_1d():
    return Header.fromstring(HEADER_TIME_1D, sep="\n")


@pytest.fixture
def time_1d_wcs(header_time_1d):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(header_time_1d)


def test_1d_sliced_low_level(time_1d_wcs):
    sll = SlicedLowLevelWCS(time_1d_wcs, np.s_[10:20])
    world = sll.pixel_to_world_values([1, 2])
    assert isinstance(world, np.ndarray)
    assert np.allclose(world, [27, 29])


def validate_info_dict(result, expected):
    result_value = result.pop("value")
    expected_value = expected.pop("value")

    np.testing.assert_allclose(result_value, expected_value)
    assert result == expected


def test_dropped_dimensions():
    wcs = WCS_SPECTRAL_CUBE

    sub = SlicedLowLevelWCS(wcs, np.s_[:, :, :])

    assert sub.dropped_world_dimensions == {}

    sub = SlicedLowLevelWCS(wcs, np.s_[:, 2:5, :])

    assert sub.dropped_world_dimensions == {}

    sub = SlicedLowLevelWCS(wcs, np.s_[:, 0])

    waocomp = sub.dropped_world_dimensions.pop("world_axis_object_components")
    assert len(waocomp) == 1 and waocomp[0][0] == "spectral" and waocomp[0][1] == 0
    waocls = sub.dropped_world_dimensions.pop("world_axis_object_classes")
    assert (
        len(waocls) == 1
        and "spectral" in waocls
        and waocls["spectral"][0] == u.Quantity
    )
    validate_info_dict(
        sub.dropped_world_dimensions,
        {
            "value": [0.5],
            "world_axis_physical_types": ["em.freq"],
            "world_axis_names": ["Frequency"],
            "world_axis_units": ["Hz"],
            "serialized_classes": False,
        },
    )

    sub = SlicedLowLevelWCS(wcs, np.s_[:, 0, 0])

    waocomp = sub.dropped_world_dimensions.pop("world_axis_object_components")
    assert len(waocomp) == 1 and waocomp[0][0] == "spectral" and waocomp[0][1] == 0
    waocls = sub.dropped_world_dimensions.pop("world_axis_object_classes")
    assert (
        len(waocls) == 1
        and "spectral" in waocls
        and waocls["spectral"][0] == u.Quantity
    )
    validate_info_dict(
        sub.dropped_world_dimensions,
        {
            "value": [0.5],
            "world_axis_physical_types": ["em.freq"],
            "world_axis_names": ["Frequency"],
            "world_axis_units": ["Hz"],
            "serialized_classes": False,
        },
    )

    sub = SlicedLowLevelWCS(wcs, np.s_[0, :, 0])

    dwd = sub.dropped_world_dimensions
    wao_classes = dwd.pop("world_axis_object_classes")
    validate_info_dict(
        dwd,
        {
            "value": [12.86995801, 20.49217541],
            "world_axis_physical_types": ["pos.galactic.lat", "pos.galactic.lon"],
            "world_axis_names": ["Latitude", "Longitude"],
            "world_axis_units": ["deg", "deg"],
            "serialized_classes": False,
            "world_axis_object_components": [
                ("celestial", 1, "spherical.lat.degree"),
                ("celestial", 0, "spherical.lon.degree"),
            ],
        },
    )

    assert wao_classes["celestial"][0] is SkyCoord
    assert wao_classes["celestial"][1] == ()
    assert isinstance(wao_classes["celestial"][2]["frame"], Galactic)
    assert wao_classes["celestial"][2]["unit"] is u.deg

    sub = SlicedLowLevelWCS(wcs, np.s_[5, :5, 12])

    dwd = sub.dropped_world_dimensions
    wao_classes = dwd.pop("world_axis_object_classes")
    validate_info_dict(
        dwd,
        {
            "value": [11.67648267, 21.01921192],
            "world_axis_physical_types": ["pos.galactic.lat", "pos.galactic.lon"],
            "world_axis_names": ["Latitude", "Longitude"],
            "world_axis_units": ["deg", "deg"],
            "serialized_classes": False,
            "world_axis_object_components": [
                ("celestial", 1, "spherical.lat.degree"),
                ("celestial", 0, "spherical.lon.degree"),
            ],
        },
    )

    assert wao_classes["celestial"][0] is SkyCoord
    assert wao_classes["celestial"][1] == ()
    assert isinstance(wao_classes["celestial"][2]["frame"], Galactic)
    assert wao_classes["celestial"][2]["unit"] is u.deg


def test_dropped_dimensions_4d(cube_4d_fitswcs):
    sub = SlicedLowLevelWCS(cube_4d_fitswcs, np.s_[:, 12, 5, 5])

    dwd = sub.dropped_world_dimensions
    wao_classes = dwd.pop("world_axis_object_classes")
    wao_components = dwd.pop("world_axis_object_components")

    validate_info_dict(
        dwd,
        {
            "value": [4.0e00, -2.0e00, 1.0e10],
            "world_axis_physical_types": ["pos.eq.ra", "pos.eq.dec", "em.freq"],
            "world_axis_names": ["Right Ascension", "Declination", "Frequency"],
            "world_axis_units": ["deg", "deg", "Hz"],
            "serialized_classes": False,
        },
    )

    assert wao_classes["celestial"][0] is SkyCoord
    assert wao_classes["celestial"][1] == ()
    assert isinstance(wao_classes["celestial"][2]["frame"], ICRS)
    assert wao_classes["celestial"][2]["unit"] is u.deg
    assert wao_classes["spectral"][0:3] == (u.Quantity, (), {})

    assert wao_components[0] == ("celestial", 0, "spherical.lon.degree")
    assert wao_components[1] == ("celestial", 1, "spherical.lat.degree")
    assert wao_components[2][0:2] == ("spectral", 0)

    sub = SlicedLowLevelWCS(cube_4d_fitswcs, np.s_[12, 12])

    dwd = sub.dropped_world_dimensions
    wao_classes = dwd.pop("world_axis_object_classes")
    wao_components = dwd.pop("world_axis_object_components")
    validate_info_dict(
        dwd,
        {
            "value": [1.0e10, 5.0e00],
            "world_axis_physical_types": ["em.freq", "time"],
            "world_axis_names": ["Frequency", "Time"],
            "world_axis_units": ["Hz", "s"],
            "serialized_classes": False,
        },
    )
    assert wao_components[0][0:2] == ("spectral", 0)
    assert wao_components[1][0] == "time"
    assert wao_components[1][1] == 0

    assert wao_classes["spectral"][0:3] == (u.Quantity, (), {})
    assert wao_classes["time"][0:3] == (Time, (), {})


def test_pixel_to_world_values_different_int_types():
    int_sliced = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, np.s_[:, 0, :])
    np64_sliced = SlicedLowLevelWCS(WCS_SPECTRAL_CUBE, np.s_[:, np.int64(0), :])
    pixel_arrays = ([0, 1], [0, 1])
    for int_coord, np64_coord in zip(
        int_sliced.pixel_to_world_values(*pixel_arrays),
        np64_sliced.pixel_to_world_values(*pixel_arrays),
    ):
        assert all(int_coord == np64_coord)


COUPLED_WCS_HEADER = {
    "WCSAXES": 3,
    "CRPIX1": (100 + 1) / 2,
    "CRPIX2": (25 + 1) / 2,
    "CRPIX3": 1.0,
    "PC1_1": 0.0,
    "PC1_2": -1.0,
    "PC1_3": 0.0,
    "PC2_1": 1.0,
    "PC2_2": 0.0,
    "PC2_3": -1.0,
    "CDELT1": 5,
    "CDELT2": 5,
    "CDELT3": 0.055,
    "CUNIT1": "arcsec",
    "CUNIT2": "arcsec",
    "CUNIT3": "Angstrom",
    "CTYPE1": "HPLN-TAN",
    "CTYPE2": "HPLT-TAN",
    "CTYPE3": "WAVE",
    "CRVAL1": 0.0,
    "CRVAL2": 0.0,
    "CRVAL3": 1.05,
}


def test_coupled_world_slicing():
    fits_wcs = WCS(header=COUPLED_WCS_HEADER)
    sl = SlicedLowLevelWCS(fits_wcs, 0)
    world = fits_wcs.pixel_to_world_values(0, 0, 0)
    out_pix = sl.world_to_pixel_values(world[0], world[1])

    assert np.allclose(out_pix[0], 0)

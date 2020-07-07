import numpy as np
from numpy.testing import assert_allclose

import pytest
from pytest import raises

from astropy import units as u

from astropy.wcs import WCS
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs.wcsapi.utils import deserialize_class, wcs_info_str


@pytest.fixture
def axis_correlation_matrix():
    return _axis_correlation_matrix()


def _axis_correlation_matrix():
    shape = (4, 4)
    acm = np.zeros(shape, dtype=bool)
    for i in range(min(shape)):
        acm[i, i] = True
    acm[0, 1] = True
    acm[1, 0] = True
    acm[-1, 0] = True
    return acm


@pytest.fixture
def test_wcs():
    return TestWCS()


class TestWCS():
    def __init__(self):
        self.world_axis_physical_types = [
            'custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat', 'em.wl', 'time']
        self.axis_correlation_matrix = _axis_correlation_matrix()


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


def test_convert_between_array_and_pixel_axes():
    test_input = np.array([1, 4, -2])
    naxes = 5
    expected = np.array([3, 0, 1])
    output = utils.wcs.convert_between_array_and_pixel_axes(test_input, naxes)
    assert all(output == expected)


def test_pixel_axis_to_world_axes(axis_correlation_matrix):
    output = utils.wcs.pixel_axis_to_world_axes(0, axis_correlation_matrix)
    expected = np.array([0, 1, 3])
    assert all(output == expected)


def test_world_axis_to_pixel_axes(axis_correlation_matrix):
    output = utils.wcs.world_axis_to_pixel_axes(1, axis_correlation_matrix)
    expected = np.array([0, 1])
    assert all(output == expected)


def test_pixel_axis_to_physical_types(test_wcs):
    output = utils.wcs.pixel_axis_to_physical_types(0, test_wcs)
    expected = np.array(['custom:pos.helioprojective.lon',
                         'custom:pos.helioprojective.lat', 'time'])
    assert all(output == expected)


def test_physical_type_to_pixel_axes(test_wcs):
    output = utils.wcs.physical_type_to_pixel_axes('lon', test_wcs)
    expected = np.array([0, 1])
    assert all(output == expected)


@pytest.mark.parametrize("test_input,expected", [('wl', 2), ('em.wl', 2)])
def test_physical_type_to_world_axis(test_input, expected):
    world_axis_physical_types = ['custom:pos.helioprojective.lon',
                                 'custom:pos.helioprojective.lat', 'em.wl', 'time']
    output = utils.wcs.physical_type_to_world_axis(test_input, world_axis_physical_types)
    assert output == expected


def test_get_dependent_pixel_axes(axis_correlation_matrix):
    output = utils.wcs.get_dependent_pixel_axes(0, axis_correlation_matrix)
    expected = np.array([0, 1, 3])
    assert all(output == expected)


def test_get_dependent_array_axes(axis_correlation_matrix):
    output = utils.wcs.get_dependent_array_axes(3, axis_correlation_matrix)
    expected = np.array([0, 2, 3])
    assert all(output == expected)


def test_get_dependent_world_axes(axis_correlation_matrix):
    output = utils.wcs.get_dependent_world_axes(3, axis_correlation_matrix)
    expected = np.array([0, 3])
    assert all(output == expected)


def test_get_dependent_physical_types(test_wcs):
    output = utils.wcs.get_dependent_physical_types("time", test_wcs)
    expected = np.array(['custom:pos.helioprojective.lon', 'time'])
    assert all(output == expected)

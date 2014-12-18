# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from numpy.testing import assert_allclose

from ...tests.helper import pytest
from ..utils import extract_array, add_array, subpixel_indices

test_positions = [(10.52, 3.12), (5.62, 12.97), (31.33, 31.77),
                  (0.46, 0.94), (20.45, 12.12), (42.24, 24.42)]

test_position_indices = [(0, 3), (0, 2), (4, 1),
                         (4, 2), (4, 3), (3, 4)]

test_slices = [slice(10.52, 3.12), slice(5.62, 12.97),
               slice(31.33, 31.77), slice(0.46, 0.94),
               slice(20.45, 12.12), slice(42.24, 24.42)]

subsampling = 5


def test_extract_array():
    """
    Test extract_array utility function.

    Test by extracting an array of ones out of an array of zeros.
    """
    large_test_array = np.zeros((11, 11))
    small_test_array = np.ones((5, 5))
    large_test_array[3:8, 3:8] = small_test_array
    extracted_array = extract_array(large_test_array, (5, 5), (5, 5))
    assert np.all(extracted_array == small_test_array)


def test_add_array_odd_shape():
    """
    Test add_array utility function.

    Test by adding an array of ones out of an array of zeros.
    """
    large_test_array = np.zeros((11, 11))
    small_test_array = np.ones((5, 5))
    large_test_array_ref = large_test_array.copy()
    large_test_array_ref[3:8, 3:8] += small_test_array

    added_array = add_array(large_test_array, small_test_array, (5, 5))
    assert np.all(added_array == large_test_array_ref)


def test_add_array_even_shape():
    """
    Test add_array_2D utility function.

    Test by adding an array of ones out of an array of zeros.
    """
    large_test_array = np.zeros((11, 11))
    small_test_array = np.ones((4, 4))
    large_test_array_ref = large_test_array.copy()
    large_test_array_ref[0:2, 0:2] += small_test_array[2:4, 2:4]

    added_array = add_array(large_test_array, small_test_array, (0, 0))
    assert np.all(added_array == large_test_array_ref)


@pytest.mark.parametrize(('position', 'subpixel_index'),
                         zip(test_positions, test_position_indices))
def test_subpixel_indices(position, subpixel_index):
    """
    Test subpixel_indices utility function.

    Test by asserting that the function returns correct results for
    given test values.
    """
    assert np.all(subpixel_indices(position, subsampling) == subpixel_index)

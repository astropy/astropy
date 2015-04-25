# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from numpy.testing import assert_allclose

from ...tests.helper import pytest
from ..utils import (extract_array, add_array, subpixel_indices,
                     block_reduce, block_replicate)


try:
    import skimage
    HAS_SKIMAGE = True
except ImportError:
    HAS_SKIMAGE = False


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


@pytest.mark.skipif('not HAS_SKIMAGE')
class TestBlockReduce(object):
    def test_1d(self):
        """Test 1D array."""
        data = np.arange(4)
        expected = np.array([1, 5])
        result = block_reduce(data, 2)
        assert np.all(result == expected)

    def test_1d_mean(self):
        """Test 1D array with func=np.mean."""
        data = np.arange(4)
        block_size = 2.
        expected = block_reduce(data, block_size, func=np.sum) / block_size
        result_mean = block_reduce(data, block_size, func=np.mean)
        assert np.all(result_mean == expected)

    def test_2d(self):
        """Test 2D array."""
        data = np.arange(4).reshape(2, 2)
        expected = np.array([[6]])
        result = block_reduce(data, 2)
        assert np.all(result == expected)

    def test_2d_mean(self):
        """Test 2D array with func=np.mean."""
        data = np.arange(4).reshape(2, 2)
        block_size = 2.
        expected = (block_reduce(data, block_size, func=np.sum) /
                    block_size**2)
        result = block_reduce(data, block_size, func=np.mean)
        assert np.all(result == expected)

    def test_2d_trim(self):
        """
        Test trimming of 2D array when size is not perfectly divisible
        by block_size.
        """

        data1 = np.arange(15).reshape(5, 3)
        result1 = block_reduce(data1, 2)
        data2 = data1[0:4, 0:2]
        result2 = block_reduce(data2, 2)
        assert np.all(result1 == result2)

    def test_block_size_broadcasting(self):
        """Test scalar block_size broadcasting."""
        data = np.arange(16).reshape(4, 4)
        result1 = block_reduce(data, 2)
        result2 = block_reduce(data, (2, 2))
        assert np.all(result1 == result2)

    def test_block_size_len(self):
        """Test block_size length."""
        data = np.ones((2, 2))
        with pytest.raises(ValueError):
            block_reduce(data, (2, 2, 2))


@pytest.mark.skipif('not HAS_SKIMAGE')
class TestBlockReplicate(object):
    def test_1d(self):
        """Test 1D array."""
        data = np.arange(2)
        expected = np.array([0, 0, 0.5, 0.5])
        result = block_replicate(data, 2)
        assert np.all(result == expected)

    def test_1d_conserve_sum(self):
        """Test 1D array with conserve_sum=False."""
        data = np.arange(2)
        block_size = 2.
        expected = block_replicate(data, block_size) * block_size
        result = block_replicate(data, block_size, conserve_sum=False)
        assert np.all(result == expected)

    def test_2d(self):
        """Test 2D array."""
        data = np.arange(2).reshape(2, 1)
        expected = np.array([[0, 0], [0, 0], [0.25, 0.25], [0.25, 0.25]])
        result = block_replicate(data, 2)
        assert np.all(result == expected)

    def test_2d_conserve_sum(self):
        """Test 2D array with conserve_sum=False."""
        data = np.arange(6).reshape(2, 3)
        block_size = 2.
        expected = block_replicate(data, block_size) * block_size**2
        result = block_replicate(data, block_size, conserve_sum=False)
        assert np.all(result == expected)

    def test_block_size_broadcasting(self):
        """Test scalar block_size broadcasting."""
        data = np.arange(4).reshape(2, 2)
        result1 = block_replicate(data, 2)
        result2 = block_replicate(data, (2, 2))
        assert np.all(result1 == result2)

    def test_block_size_len(self):
        """Test block_size length."""
        data = np.arange(5)
        with pytest.raises(ValueError):
            block_replicate(data, (2, 2))

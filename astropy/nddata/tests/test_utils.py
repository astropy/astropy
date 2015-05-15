# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np

from ...tests.helper import pytest
from ..utils import (extract_array, add_array, subpixel_indices,
                     block_reduce, block_replicate,
                     overlap_slices, NoOverlapError, PartialOverlapError)

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

test_pos_bad = [(-1, -4), (-1, 0), (6, 2), (6, 6)]


def test_slices_different_dim():
    '''Overlap from arrays with different number of dim is undefined.'''
    with pytest.raises(ValueError) as e:
        temp = overlap_slices((4, 5, 6), (1, 2), (0, 0))
    assert "the same number of dimensions" in str(e.value)


def test_slices_pos_different_dim():
    '''Position must have same dim as arrays.'''
    with pytest.raises(ValueError) as e:
        temp = overlap_slices((4, 5), (1, 2), (0, 0, 3))
    assert "the same number of dimensions" in str(e.value)


@pytest.mark.parametrize('pos', test_pos_bad)
def test_slices_no_overlap(pos):
    '''If there is no overlap between arrays, an error should be raised.'''
    with pytest.raises(NoOverlapError):
        temp = overlap_slices((5, 5), (2, 2), pos)


def test_slices_partial_overlap():
    '''Compute a slice for partially overlapping arrays.'''
    temp = overlap_slices((5,), (3,), (0,))
    assert temp == ((slice(0, 2, None),), (slice(1, 3, None),))

    temp = overlap_slices((5,), (3,), (0,), mode='partial')
    assert temp == ((slice(0, 2, None),), (slice(1, 3, None),))

    for pos in [0, 4]:
        with pytest.raises(PartialOverlapError) as e:
            temp = overlap_slices((5,), (3,), (pos,), mode='strict')
        assert 'Arrays overlap only partially.' in str(e.value)


def test_slices_overlap_wrong_mode():
    '''Call overlap_slices with non-existing mode.'''
    with pytest.raises(ValueError) as e:
        temp = overlap_slices((5,), (3,), (0,), mode='full')
    assert "Mode can be only" in str(e.value)


def test_extract_array_wrong_mode():
    '''Call extract_array with non-existing mode.'''
    with pytest.raises(ValueError) as e:
        temp = extract_array(np.arange(4), (2, ), (0, ), mode='full')
    assert "Valid modes are 'partial', 'trim', and 'strict'." == str(e.value)


def test_extract_array_1d_even():
    '''Extract 1 d arrays.

    All dimensions are treated the same, so we can test in 1 dim.
    '''
    assert np.all(extract_array(np.arange(4), (2, ), (0, ), fill_value=-99) == np.array([-99, 0]))
    for i in [1, 2, 3]:
        assert np.all(extract_array(np.arange(4), (2, ), (i, )) == np.array([i -1 , i]))
    assert np.all(extract_array(np.arange(4.), (2, ), (4, ), fill_value=np.inf) == np.array([3, np.inf]))


def test_extract_array_1d_odd():
    '''Extract 1 d arrays.

    All dimensions are treated the same, so we can test in 1 dim.
    The first few lines test the most error-prone part: Extraction of an
    array on the boundaries.
    Additional tests (e.g. dtype of return array) are done for the last
    case only.
    '''
    assert np.all(extract_array(np.arange(4), (3,), (-1, ), fill_value=-99) == np.array([-99, -99, 0]))
    assert np.all(extract_array(np.arange(4), (3,), (0, ), fill_value=-99) == np.array([-99, 0, 1]))
    for i in [1,2]:
        assert np.all(extract_array(np.arange(4), (3,), (i, )) == np.array([i-1, i, i+1]))
    assert np.all(extract_array(np.arange(4), (3,), (3, ), fill_value=-99) == np.array([2, 3, -99]))
    arrayin = np.arange(4.)
    extracted = extract_array(arrayin, (3,), (4, ))
    assert extracted[0] == 3
    assert np.isnan(extracted[1]) # since I cannot use `==` to test for nan
    assert extracted.dtype == arrayin.dtype


def test_extract_array_1d():
    """In 1d, shape can be int instead of tuple"""
    assert np.all(extract_array(np.arange(4), 3, (-1, ), fill_value=-99) == np.array([-99, -99, 0]))
    assert np.all(extract_array(np.arange(4), 3, -1, fill_value=-99) == np.array([-99, -99, 0]))


def test_extract_Array_float():
    """integer is at bin center"""
    for a in np.arange(2.51, 3.49, 0.1):
        assert np.all(extract_array(np.arange(5), 3, a) == np.array([2, 3, 4]))


def test_extract_array_1d_trim():
    '''Extract 1 d arrays.

    All dimensions are treated the same, so we can test in 1 dim.
    '''
    assert np.all(extract_array(np.arange(4), (2, ), (0, ), mode='trim') == np.array([0]))
    for i in [1, 2, 3]:
        assert np.all(extract_array(np.arange(4), (2, ), (i, ), mode='trim') == np.array([i -1 , i]))
    assert np.all(extract_array(np.arange(4.), (2, ), (4, ), mode='trim') == np.array([3]))


@pytest.mark.parametrize('mode', ['partial', 'trim', 'strict'])
def test_extract_array_easy(mode):
    """
    Test extract_array utility function.

    Test by extracting an array of ones out of an array of zeros.
    """
    large_test_array = np.zeros((11, 11))
    small_test_array = np.ones((5, 5))
    large_test_array[3:8, 3:8] = small_test_array
    extracted_array = extract_array(large_test_array, (5, 5), (5, 5), mode=mode)
    assert np.all(extracted_array == small_test_array)


def test_extract_array_return_pos():
    '''Check that the return position is calculated correctly.

    The result will differ by mode. All test here are done in 1d because it's
    easier to construct correct test cases.
    '''
    large_test_array = np.arange(5)
    for i in np.arange(-1, 6):
        extracted, new_pos = extract_array(large_test_array, 3, i,
                                           mode='partial', return_position=True)
        assert new_pos == (1, )
    # Now check an array with an even number
    for i, expected in zip([1.49, 1.51, 3], [1.49, 0.51, 1]):
        extracted, new_pos = extract_array(large_test_array, (2,), (i,),
                                           mode='strict', return_position=True)
        assert new_pos == (expected, )
    # For mode='trim' the answer actually depends
    for i, expected in zip(np.arange(-1, 6), (-1, 0, 1, 1, 1, 1, 1)):
        extracted, new_pos = extract_array(large_test_array, (3,), (i,),
                                           mode='trim', return_position=True)
        assert new_pos == (expected, )


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

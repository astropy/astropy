# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np

from ...tests.helper import pytest
from ..utils import extract_array, add_array, subpixel_indices, \
                    overlap_slices, NoOverlapError, PartialOverlapError

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
    assert "Mode can only be" in str(e.value)


def test_extract_array_wrong_mode():
    '''Call extract_array with non-existing mode.'''
    with pytest.raises(ValueError) as e:
        temp = extract_array(np.arange(4), (2, ), (0, ), mode='full')
    assert "Valid modes are 'partial', 'trim', and 'strict'" == str(e.value)


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

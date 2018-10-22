# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from ...tests.helper import catch_warnings
from ..bitmask import bitfield_to_boolean_mask


def test_bitfield_not_integer():
    with pytest.raises(TypeError):
        bitfield_to_boolean_mask(np.random.random((10, 10)))


def test_bitfield_negative_flags():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(ValueError):
        bitfield_to_boolean_mask(bm, [-1])


def test_bitfield_non_poweroftwo_flags():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(ValueError):
        bitfield_to_boolean_mask(bm, [3])


def test_bitfield_flipbits_when_no_bits():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(TypeError):
        bitfield_to_boolean_mask(bm, None, flip_bits=1)


def test_bitfield_flipbits_when_stringbits():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(TypeError):
        bitfield_to_boolean_mask(bm, '3', flip_bits=1)


def test_bitfield_string_flag_flip_not_start_of_string():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(ValueError):
        bitfield_to_boolean_mask(bm, '1, ~4')


def test_bitfield_string_flag_unbalanced_parens():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(ValueError):
        bitfield_to_boolean_mask(bm, '(1, 4))')


def test_bitfield_string_flag_wrong_positioned_parens():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(ValueError):
        bitfield_to_boolean_mask(bm, '((1, )4)')


def test_bitfield_string_flag_empty():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(ValueError):
        bitfield_to_boolean_mask(bm, '~')


def test_bitfield_flag_non_integer():
    bm = np.random.randint(0, 10, (10, 10))
    with pytest.raises(TypeError):
        bitfield_to_boolean_mask(bm, [1.3])


def test_bitfield_duplicate_flag_throws_warning():
    bm = np.random.randint(0, 10, (10, 10))
    with catch_warnings(UserWarning) as w:
        bitfield_to_boolean_mask(bm, [1, 1])
    assert len(w)


def test_bitfield_none_identical_to_strNone():
    bm = np.random.randint(0, 10, (10, 10))
    m1 = bitfield_to_boolean_mask(bm, None)
    m2 = bitfield_to_boolean_mask(bm, 'None')
    np.testing.assert_array_equal(m1, m2)


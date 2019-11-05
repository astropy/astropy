"""
A module containing unit tests for the `bitmask` modue.

Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
import warnings
import numpy as np
import pytest

from astropy.nddata import bitmask


MAX_INT_TYPE = np.maximum_sctype(np.int_)
MAX_UINT_TYPE = np.maximum_sctype(np.uint)
MAX_UINT_FLAG = np.left_shift(
    MAX_UINT_TYPE(1),
    MAX_UINT_TYPE(np.iinfo(MAX_UINT_TYPE).bits - 1)
)
MAX_INT_FLAG = np.left_shift(
    MAX_INT_TYPE(1),
    MAX_INT_TYPE(np.iinfo(MAX_INT_TYPE).bits - 2)
)
SUPER_LARGE_FLAG = 1 << np.iinfo(MAX_UINT_TYPE).bits
EXTREME_TEST_DATA = np.array([
        0, 1, 1 + 1 << 2, MAX_INT_FLAG, ~0, MAX_INT_TYPE(MAX_UINT_FLAG),
        1 + MAX_INT_TYPE(MAX_UINT_FLAG)
], dtype=MAX_INT_TYPE)


@pytest.mark.parametrize('flag', [0, -1])
def test_nonpositive_not_a_bit_flag(flag):
    assert not bitmask._is_bit_flag(n=flag)


@pytest.mark.parametrize('flag', [
    1, MAX_UINT_FLAG, int(MAX_UINT_FLAG), SUPER_LARGE_FLAG
])
def test_is_bit_flag(flag):
    assert bitmask._is_bit_flag(n=flag)


@pytest.mark.parametrize('number', [0, 1, MAX_UINT_FLAG, SUPER_LARGE_FLAG])
def test_is_int(number):
    assert bitmask._is_int(number)


@pytest.mark.parametrize('number', ['1', True, 1.0])
def test_nonint_is_not_an_int(number):
    assert not bitmask._is_int(number)


@pytest.mark.parametrize('flag,flip,expected', [
    (3, None, 3),
    (3, True, -4),
    (3, False, 3),
    ([1, 2], False, 3),
    ([1, 2], True, -4)
])
def test_interpret_valid_int_bit_flags(flag, flip, expected):
    assert(
        bitmask.interpret_bit_flags(bit_flags=flag, flip_bits=flip) == expected
    )


@pytest.mark.parametrize('flag', [None, ' ', 'None', 'Indef'])
def test_interpret_none_bit_flags_as_None(flag):
    assert bitmask.interpret_bit_flags(bit_flags=flag) is None


@pytest.mark.parametrize('flag,expected', [
    ('1', 1),
    ('~-1', ~(-1)),
    ('~1', ~1),
    ('1,2', 3),
    ('1+2', 3),
    ('(1,2)', 3),
    ('(1+2)', 3),
    ('~1,2', ~3),
    ('~1+2', ~3),
    ('~(1,2)', ~3),
    ('~(1+2)', ~3)
])
def test_interpret_valid_str_bit_flags(flag, expected):
    assert(
        bitmask.interpret_bit_flags(bit_flags=flag) == expected
    )


@pytest.mark.parametrize('flag,flip', [
    (None, True),
    (' ', True),
    ('None', True),
    ('Indef', True),
    (None, False),
    (' ', False),
    ('None', False),
    ('Indef', False),
    ('1', True),
    ('1', False)
])
def test_interpret_None_or_str_and_flip_incompatibility(flag, flip):
    with pytest.raises(TypeError):
        bitmask.interpret_bit_flags(bit_flags=flag, flip_bits=flip)


@pytest.mark.parametrize('flag', [True, 1.0, [1.0], object])
def test_interpret_wrong_flag_type(flag):
    with pytest.raises(TypeError):
        bitmask.interpret_bit_flags(bit_flags=flag)


@pytest.mark.parametrize('flag', ['SOMETHING', '1.0,2,3'])
def test_interpret_wrong_string_int_format(flag):
    with pytest.raises(ValueError):
        bitmask.interpret_bit_flags(bit_flags=flag)


def test_interpret_duplicate_flag_warning():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        assert bitmask.interpret_bit_flags([2, 4, 4]) == 6
        assert len(w)
        assert issubclass(w[-1].category, UserWarning)
        assert "Duplicate" in str(w[-1].message)


@pytest.mark.parametrize('flag', [[1, 2, 3], '1, 2, 3'])
def test_interpret_non_flag(flag):
    with pytest.raises(ValueError):
        bitmask.interpret_bit_flags(bit_flags=flag)


def test_interpret_allow_single_value_str_nonflags():
    assert bitmask.interpret_bit_flags(bit_flags=str(3)) == 3


@pytest.mark.parametrize('flag', [
    '~',
    '( )',
    '(~1,2)',
    '~(1,2',
    '1,~2',
    '1,(2,4)',
    '1,2+4',
    '1+4,2'
])
def test_interpret_bad_str_syntax(flag):
    with pytest.raises(ValueError):
        bitmask.interpret_bit_flags(bit_flags=flag)


def test_bitfield_must_be_integer_check():
    with pytest.raises(TypeError):
        bitmask.bitfield_to_boolean_mask(1.0, 1)


@pytest.mark.parametrize('data,flags,flip,goodval,dtype,ref', [
    (EXTREME_TEST_DATA, None, None, True, np.bool_,
     EXTREME_TEST_DATA.size * [1]),
    (EXTREME_TEST_DATA, None, None, False, np.bool_,
     EXTREME_TEST_DATA.size * [0]),
    (EXTREME_TEST_DATA, [1, MAX_UINT_FLAG], False, True, np.bool_,
     [1, 1, 0, 0, 0, 1, 1]),
    (EXTREME_TEST_DATA, None, None, True, np.bool_,
     EXTREME_TEST_DATA.size * [1]),
    (EXTREME_TEST_DATA, [1, MAX_UINT_FLAG], False, False, np.bool_,
     [0, 0, 1, 1, 1, 0, 0]),
    (EXTREME_TEST_DATA, [1, MAX_UINT_FLAG], True, True, np.int8,
     [1, 0, 1, 1, 0, 0, 0])
])
def test_bitfield_to_boolean_mask(data, flags, flip, goodval, dtype, ref):
    mask = bitmask.bitfield_to_boolean_mask(
        bitfield=data,
        ignore_flags=flags,
        flip_bits=flip,
        good_mask_value=goodval,
        dtype=dtype
    )

    assert(mask.dtype == dtype)
    assert np.all(mask == ref)

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
    ('1|2', 3),
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


@pytest.mark.parametrize('flag,expected', [
    ('CR', 1),
    ('~CR', ~1),
    ('CR|HOT', 3),
    ('CR,HOT', 3),
    ('CR+HOT', 3),
    (['CR', 'HOT'], 3),
    ('(CR,HOT)', 3),
    ('(HOT+CR)', 3),
    ('~HOT,CR', ~3),
    ('~CR+HOT', ~3),
    ('~(HOT,CR)', ~3),
    ('~(HOT|CR)', ~3),
    ('~(CR+HOT)', ~3)
])
def test_interpret_valid_mnemonic_bit_flags(flag, expected):
    flagmap = bitmask.extend_bit_flag_map('DetectorMap', CR=1, HOT=2)

    assert(
        bitmask.interpret_bit_flags(bit_flags=flag, flag_name_map=flagmap)
        == expected
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
    '1+4,2',
    '1|4+2'
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


@pytest.mark.parametrize('flag', [(4, 'flag1'), 8])
def test_bitflag(flag):
    f = bitmask.BitFlag(flag)
    if isinstance(flag, tuple):
        assert f == flag[0]
        assert f.__doc__ == flag[1]

        f = bitmask.BitFlag(*flag)
        assert f == flag[0]
        assert f.__doc__ == flag[1]

    else:
        assert f == flag


def test_bitflag_docs2():
    with pytest.raises(ValueError):
        bitmask.BitFlag((1, 'docs1'), 'docs2')


@pytest.mark.parametrize('flag', [0, 3])
def test_bitflag_not_pow2(flag):
    with pytest.raises(bitmask.InvalidBitFlag):
        bitmask.BitFlag(flag, 'custom flag')


@pytest.mark.parametrize('flag', [0.0, True, '1'])
def test_bitflag_not_int_flag(flag):
    with pytest.raises(bitmask.InvalidBitFlag):
        bitmask.BitFlag((flag, 'custom flag'))


@pytest.mark.parametrize('caching', [True, False])
def test_basic_map(monkeypatch, caching):
    monkeypatch.setattr(bitmask, '_ENABLE_BITFLAG_CACHING', False)

    class ObservatoryDQMap(bitmask.BitFlagNameMap):
        _not_a_flag = 1
        CR = 1, 'cosmic ray'
        HOT = 2
        DEAD = 4

    class DetectorMap(ObservatoryDQMap):
        __version__ = '1.0'
        _not_a_flag = 181
        READOUT_ERR = 16

    assert ObservatoryDQMap.cr == 1
    assert ObservatoryDQMap.cr.__doc__ == 'cosmic ray'
    assert DetectorMap.READOUT_ERR == 16


@pytest.mark.parametrize('caching', [True, False])
def test_extend_map(monkeypatch, caching):
    monkeypatch.setattr(bitmask, '_ENABLE_BITFLAG_CACHING', caching)

    class ObservatoryDQMap(bitmask.BitFlagNameMap):
        CR = 1
        HOT = 2
        DEAD = 4

    DetectorMap = bitmask.extend_bit_flag_map(
        'DetectorMap', ObservatoryDQMap,
        __version__='1.0',
        DEAD=4,
        READOUT_ERR=16
    )

    assert DetectorMap.CR == 1
    assert DetectorMap.readout_err == 16


@pytest.mark.parametrize('caching', [True, False])
def test_extend_map_redefine_flag(monkeypatch, caching):
    monkeypatch.setattr(bitmask, '_ENABLE_BITFLAG_CACHING', caching)

    class ObservatoryDQMap(bitmask.BitFlagNameMap):
        CR = 1
        HOT = 2
        DEAD = 4

    with pytest.raises(AttributeError):
        bitmask.extend_bit_flag_map(
            'DetectorMap',
            ObservatoryDQMap,
            __version__='1.0',
            DEAD=32
        )

    with pytest.raises(AttributeError):
        bitmask.extend_bit_flag_map(
            'DetectorMap',
            ObservatoryDQMap,
            __version__='1.0',
            DEAD=32,
            dead=64
        )


@pytest.mark.parametrize('caching', [True, False])
def test_map_redefine_flag(monkeypatch, caching):
    monkeypatch.setattr(bitmask, '_ENABLE_BITFLAG_CACHING', caching)

    class ObservatoryDQMap(bitmask.BitFlagNameMap):
        _not_a_flag = 8
        CR = 1
        HOT = 2
        DEAD = 4

    with pytest.raises(AttributeError):
        class DetectorMap1(ObservatoryDQMap):
            __version__ = '1.0'
            CR = 16

    with pytest.raises(AttributeError):
        class DetectorMap2(ObservatoryDQMap):
            SHADE = 8
            _FROZEN = 16

        DetectorMap2.novel = 32

    with pytest.raises(AttributeError):
        bitmask.extend_bit_flag_map(
            'DetectorMap', ObservatoryDQMap,
            READOUT_ERR=16,
            SHADE=32,
            readout_err=128
        )


def test_map_cant_modify_version():
    class ObservatoryDQMap(bitmask.BitFlagNameMap):
        __version__ = '1.2.3'
        CR = 1

    assert ObservatoryDQMap.__version__ == '1.2.3'
    assert ObservatoryDQMap.CR == 1

    with pytest.raises(AttributeError):
        ObservatoryDQMap.__version__ = '3.2.1'


@pytest.mark.parametrize('flag', [0, 3])
def test_map_not_bit_flag(flag):
    with pytest.raises(ValueError):
        bitmask.extend_bit_flag_map('DetectorMap', DEAD=flag)

    with pytest.raises(ValueError):
        class DetectorMap(bitmask.BitFlagNameMap):
            DEAD=flag


@pytest.mark.parametrize('flag', [0.0, True, '1'])
def test_map_not_int_flag(flag):
    with pytest.raises(bitmask.InvalidBitFlag):
        bitmask.extend_bit_flag_map('DetectorMap', DEAD=flag)

    with pytest.raises(bitmask.InvalidBitFlag):
        class ObservatoryDQMap(bitmask.BitFlagNameMap):
            CR = flag


def test_map_access_undefined_flag():
    DetectorMap = bitmask.extend_bit_flag_map('DetectorMap', DEAD=1)

    with pytest.raises(AttributeError):
        DetectorMap.DEAD1

    with pytest.raises(AttributeError):
        DetectorMap['DEAD1']


def test_map_delete_flag():
    DetectorMap = bitmask.extend_bit_flag_map('DetectorMap', DEAD=1)

    with pytest.raises(AttributeError):
        del DetectorMap.DEAD1

    with pytest.raises(AttributeError):
        del DetectorMap['DEAD1']


def test_map_repr():
    DetectorMap = bitmask.extend_bit_flag_map('DetectorMap', DEAD=1)
    assert repr(DetectorMap) == "<BitFlagNameMap 'DetectorMap'>"


def test_map_add_flags():
    map1 = bitmask.extend_bit_flag_map('DetectorMap', CR=1)

    map2 = map1 + {'HOT': 2, 'DEAD': (4, 'a really dead pixel')}
    assert map2.CR == 1
    assert map2.HOT == 2
    assert map2.DEAD.__doc__ == 'a really dead pixel'
    assert map2.DEAD == 4

    map2 = map1 + [('HOT', 2), ('DEAD', 4)]
    assert map2.CR == 1
    assert map2.HOT == 2

    map2 = map1 + ('HOT', 2)
    assert map2.CR == 1
    assert map2.HOT == 2

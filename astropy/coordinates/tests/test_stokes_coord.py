import numpy as np
import pytest
from numpy.testing import assert_equal

from astropy import units as u
from astropy.coordinates.stokes_coord import (
    StokesCoord,
    StokesSymbol,
    custom_stokes_symbol_mapping,
)
from astropy.utils import unbroadcast


def test_scalar():
    sk = StokesCoord(2)
    assert repr(sk) == "<StokesCoord 'Q'>"
    assert sk.value == 2.0
    assert sk.symbol == "Q"


def test_vector():
    # This also checks that floats are rounded when converting
    # to strings
    values = [1.2, 1.8, 2.0, 2.2, 2.8]
    sk = StokesCoord(values)
    assert_equal(sk.value, values)
    assert_equal(sk.symbol, np.array(["I", "Q", "Q", "Q", "U"]))
    assert repr(sk) == "<StokesCoord ['I', 'Q', 'Q', 'Q', 'U']>"


def test_vector_list_init():
    sk = StokesCoord(["I", "Q", "Q", "U", "U"])
    assert repr(sk) == "<StokesCoord ['I', 'Q', 'Q', 'U', 'U']>"
    assert_equal(sk.symbol, np.array(["I", "Q", "Q", "U", "U"]))


def test_unit():
    StokesCoord(1, unit=u.one)
    for unit in [u.radian, u.deg, u.Hz]:
        with pytest.raises(
            u.UnitsError, match="unit should not be specified explicitly"
        ):
            StokesCoord(1, unit=unit)


def test_undefined():
    sk = StokesCoord(np.arange(-10, 7))
    assert_equal(
        sk.symbol,
        np.array(
            [
                "?",
                "?",
                "YX",
                "XY",
                "YY",
                "XX",
                "LR",
                "RL",
                "LL",
                "RR",
                "?",
                "I",
                "Q",
                "U",
                "V",
                "?",
                "?",
            ]
        ),
    )


def test_undefined_init():
    with pytest.raises(Exception, match="Unknown stokes symbols.*Spam"):
        StokesCoord("Spam")


def test_custom_symbol_mapping():
    custom_mapping = {
        10000: StokesSymbol("A"),
        10001: StokesSymbol("B"),
        10002: StokesSymbol("C"),
        10003: StokesSymbol("D"),
    }

    # Check that we can supply a custom mapping
    with custom_stokes_symbol_mapping(custom_mapping):
        values = [0.6, 1.7, 10000.1, 10002.4]
        sk1 = StokesCoord(values)
        assert repr(sk1) == "<StokesCoord ['I', 'Q', 'A', 'C']>"
        assert_equal(sk1.value, values)
        assert_equal(sk1.symbol, np.array(["I", "Q", "A", "C"]))

    # Check that the mapping is not active outside the context manager
    assert_equal(sk1.symbol, np.array(["I", "Q", "?", "?"]))

    # But not for new StokesCoords
    sk2 = StokesCoord(values)
    assert_equal(sk2.symbol, np.array(["I", "Q", "?", "?"]))


def test_custom_symbol_mapping_overlap():
    # Make a custom mapping that overlaps with some of the existing values

    custom_mapping = {
        3: StokesSymbol("A"),
        4: StokesSymbol("B"),
        5: StokesSymbol("C"),
        6: StokesSymbol("D"),
    }

    with custom_stokes_symbol_mapping(custom_mapping):
        sk = StokesCoord(np.arange(1, 7))
        assert_equal(sk.symbol, np.array(["I", "Q", "A", "B", "C", "D"]))


def test_custom_symbol_mapping_replace():
    # Check that we can replace the mapping completely

    custom_mapping = {
        3: StokesSymbol("A"),
        4: StokesSymbol("B"),
        5: StokesSymbol("C"),
        6: StokesSymbol("D"),
    }

    with custom_stokes_symbol_mapping(custom_mapping, replace=True):
        sk = StokesCoord(np.arange(1, 7))
        assert_equal(sk.symbol, np.array(["?", "?", "A", "B", "C", "D"]))


def test_comparison_scalar():
    sk = StokesCoord(np.arange(1, 6))
    assert_equal(np.equal(sk, "Q"), [False, True, False, False, False])
    assert_equal(np.equal("Q", sk), [False, True, False, False, False])
    assert_equal("Q" == sk, [False, True, False, False, False])
    assert_equal(sk == 1, [True, False, False, False, False])
    assert_equal(sk == "Q", [False, True, False, False, False])
    assert_equal(sk == "?", [False, False, False, False, True])


def test_comparison_vector():
    sk = StokesCoord(np.arange(1, 6))
    assert_equal(
        sk == np.array(["I", "Q", "I", "I", "Q"]), [True, True, False, False, False]
    )


def test_comparison_other_coord():
    sk1 = StokesCoord(np.arange(1, 6))
    sk2 = StokesCoord("I")
    assert_equal(sk1 == sk2, [True, False, False, False, False])
    sk3 = StokesCoord(np.repeat(2, 5))
    assert_equal(sk1 == sk3, [False, True, False, False, False])


def test_efficient():
    # Make sure that if we pass a broadcasted array in we get a broadcasted
    # arrayt of symbols.

    values = np.broadcast_to(np.arange(1, 5), (512, 256, 4))
    sk = StokesCoord(values)
    assert sk.symbol.shape == (512, 256, 4)
    assert unbroadcast(sk.value).shape == (1, 1, 4)
    assert unbroadcast(sk.symbol).shape == (1, 1, 4)
    assert_equal(unbroadcast(sk.symbol)[0, 0], np.array(["I", "Q", "U", "V"]))

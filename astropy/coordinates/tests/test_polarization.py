import numpy as np
import pytest
from numpy.testing import assert_equal

from astropy.coordinates.polarization import (
    StokesCoord,
    StokesSymbol,
    custom_stokes_symbol_mapping,
)
from astropy.table import Table, vstack
from astropy.utils import unbroadcast


def test_scalar():
    sk = StokesCoord(2)
    assert repr(sk) == "StokesCoord('Q')" == str(sk)
    assert sk.value == 2.0
    assert sk.symbol == "Q"


def test_vector():
    # This also checks that floats are rounded when converting
    # to strings
    values = [1.2, 1.8, 2.0, 2.2, 2.8]
    sk = StokesCoord(values)
    assert_equal(sk.value, values)
    assert_equal(sk.symbol, np.array(["I", "Q", "Q", "Q", "U"]))
    assert repr(sk) == "StokesCoord(['I', 'Q', 'Q', 'Q', 'U'])" == str(sk)


def test_vector_list_init():
    sk = StokesCoord(["I", "Q", "Q", "U", "U"])
    assert repr(sk) == "StokesCoord(['I', 'Q', 'Q', 'U', 'U'])" == str(sk)
    assert_equal(sk.symbol, np.array(["I", "Q", "Q", "U", "U"]))


def test_undefined():
    sk = StokesCoord(np.arange(-10, 7))
    assert_equal(sk.symbol, "? ? YX XY YY XX LR RL LL RR ? I Q U V ? ?".split())  # noqa: SIM905


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
        assert repr(sk1) == "StokesCoord(['I', 'Q', 'A', 'C'])" == str(sk1)
        assert_equal(sk1.value, values)
        assert_equal(sk1.symbol, np.array(["I", "Q", "A", "C"]))

    # Check that the mapping is not active outside the context manager
    assert_equal(sk1.symbol, np.array(["I", "Q", "?", "?"]))

    # Also not for new StokesCoords
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
    assert_equal("Q" == sk, [False, True, False, False, False])
    assert_equal(sk == 1, [True, False, False, False, False])
    assert_equal(sk == "Q", [False, True, False, False, False])


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
    # array of symbols.

    values = np.broadcast_to(np.arange(1, 5, dtype=float), (512, 256, 4))
    sk = StokesCoord(values, copy=False)

    assert sk.symbol.shape == (512, 256, 4)
    assert unbroadcast(sk.value).shape == (4,)
    assert unbroadcast(sk.symbol).shape == (4,)
    assert_equal(unbroadcast(sk.symbol), np.array(["I", "Q", "U", "V"]))


def test_broadcast_to():
    sk = StokesCoord(np.arange(1, 5, dtype=int), copy=False)
    sk2 = np.broadcast_to(sk, (512, 256, 4))

    assert sk2.symbol.shape == (512, 256, 4)
    assert unbroadcast(sk2.value).shape == (4,)
    assert unbroadcast(sk2.symbol).shape == (4,)
    assert_equal(unbroadcast(sk.symbol), np.array(["I", "Q", "U", "V"]))


def test_table_vstack_stokes():
    sk = StokesCoord(np.arange(1, 5, dtype=int), copy=False)
    tt = Table([sk])
    assert isinstance(tt["col0"], StokesCoord)
    assert np.allclose(tt["col0"].value, np.arange(1, 5, dtype=int))

    sk2 = StokesCoord([1, 2, 2, 2, 4, 5])
    tt2 = Table([sk2])
    assert isinstance(tt2["col0"], StokesCoord)
    assert np.allclose(tt2["col0"].value, sk2.value)

    tt3 = vstack([tt, tt2])
    assert isinstance(tt3["col0"], StokesCoord)
    assert len(tt3) == 10
    assert np.allclose(tt3["col0"].value, np.array([1, 2, 3, 4, 1, 2, 2, 2, 4, 5]))


def test_init_copy():
    input = np.arange(1, 5, dtype=int)
    sk1 = StokesCoord(input, copy=False)
    assert sk1._data is input

    skc = StokesCoord(input, copy=True)
    assert skc._data is not input

    sk2 = StokesCoord(sk1)
    assert sk1._data is sk2._data

    sk3 = StokesCoord(sk1, copy=True)
    assert sk1._data is not sk3._data
    assert np.allclose(sk1._data, sk3._data)


def test_init_error():
    with pytest.raises(ValueError, match="object array"):
        StokesCoord(None)

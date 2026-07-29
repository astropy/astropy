import operator
import pickle

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from scaler import Scaler


@pytest.mark.parametrize(
    "dtype",
    [
        float,
        np.float32,
        np.float64,
        complex,
        np.complex64,
        np.complex128,
        int,
        np.int32,
        np.uint64,
    ],
)
def test_with_scalar(dtype):
    a = dtype(5)
    if issubclass(dtype, (complex, np.complexfloating)):
        a = a - 4.0j

    assert isinstance(a, dtype)
    assert Scaler(10.0)(a) == a * 10.0


@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("swapbyteorder", [False, True])
@pytest.mark.parametrize(
    "dtype",
    [
        np.float32,
        np.float64,
        np.complex64,
        np.complex128,
        np.int32,
        np.uint64,
    ],
)
def test_with_array(dtype, swapbyteorder, order):
    dtype = np.dtype(dtype)
    if swapbyteorder:
        dtype = dtype.newbyteorder()
    a = np.arange(10, dtype=dtype).reshape(2, 5)
    if dtype.kind == "c":
        a -= 4.0j
    a = a.copy(order=order)

    assert a.dtype == dtype
    if order == "C":
        assert a.flags["C_CONTIGUOUS"] and not a.flags["F_CONTIGUOUS"]
    else:
        assert a.flags["F_CONTIGUOUS"] and not a.flags["C_CONTIGUOUS"]
    got = Scaler(10.0)(a)
    exp = a * 10.0
    assert_array_equal(got, exp, strict=True)


def test_with_list():
    a = [[1.0, 2.0], [3.0, 4.0]]
    got = Scaler(10.0)(a)
    exp = np.array(a) * 10.0
    assert got.dtype == exp.dtype
    assert_array_equal(got, exp)


@pytest.mark.parametrize("item", [["abc"], [{}]])
def test_with_wrong_items(item):
    sc = Scaler(10.0)
    with pytest.raises(ValueError, match="compatible.*convertible"):
        sc(item)


def test_repr():
    sc = Scaler(10.0)
    assert repr(sc) == "Scaler(10.0)"


def test_pickle():
    sc = Scaler(10.0)
    p = pickle.dumps(sc)
    got = pickle.loads(p)
    assert type(got) is Scaler
    assert got.factor == sc.factor


def test_equality():
    sc1 = Scaler(10.0)
    sc1_2 = Scaler(10.0)
    assert sc1_2 == sc1
    assert not (sc1_2 != sc1)  # noqa: SIM202
    sc2 = Scaler(20.0)
    assert sc1 != sc2
    assert not (sc1 == sc2)


def test_hash():
    sc1 = Scaler(10.0)
    assert hash(sc1) != hash(10.0)
    sc1_2 = Scaler(10.0)
    assert hash(sc1_2) == hash(sc1)
    sc2 = Scaler(20.0)
    assert hash(sc2) != hash(sc1)


@pytest.mark.parametrize("op", [operator.lt, operator.le, operator.gt, operator.ge])
def test_comparison_impossible(op):
    sc = Scaler(10.0)
    with pytest.raises(TypeError, match="not supported"):
        op(sc, sc)

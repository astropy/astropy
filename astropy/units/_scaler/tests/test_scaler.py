import operator
import pickle

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from scaler import Scaler


NUMPY_SCALAR_TYPES = [np.float32, np.float64, np.complex64, np.complex128]
SCALAR_TYPES = NUMPY_SCALAR_TYPES + [float, complex, int, np.int32, np.uint64]

FLOATING_DTYPES = ["f4", "f8", "c8", "c16"]
DTYPES = FLOATING_DTYPES + ["i4", "u8"]


class TestScalar:
    @pytest.mark.parametrize("scalar_type", SCALAR_TYPES)
    def test_scalar(self, scalar_type):
        a = scalar_type(5)
        if issubclass(scalar_type, (complex, np.complexfloating)):
            a = a - 4.0j

        assert isinstance(a, scalar_type)
        scaler = Scaler(10.0)
        got = scaler(a)
        exp = a * 10.0
        assert got == exp

    @pytest.mark.parametrize("scalar_type", NUMPY_SCALAR_TYPES)
    def test_scalar_overflow(self, scalar_type):
        a = scalar_type(np.finfo(scalar_type).max)
        if issubclass(scalar_type, (complex, np.complexfloating)):
            a = a - 4.0j

        assert isinstance(a, scalar_type)
        scaler = Scaler(1000.0)
        with pytest.warns(RuntimeWarning, match="overflow.*scalar multiply"):
            scaler(a)


@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("swapbyteorder", [False, True])
class TestArray:
    def get_dtype(self, dtype, swapbyteorder):
        dtype = np.dtype(dtype)
        return dtype.newbyteorder() if swapbyteorder else dtype

    def get_array(self, data, dtype, order, imag=-4j):
        if dtype.kind == "c":
            data = np.asanyarray(data) + np.asanyarray(imag)
        data = np.array(data).astype(dtype, order=order)
        assert data.dtype == dtype
        return data

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_array(self, dtype, swapbyteorder, order):
        dtype = self.get_dtype(dtype, swapbyteorder)
        a = self.get_array(np.arange(10).reshape(2, 5), dtype, order)
        if order == "C":
            assert a.flags["C_CONTIGUOUS"] and not a.flags["F_CONTIGUOUS"]
        else:
            assert a.flags["F_CONTIGUOUS"] and not a.flags["C_CONTIGUOUS"]
        scaler = Scaler(10.0)
        got = scaler(a)
        exp = a * 10.0
        assert_array_equal(got, exp, strict=True)

    @pytest.mark.parametrize("dtype", FLOATING_DTYPES)
    def test_array_overflow(self, dtype, swapbyteorder, order):
        dtype = self.get_dtype(dtype, swapbyteorder)
        a = self.get_array([1.0, np.finfo(dtype).max], dtype, order)
        scaler = Scaler(1000.0)
        with pytest.warns(RuntimeWarning, match="overflow.*multiply"):
            scaler(a)


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


def test_scale_attribute():
    sc = Scaler(10.0)
    assert sc.factor == 10.0


def test_init_errors():
    with pytest.raises(TypeError, match=r"exactly 1 positional.*\(0 given\)"):
        Scaler()
    with pytest.raises(TypeError, match=r"exactly 1 positional.*\(0 given\)"):
        Scaler(b=20.0)
    with pytest.raises(TypeError, match=r"at most 1 argument \(2 given\)"):
        Scaler(10.0, 20.0)
    with pytest.raises(TypeError, match=r"at most 1 argument \(2 given\)"):
        Scaler(10.0, b=20.0)
    with pytest.raises(TypeError, match=r"at most 1 argument \(3 given\)"):
        Scaler(10.0, 20.0, b=20.0)
    with pytest.raises(TypeError, match=r"at most 1 argument \(3 given\)"):
        Scaler(10.0, b=20.0, c=20.0)
    with pytest.raises(TypeError, match="must be real number, not"):
        Scaler(1 + 1j)
    with pytest.raises(TypeError, match="must be real number, not"):
        Scaler({})


def test_call_errors():
    sc = Scaler(10.0)
    with pytest.raises(TypeError, match="takes 1 positional.*0 were given"):
        sc()
    with pytest.raises(TypeError, match="takes 1 positional.*2 were given"):
        sc(10.0, 20.0)
    with pytest.raises(TypeError, match="got an unexpected.*'b'"):
        sc(b=20.0)
    with pytest.raises(TypeError, match="got an unexpected.*'b'"):
        sc(10.0, b=20.0)
    with pytest.raises(TypeError, match="got an unexpected.*'b'"):
        sc(10.0, 20.0, b=20.0)


def test_unity_scaler_is_singleton():
    sc1a = Scaler(1.0)
    sc2a = Scaler(2.0)
    sc1b = Scaler(1.0)
    sc2b = Scaler(2.0)
    assert sc1a is sc1b
    assert sc1a is not sc2a
    assert sc2b is not sc2a


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

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
def test_with_array(dtype, order):
    a = np.arange(10, dtype=dtype).reshape(2, 5)
    if issubclass(dtype, np.complexfloating):
        a = a - 4.0j
    a = a.copy(order=order)

    assert a.dtype == dtype
    if order == "C":
        assert a.flags["C_CONTIGUOUS"] and not a.flags["F_CONTIGUOUS"]
    else:
        assert a.flags["F_CONTIGUOUS"] and not a.flags["C_CONTIGUOUS"]
    got = Scaler(10.0)(a)
    exp = a * 10.0
    assert got.dtype == exp.dtype
    assert_array_equal(got, exp)

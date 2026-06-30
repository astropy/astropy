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
def test_with_array(dtype):
    a = np.arange(10, dtype=dtype)
    if issubclass(dtype, np.complexfloating):
        a = a - 4.0j

    assert a.dtype == dtype
    got = Scaler(10.0)(a)
    exp = a * 10.0
    assert_array_equal(got, exp)

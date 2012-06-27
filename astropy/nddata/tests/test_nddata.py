# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from ..nddata import NDData

np.random.seed(12345)


def test_nddata_empty():

    with pytest.raises(TypeError) as exc:
        NDData()  # empty initializer should fail
    assert exc.value.args[0] == '__init__() takes at least 2 arguments (1 given)'


def test_nddata_simple():

    nd = NDData(np.random.random((10, 10)))
    assert nd.shape == (10, 10)
    assert nd.size == 100
    assert nd.dtype == np.dtype(float)


def test_nddata_error_valid_scalar():
    NDData(np.random.random((10, 10)), error=1.)


def test_nddata_error_valid_1d():
    NDData(np.random.random((10, 10)), error=np.ones((10,)))


def test_nddata_error_valid_2d():
    NDData(np.random.random((10, 10)), error=np.ones((10, 10)))


def test_nddata_error_valid_3d():
    NDData(np.random.random((10, 10)), error=np.ones((3, 10, 10)))


def test_nddata_error_invalid():
    with pytest.raises(ValueError) as exc:
        NDData(np.random.random((10, 10)), error=np.ones((3, 10)))
    assert exc.value.args[0] == 'dimensions of `error` do not match data'


def test_nddata_mask_valid():
    NDData(np.random.random((10, 10)), mask=np.random.random((10, 10)) > 0.5)


def test_nddata_mask_invalid_1():
    with pytest.raises(ValueError) as exc:
        NDData(np.random.random((10, 10)), mask=np.random.random((10)) > 0.5)
    assert exc.value.args[0] == 'dimensions of `mask` do not match data'

def test_nddata_mask_invalid_2():
    with pytest.raises(ValueError) as exc:
        NDData(np.random.random((10, 10)), mask=np.random.random((5, 5)) > 0.5)
    assert exc.value.args[0] == 'dimensions of `mask` do not match data'


def test_nddata_mask_invalid_3():
    with pytest.raises(ValueError) as exc:
        NDData(np.random.random((10, 10)), mask=np.random.random((3, 10, 10)) > 0.5)
    assert exc.value.args[0] == 'dimensions of `mask` do not match data'


def test_nddata_flags_valid_1d():
    NDData(np.random.random((10, 10)), flags=np.ones((10,)))


def test_nddata_flags_valid_2d():
    NDData(np.random.random((10, 10)), flags=np.ones((10, 10)))


def test_nddata_flags_valid_3d():
    NDData(np.random.random((10, 10)), flags=np.ones((3, 10, 10)))


def test_nddata_flags_invalid():
    with pytest.raises(ValueError) as exc:
        NDData(np.random.random((10, 10)), flags=np.ones((3, 10)))
    assert exc.value.args[0] == 'dimensions of `flags` do not match data'


def test_nddata_copy():
    a = np.ones((10, 10))
    nd_copy = NDData(a, copy=True)
    nd_ref = NDData(a, copy=False)
    a[0, 0] = 0
    assert nd_ref.data[0, 0] == 0
    assert nd_copy.data[0, 0] == 1


def test_nddata_conversion():
    nd = NDData([[1, 2, 3], [4, 5, 6]])
    assert nd.size == 6
    assert nd.dtype == np.dtype(int)


def test_boolmask():
    boolmask = np.random.random((10, 10)) > 0.5  # random mask that boolmask should look like
    nd1 = NDData(np.random.random((10, 10)), mask=~boolmask)  # mask False where valid
    assert np.all(nd1.boolmask == boolmask)

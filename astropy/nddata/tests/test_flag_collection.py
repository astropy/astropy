# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np
import pytest

from astropy.nddata import FlagCollection


def test_init():
    FlagCollection(shape=(1, 2, 3))


def test_init_noshape():
    with pytest.raises(Exception) as exc:
        FlagCollection()
    assert (
        exc.value.args[0]
        == "FlagCollection should be initialized with the shape of the data"
    )


def test_init_notiterable():
    with pytest.raises(Exception) as exc:
        FlagCollection(shape=1.0)
    assert exc.value.args[0] == "FlagCollection shape should be an iterable object"


def test_setitem():
    f = FlagCollection(shape=(1, 2, 3))
    f["a"] = np.ones((1, 2, 3)).astype(float)
    f["b"] = np.ones((1, 2, 3)).astype(int)
    f["c"] = np.ones((1, 2, 3)).astype(bool)
    f["d"] = np.ones((1, 2, 3)).astype(str)


@pytest.mark.parametrize("value", [1, 1.0, "spam", [1, 2, 3], (1.0, 2.0, 3.0)])
def test_setitem_invalid_type(value):
    f = FlagCollection(shape=(1, 2, 3))
    with pytest.raises(Exception) as exc:
        f["a"] = value
    assert exc.value.args[0] == "flags should be given as a Numpy array"


def test_setitem_invalid_shape():
    f = FlagCollection(shape=(1, 2, 3))
    with pytest.raises(ValueError) as exc:
        f["a"] = np.ones((3, 2, 1))
    assert exc.value.args[0].startswith("flags array shape")
    assert exc.value.args[0].endswith("does not match data shape (1, 2, 3)")

# Licensed under a 3-clause BSD style license - see LICENSE.rst
from numpy.testing.utils import assert_allclose, assert_equal
from ...tests.helper import pytest
from .. import models

def test_Gaussian1DModel():
    # TODO: add many more tests
    # see https://github.com/astropy/astropy/issues/1077
    m = models.Gaussian1DModel(42, 43, 44)
    assert m(43) == 42

@pytest.mark.xfail
def test_Gaussian2DModel():
    # TODO: add many more tests
    # see https://github.com/astropy/astropy/issues/1077
    m = models.Gaussian2DModel(42, 43, 44, 45, 46)
    assert m(43, 44) == 42

def test_ShiftModel():
    # Shift by a scalar
    m = models.ShiftModel(42)
    assert m(0) == 42
    assert_equal(m([1, 2]), [43, 44])

    # Shift by a list
    m = models.ShiftModel([42, 43])
    assert_equal(m(0), [42, 43])
    assert_equal(m([1, 2]), [[ 43,  44], [ 44,  45]])

def test_ScaleModel():
    # Scale by a scalar
    m = models.ScaleModel(42)
    assert m(0) == 0
    assert_equal(m([1, 2]), [42, 84])

    # Scale by a list
    m = models.ScaleModel([42, 43])
    assert_equal(m(0), [0, 0])
    assert_equal(m([1, 2]), [[ 42,  43], [ 84,  86]])

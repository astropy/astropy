# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_almost_equal_nulp

from ..biweight import (biweight_location, biweight_scale,
                        biweight_midvariance, biweight_midcovariance,
                        biweight_midcorrelation)
from ...tests.helper import catch_warnings
from ...utils.misc import NumpyRNGContext


def test_biweight_location():
    with NumpyRNGContext(12345):
        # test that it runs
        randvar = np.random.randn(10000)
        cbl = biweight_location(randvar)
        assert abs(cbl - 0) < 1e-2


def test_biweight_location_constant():
    cbl = biweight_location(np.ones((10, 5)))
    assert cbl == 1.


def test_biweight_location_constant_axis_2d():
    shape = (10, 5)
    data = np.ones(shape)
    cbl = biweight_location(data, axis=0)
    assert_allclose(cbl, np.ones(shape[1]))
    cbl = biweight_location(data, axis=1)
    assert_allclose(cbl, np.ones(shape[0]))

    val1 = 100.
    val2 = 2.
    data = np.arange(50).reshape(10, 5)
    data[2] = val1
    data[7] = val2
    cbl = biweight_location(data, axis=1)
    assert_allclose(cbl[2], val1)
    assert_allclose(cbl[7], val2)


def test_biweight_location_constant_axis_3d():
    shape = (10, 5, 2)
    data = np.ones(shape)
    cbl = biweight_location(data, axis=0)
    assert_allclose(cbl, np.ones((shape[1], shape[2])))
    cbl = biweight_location(data, axis=1)
    assert_allclose(cbl, np.ones((shape[0], shape[2])))
    cbl = biweight_location(data, axis=2)
    assert_allclose(cbl, np.ones((shape[0], shape[1])))


def test_biweight_location_small():
    cbl = biweight_location([1, 3, 5, 500, 2])
    assert abs(cbl - 2.745) < 1e-3


def test_biweight_location_axis():
    """Test a 2D array with the axis keyword."""
    with NumpyRNGContext(12345):
        ny = 100
        nx = 200
        data = np.random.normal(5, 2, (ny, nx))

        bw = biweight_location(data, axis=0)
        bwi = []
        for i in range(nx):
            bwi.append(biweight_location(data[:, i]))
        bwi = np.array(bwi)
        assert_allclose(bw, bwi)

        bw = biweight_location(data, axis=1)
        bwi = []
        for i in range(ny):
            bwi.append(biweight_location(data[i, :]))
        bwi = np.array(bwi)
        assert_allclose(bw, bwi)


def test_biweight_location_axis_3d():
    """Test a 3D array with the axis keyword."""
    with NumpyRNGContext(12345):
        nz = 3
        ny = 4
        nx = 5
        data = np.random.normal(5, 2, (nz, ny, nx))
        bw = biweight_location(data, axis=0)
        assert bw.shape == (ny, nx)

        y = 0
        bwi = []
        for i in range(nx):
            bwi.append(biweight_location(data[:, y, i]))
        bwi = np.array(bwi)
        assert_allclose(bw[y], bwi)


def test_biweight_scale():
    # NOTE:  biweight_scale is covered by biweight_midvariance tests
    data = [1, 3, 5, 500, 2]
    scl = biweight_scale(data)
    var = biweight_midvariance(data)
    assert_allclose(scl, np.sqrt(var))


def test_biweight_midvariance():
    with NumpyRNGContext(12345):
        # test that it runs
        randvar = np.random.randn(10000)
        var = biweight_midvariance(randvar)
        assert_allclose(var, 1.0, rtol=0.02)


def test_biweight_midvariance_small():
    data = [1, 3, 5, 500, 2]
    var = biweight_midvariance(data)
    assert_allclose(var, 2.9238456)    # verified with R

    var = biweight_midvariance(data, modify_sample_size=True)
    assert_allclose(var, 2.3390765)


def test_biweight_midvariance_5127():
    # test a regression introduced in #5127
    rand = np.random.RandomState(12345)
    data = rand.normal(loc=0., scale=20., size=(100, 100))
    var = biweight_midvariance(data)
    assert_allclose(var, 406.86938710817344)    # verified with R


def test_biweight_midvariance_axis():
    """Test a 2D array with the axis keyword."""
    with NumpyRNGContext(12345):
        ny = 100
        nx = 200
        data = np.random.normal(5, 2, (ny, nx))

        bw = biweight_midvariance(data, axis=0)
        bwi = []
        for i in range(nx):
            bwi.append(biweight_midvariance(data[:, i]))
        bwi = np.array(bwi)
        assert_allclose(bw, bwi)

        bw = biweight_midvariance(data, axis=1)
        bwi = []
        for i in range(ny):
            bwi.append(biweight_midvariance(data[i, :]))
        bwi = np.array(bwi)
        assert_allclose(bw, bwi)


def test_biweight_midvariance_axis_3d():
    """Test a 3D array with the axis keyword."""
    with NumpyRNGContext(12345):
        nz = 3
        ny = 4
        nx = 5
        data = np.random.normal(5, 2, (nz, ny, nx))
        bw = biweight_midvariance(data, axis=0)
        assert bw.shape == (ny, nx)

        y = 0
        bwi = []
        for i in range(nx):
            bwi.append(biweight_midvariance(data[:, y, i]))
        bwi = np.array(bwi)
        assert_allclose(bw[y], bwi)


def test_biweight_midvariance_constant_axis():
    bw = biweight_midvariance(np.ones((10, 5)))
    assert bw == 0.0


def test_biweight_midvariance_constant_axis_2d():
    shape = (10, 5)
    data = np.ones(shape)
    cbl = biweight_midvariance(data, axis=0)
    assert_allclose(cbl, np.zeros(shape[1]))
    cbl = biweight_midvariance(data, axis=1)
    assert_allclose(cbl, np.zeros(shape[0]))

    data = np.arange(50).reshape(10, 5)
    data[2] = 100.
    data[7] = 2.
    bw = biweight_midvariance(data, axis=1)
    assert_allclose(bw[2], 0.)
    assert_allclose(bw[7], 0.)


def test_biweight_midvariance_constant_axis_3d():
    shape = (10, 5, 2)
    data = np.ones(shape)
    cbl = biweight_midvariance(data, axis=0)
    assert_allclose(cbl, np.zeros((shape[1], shape[2])))
    cbl = biweight_midvariance(data, axis=1)
    assert_allclose(cbl, np.zeros((shape[0], shape[2])))
    cbl = biweight_midvariance(data, axis=2)
    assert_allclose(cbl, np.zeros((shape[0], shape[1])))


def test_biweight_midcovariance_1d():
    d = [0, 1, 2]
    cov = biweight_midcovariance(d)
    var = biweight_midvariance(d)
    assert_allclose(cov, [[var]])


def test_biweight_midcovariance_2d():
    d = [[0, 1, 2], [2, 1, 0]]
    cov = biweight_midcovariance(d)
    val = 0.70121809
    assert_allclose(cov, [[val, -val], [-val, val]])    # verified with R

    d = [[5, 1, 10], [500, 5, 2]]
    cov = biweight_midcovariance(d)
    assert_allclose(cov, [[14.54159077, -7.79026256],    # verified with R
                          [-7.79026256, 6.92087252]])

    cov = biweight_midcovariance(d, modify_sample_size=True)
    assert_allclose(cov, [[14.54159077, -5.19350838],
                          [-5.19350838, 4.61391501]])


def test_biweight_midcovariance_constant():
    data = np.ones((3, 10))
    cov = biweight_midcovariance(data)
    assert_allclose(cov, np.zeros((3, 3)))


def test_biweight_midcovariance_midvariance():
    """
    Test that biweight_midcovariance diagonal elements agree with
    biweight_midvariance.
    """

    rng = np.random.RandomState(1)
    d = rng.normal(0, 2, size=(100, 3))
    cov = biweight_midcovariance(d)
    var = [biweight_midvariance(a) for a in d]
    assert_allclose(cov.diagonal(), var)

    cov2 = biweight_midcovariance(d, modify_sample_size=True)
    var2 = [biweight_midvariance(a, modify_sample_size=True)
            for a in d]
    assert_allclose(cov2.diagonal(), var2)


def test_midcovariance_shape():
    """
    Test that biweight_midcovariance raises error with a 3D array.
    """

    d = np.ones(27).reshape(3, 3, 3)
    with pytest.raises(ValueError) as e:
        biweight_midcovariance(d)
    assert 'The input array must be 2D or 1D.' in str(e.value)


def test_midcovariance_M_shape():
    """
    Test that biweight_midcovariance raises error when M is not a scalar
    or 1D array.
    """

    d = [0, 1, 2]
    M = [[0, 1], [2, 3]]
    with pytest.raises(ValueError) as e:
        biweight_midcovariance(d, M=M)
    assert 'M must be a scalar or 1D array.' in str(e.value)


def test_biweight_midcovariance_symmetric():
    """
    Regression test to ensure that midcovariance matrix is symmetric
    when ``modify_sample_size=True`` (see #5972).
    """

    rng = np.random.RandomState(1)
    d = rng.gamma(2, 2, size=(3, 500))
    cov = biweight_midcovariance(d)
    assert_array_almost_equal_nulp(cov, cov.T, nulp=5)

    cov = biweight_midcovariance(d, modify_sample_size=True)
    assert_array_almost_equal_nulp(cov, cov.T, nulp=5)


def test_biweight_midcorrelation():
    x = [0, 1, 2]
    y = [2, 1, 0]
    assert_allclose(biweight_midcorrelation(x, x), 1.0)
    assert_allclose(biweight_midcorrelation(x, y), -1.0)

    x = [5, 1, 10, 12.4, 13.2]
    y = [500, 5, 2, 7.1, 0.9]
    # verified with R
    assert_allclose(biweight_midcorrelation(x, y), -0.14411038976763313)


def test_biweight_midcorrelation_inputs():
    a1 = np.ones((3, 3))
    a2 = np.ones(5)
    a3 = np.ones(7)

    with pytest.raises(ValueError) as e:
        biweight_midcorrelation(a1, a2)
        assert 'x must be a 1D array.' in str(e.value)

    with pytest.raises(ValueError) as e:
        biweight_midcorrelation(a2, a1)
        assert 'y must be a 1D array.' in str(e.value)

    with pytest.raises(ValueError) as e:
        biweight_midcorrelation(a2, a3)
        assert 'x and y must have the same shape.' in str(e.value)


def test_biweight_32bit_runtime_warnings():
    """Regression test for #6905."""
    with NumpyRNGContext(12345):
        data = np.random.random(100).astype(np.float32)
        data[50] = 30000.

        with catch_warnings(RuntimeWarning) as warning_lines:
            biweight_scale(data)
            assert len(warning_lines) == 0

        with catch_warnings(RuntimeWarning) as warning_lines:
            biweight_midvariance(data)
            assert len(warning_lines) == 0

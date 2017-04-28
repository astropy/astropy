# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.random import randn, normal
from numpy.testing import assert_equal
from numpy.testing.utils import assert_allclose

from ..biweight import (biweight_location, biweight_midvariance,
                        biweight_midcovariance)
from ...tests.helper import pytest
from ...extern.six.moves import range
from ...utils.misc import NumpyRNGContext


def test_biweight_location():
    with NumpyRNGContext(12345):
        # test that it runs
        randvar = randn(10000)
        cbl = biweight_location(randvar)
        assert abs(cbl - 0) < 1e-2


def test_biweight_location_small():
    cbl = biweight_location([1, 3, 5, 500, 2])
    assert abs(cbl - 2.745) < 1e-3


def test_biweight_location_axis():
    """Test a 2D array with the axis keyword."""
    with NumpyRNGContext(12345):
        ny = 100
        nx = 200
        data = normal(5, 2, (ny, nx))

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
        data = normal(5, 2, (nz, ny, nx))
        bw = biweight_location(data, axis=0)
        assert bw.shape == (ny, nx)

        y = 0
        bwi = []
        for i in range(nx):
            bwi.append(biweight_location(data[:, y, i]))
        bwi = np.array(bwi)
        assert_allclose(bw[y], bwi)


def test_biweight_midvariance():
    with NumpyRNGContext(12345):
        # test that it runs
        randvar = randn(10000)
        scl = biweight_midvariance(randvar)
        assert_allclose(scl, 1.0, rtol=0.02)


def test_biweight_midvariance_small():
    data = [1, 3, 5, 500, 2]
    scl = biweight_midvariance([1, 3, 5, 500, 2])
    assert_allclose(scl, 2.9238456)    # verified with R

    scl = biweight_midvariance(data, modify_sample_size=True)
    assert_allclose(scl, 2.3390765)


def test_biweight_midvariance_5127():
    # test a regression introduced in #5127
    rand = np.random.RandomState(12345)
    data = rand.normal(loc=0., scale=20., size=(100, 100))
    scl = biweight_midvariance(data)
    assert_allclose(scl, 406.86938710817344)    # verified with R


def test_biweight_midvariance_axis():
    """Test a 2D array with the axis keyword."""
    with NumpyRNGContext(12345):
        ny = 100
        nx = 200
        data = normal(5, 2, (ny, nx))

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
        data = normal(5, 2, (nz, ny, nx))
        bw = biweight_midvariance(data, axis=0)
        assert bw.shape == (ny, nx)

        y = 0
        bwi = []
        for i in range(nx):
            bwi.append(biweight_midvariance(data[:, y, i]))
        bwi = np.array(bwi)
        assert_allclose(bw[y], bwi)


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


def test_biweight_midcovariance_symmetric():
    """
    Regression test to ensure that midcovariance matrix is symmetric
    when ``modify_sample_size=True`` (see #5972).
    """

    rng = np.random.RandomState(1)
    d = rng.gamma(2, 2, size=(3, 500))
    cov = biweight_midcovariance(d)
    assert_equal(cov, cov.T)

    cov = biweight_midcovariance(d, modify_sample_size=True)
    assert_equal(cov, cov.T)

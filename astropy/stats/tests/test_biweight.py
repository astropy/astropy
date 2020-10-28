# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import (assert_allclose, assert_array_almost_equal_nulp,
                           assert_equal)

from astropy.stats.biweight import (biweight_location, biweight_scale,
                                    biweight_midvariance,
                                    biweight_midcovariance,
                                    biweight_midcorrelation)
from astropy.utils.misc import NumpyRNGContext


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
    val3 = 5.
    data = np.arange(50).reshape(10, 5)
    data[2] = val1
    data[7] = val2
    data[8] = [val3, 0.8, val3, -0.8, val3]
    cbl = biweight_location(data, axis=1)
    assert_allclose(cbl[2], val1)
    assert_allclose(cbl[7], val2)
    assert_allclose(cbl[8], val3)


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
    bw_loc = biweight_location([1, 3, 5, 500, 2])
    assert_allclose(bw_loc, 2.7456117)


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


def test_biweight_location_axis_tuple():
    """Test a 3D array with a tuple axis keyword."""

    data = np.arange(24).reshape(2, 3, 4)
    data[0, 0] = 100.

    assert_equal(biweight_location(data, axis=0),
                 biweight_location(data, axis=(0,)))
    assert_equal(biweight_location(data, axis=-1),
                 biweight_location(data, axis=(2,)))
    assert_equal(biweight_location(data, axis=(0, 1)),
                 biweight_location(data, axis=(1, 0)))
    assert_equal(biweight_location(data, axis=(0, 2)),
                 biweight_location(data, axis=(0, -1)))
    assert_equal(biweight_location(data, axis=(0, 1, 2)),
                 biweight_location(data, axis=(2, 0, 1)))
    assert_equal(biweight_location(data, axis=(0, 1, 2)),
                 biweight_location(data, axis=None))


@pytest.mark.filterwarnings('ignore:All-NaN slice encountered')
@pytest.mark.filterwarnings('ignore:Invalid value encountered in median')
def test_biweight_location_ignore_nan():
    data1d = np.array([1, 3, 5, 500, 2, np.nan])
    data2d = np.array([data1d, data1d])

    assert np.isnan(biweight_location(data1d, ignore_nan=False))

    biw_expected = biweight_location(data1d[:-1], ignore_nan=False)
    assert_equal(biweight_location(data1d, ignore_nan=True), biw_expected)

    assert_equal(biweight_location(data2d, axis=0, ignore_nan=True),
                 data1d)
    assert_equal(biweight_location(data2d, axis=1, ignore_nan=True),
                 [biw_expected, biw_expected])


@pytest.mark.filterwarnings('ignore:All-NaN slice encountered')
@pytest.mark.filterwarnings('ignore:Invalid value encountered in median')
def test_biweight_location_nan():
    data1d = np.array([1, 3, 5, 500, 2, np.nan])
    all_nan = data1d.copy()
    all_nan[:] = np.nan
    data2d = np.array([data1d, data1d, all_nan])
    data1d_masked = np.ma.masked_invalid(data1d)
    data1d_masked.data[0] = np.nan
    data2d_masked = np.ma.masked_invalid(data2d)

    assert np.isnan(biweight_location(data1d))
    bw_loc = biweight_location(data1d_masked)
    assert not isinstance(bw_loc, np.ma.MaskedArray)
    assert np.isnan(biweight_location(data2d))

    for axis in (0, 1):
        assert np.all(np.isnan(biweight_location(data2d, axis=axis)))
        assert isinstance(biweight_location(data2d_masked, axis=axis),
                          np.ma.MaskedArray)


@pytest.mark.filterwarnings('ignore:All-NaN slice encountered')
@pytest.mark.filterwarnings('ignore:Invalid value encountered in median')
def test_biweight_location_masked():
    data1d = np.array([1, 3, 5, 500, 2, np.nan])
    data2d = np.array([data1d, data1d])

    data1d_masked = np.ma.masked_invalid(data1d)
    data2d_masked = np.ma.masked_invalid(data2d)

    assert_equal(biweight_location(data1d, ignore_nan=True),
                 biweight_location(data1d_masked))
    assert_equal(biweight_location(data2d, ignore_nan=True),
                 biweight_location(data2d_masked))

    bw_loc = biweight_location(data1d_masked)
    assert_allclose(bw_loc, 2.7456117)
    assert np.isscalar(bw_loc)

    bw_loc = biweight_location(data2d, ignore_nan=True, axis=1)
    bw_loc_masked = biweight_location(data2d_masked, axis=1)
    assert isinstance(bw_loc_masked, np.ma.MaskedArray)
    assert ~np.any(bw_loc_masked.mask)  # mask is all False
    assert_equal(bw_loc, bw_loc_masked.data)

    bw_loc = biweight_location(data2d, ignore_nan=True, axis=0)
    bw_loc_masked = biweight_location(data2d_masked, axis=0)
    assert_equal(bw_loc_masked.data[:-1], bw_loc[:-1])
    assert bw_loc_masked.mask[-1]  # last mask element is True

    data1d_masked.data[0] = np.nan  # unmasked NaN
    bw_loc = biweight_location(data1d_masked)
    assert not isinstance(bw_loc, np.ma.MaskedArray)
    assert np.isscalar(bw_loc)
    assert np.isnan(bw_loc)
    assert_equal(biweight_location(data1d_masked, ignore_nan=True),
                 biweight_location(data1d[1:], ignore_nan=True))

    # ensure that input masked array is not modified
    assert np.isnan(data1d_masked[0])


def test_biweight_scale():
    # NOTE:  biweight_scale is covered by biweight_midvariance tests
    data = [1, 3, 5, 500, 2]
    scl = biweight_scale(data)
    var = biweight_midvariance(data)
    assert_allclose(scl, np.sqrt(var))

    data = np.ma.masked_invalid([1, 3, 5, 500, 2, np.nan])
    data[0] = np.nan
    scl = biweight_scale(data, ignore_nan=True)
    var = biweight_midvariance(data, ignore_nan=True)
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
    assert_allclose(var, 2.9238456)  # verified with R

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


@pytest.mark.filterwarnings('ignore:All-NaN slice encountered')
@pytest.mark.filterwarnings('ignore:Invalid value encountered in median')
def test_biweight_midvariance_ignore_nan():
    data1d = np.array([1, 3, 5, 500, 2, np.nan])
    data2d = np.array([data1d, data1d])

    assert np.isnan(biweight_midvariance(data1d, ignore_nan=False))

    biw_var = biweight_midvariance(data1d[:-1], ignore_nan=False)
    biw_var_nonan = biweight_midvariance(data1d, ignore_nan=True)
    assert_equal(biw_var_nonan, biw_var)

    assert_equal(biweight_midvariance(data2d, axis=0, ignore_nan=True),
                 [0., 0., 0., 0., 0., np.nan])
    assert_equal(biweight_midvariance(data2d, axis=1, ignore_nan=True),
                 [biw_var_nonan, biw_var_nonan])


@pytest.mark.filterwarnings('ignore:All-NaN slice encountered')
@pytest.mark.filterwarnings('ignore:Invalid value encountered in median')
def test_biweight_scale_nan():
    data1d = np.array([1, 3, 5, 500, 2, np.nan])
    all_nan = data1d.copy()
    all_nan[:] = np.nan
    data2d = np.array([data1d, data1d, all_nan])
    data1d_masked = np.ma.masked_invalid(data1d)
    data1d_masked.data[0] = np.nan
    data2d_masked = np.ma.masked_invalid(data2d)

    assert np.isnan(biweight_scale(data1d))
    bw_scl = biweight_scale(data1d_masked)
    assert not isinstance(bw_scl, np.ma.MaskedArray)
    assert np.isnan(bw_scl)
    assert np.isnan(biweight_scale(data2d))
    assert_allclose(biweight_scale(data2d_masked), 1.709926, atol=1e-5)

    for axis in (0, 1):
        assert np.all(np.isnan(biweight_scale(data2d, axis=axis)))
        assert isinstance(biweight_scale(data2d_masked, axis=axis),
                          np.ma.MaskedArray)


@pytest.mark.filterwarnings('ignore:All-NaN slice encountered')
@pytest.mark.filterwarnings('ignore:Invalid value encountered in median')
def test_biweight_midvariance_masked():
    data1d = np.array([1, 3, 5, 500, 2, np.nan])
    data2d = np.array([data1d, data1d])

    data1d_masked = np.ma.masked_invalid(data1d)
    data2d_masked = np.ma.masked_invalid(data2d)

    assert_equal(biweight_midvariance(data1d, ignore_nan=True),
                 biweight_midvariance(data1d_masked))
    assert_equal(biweight_midvariance(data2d, ignore_nan=True),
                 biweight_midvariance(data2d_masked))

    bw_scl = biweight_midvariance(data1d_masked)
    assert_allclose(bw_scl, 2.9238456)
    assert np.isscalar(bw_scl)

    bw_loc = biweight_midvariance(data2d, ignore_nan=True, axis=1)
    bw_loc_masked = biweight_midvariance(data2d_masked, axis=1)
    assert isinstance(bw_loc_masked, np.ma.MaskedArray)
    assert ~np.any(bw_loc_masked.mask)  # mask is all False
    assert_equal(bw_loc, bw_loc_masked.data)

    bw_loc = biweight_midvariance(data2d, ignore_nan=True, axis=0)
    bw_loc_masked = biweight_midvariance(data2d_masked, axis=0)
    assert_equal(bw_loc_masked.data[:-1], bw_loc[:-1])
    assert bw_loc_masked.mask[-1]  # last mask element is True

    data1d_masked.data[0] = np.nan  # unmasked NaN
    bw_scl = biweight_midvariance(data1d_masked)
    assert not isinstance(bw_scl, np.ma.MaskedArray)
    assert np.isscalar(bw_scl)
    assert np.isnan(bw_scl)
    assert_equal(biweight_midvariance(data1d_masked, ignore_nan=True),
                 biweight_midvariance(data1d[1:], ignore_nan=True))

    # ensure that input masked array is not modified
    assert np.isnan(data1d_masked[0])


def test_biweight_scale_axis_tuple():
    """Test a 3D array with a tuple axis keyword."""

    data = np.arange(24).reshape(2, 3, 4)
    data[0, 0] = 100.

    assert_equal(biweight_scale(data, axis=0),
                 biweight_scale(data, axis=(0,)))
    assert_equal(biweight_scale(data, axis=-1),
                 biweight_scale(data, axis=(2,)))
    assert_equal(biweight_scale(data, axis=(0, 1)),
                 biweight_scale(data, axis=(1, 0)))
    assert_equal(biweight_scale(data, axis=(0, 2)),
                 biweight_scale(data, axis=(0, -1)))
    assert_equal(biweight_scale(data, axis=(0, 1, 2)),
                 biweight_scale(data, axis=(2, 0, 1)))
    assert_equal(biweight_scale(data, axis=(0, 1, 2)),
                 biweight_scale(data, axis=None))
    assert_equal(biweight_scale(data, axis=(0, 2), modify_sample_size=True),
                 biweight_scale(data, axis=(0, -1), modify_sample_size=True))


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
    data[8] = [5.0, 0.8, 5.0, -0.8, 5.0]
    bw = biweight_midvariance(data, axis=1)
    assert_allclose(bw[2], 0.)
    assert_allclose(bw[7], 0.)
    assert_allclose(bw[8], 0.)


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
    val3 = 5.0
    data[1] = [val3, 0.8, val3, -0.8, val3, val3, val3, 1.0, val3, -0.7]
    cov = biweight_midcovariance(data)
    assert_allclose(cov, np.zeros((3, 3)))

    rng = np.random.default_rng(123)
    data = rng.random((5, 5))
    val3 = 5.0
    data[1] = [val3, 0.8, val3, -0.8, val3]
    cov = biweight_midcovariance(data)
    assert_allclose(cov[1, :], 0.)
    assert_allclose(cov[:, 1], 0.)


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
        biweight_scale(data)
        biweight_midvariance(data)

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.coordinates import Angle
from astropy.stats import mad_std
from astropy.stats.sigma_clipping import (
    SigmaClip,
    SigmaClippedStats,
    sigma_clip,
    sigma_clipped_stats,
)
from astropy.table import MaskedColumn
from astropy.utils.compat import COPY_IF_NEEDED
from astropy.utils.compat.optional_deps import HAS_BOTTLENECK, HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.misc import NumpyRNGContext


def test_sigma_clip():
    # need to seed the numpy RNG to make sure we don't get some
    # amazingly flukey random number that breaks one of the tests

    with NumpyRNGContext(12345):
        # Amazing, I've got the same combination on my luggage!
        randvar = np.random.randn(10000)

        filtered_data = sigma_clip(randvar, sigma=1, maxiters=2)

        assert sum(filtered_data.mask) > 0
        assert sum(~filtered_data.mask) < randvar.size

        # this is actually a silly thing to do, because it uses the
        # standard deviation as the variance, but it tests to make sure
        # these arguments are actually doing something
        filtered_data2 = sigma_clip(randvar, sigma=1, maxiters=2, stdfunc=np.var)
        assert not np.all(filtered_data.mask == filtered_data2.mask)

        filtered_data3 = sigma_clip(randvar, sigma=1, maxiters=2, cenfunc=np.mean)
        assert not np.all(filtered_data.mask == filtered_data3.mask)

        # make sure the maxiters=None method works at all.
        filtered_data = sigma_clip(randvar, sigma=3, maxiters=None)

        # test copying
        assert filtered_data.data[0] == randvar[0]
        filtered_data.data[0] += 1.0
        assert filtered_data.data[0] != randvar[0]

        filtered_data = sigma_clip(randvar, sigma=3, maxiters=None, copy=False)
        assert filtered_data.data[0] == randvar[0]
        filtered_data.data[0] += 1.0
        assert filtered_data.data[0] == randvar[0]

        # test axis
        data = np.arange(5) + np.random.normal(0.0, 0.05, (5, 5)) + np.diag(np.ones(5))
        filtered_data = sigma_clip(data, axis=0, sigma=2.3)
        assert filtered_data.count() == 20
        filtered_data = sigma_clip(data, axis=1, sigma=2.3)
        assert filtered_data.count() == 25


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_axis_none():
    """
    For masked=False and axis=None, masked elements should be removed
    from the result.
    """
    data = np.arange(10.0)
    data[0] = 100
    result = sigma_clip(data, masked=False, axis=None)
    assert_equal(result, data[1:])


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_compare_to_scipy_sigmaclip():
    from scipy import stats

    # need to seed the numpy RNG to make sure we don't get some
    # amazingly flukey random number that breaks one of the tests
    with NumpyRNGContext(12345):
        randvar = np.random.randn(10000)

        astropyres = sigma_clip(randvar, sigma=3, maxiters=None, cenfunc=np.mean)
        scipyres = stats.sigmaclip(randvar, 3, 3)[0]

        assert astropyres.count() == len(scipyres)
        assert_equal(astropyres[~astropyres.mask].data, scipyres)


def test_sigma_clip_scalar_mask():
    """Test that the returned mask is not a scalar."""
    data = np.arange(5)
    result = sigma_clip(data, sigma=100.0, maxiters=1)
    assert result.mask.shape != ()


def test_sigma_clip_masked_array():
    """
    Regression test to check that MaskedArray masks are propagated
    correctly if cenfunc/stdfunc is not the default, axis is not None,
    and masked=False.
    """
    arr = np.ones(10).astype(float)
    arr[-2:] = 5
    arr = np.ma.masked_where(arr > 1, arr)
    out = sigma_clip(arr, cenfunc=np.mean, axis=0, masked=False)
    assert np.all(np.isnan(out[-2:]))
    assert_allclose(np.nanmean(out), 1.0)

    out = sigma_clip(arr, stdfunc=np.std, axis=0, masked=False)
    assert np.all(np.isnan(out[-2:]))
    assert_allclose(np.nanmean(out), 1.0)


def test_sigma_clip_class():
    with NumpyRNGContext(12345):
        data = np.random.randn(100)
        data[10] = 1.0e5
        sobj = SigmaClip(sigma=1, maxiters=2)
        sfunc = sigma_clip(data, sigma=1, maxiters=2)
        assert_equal(sobj(data), sfunc)


def test_sigma_clip_mean():
    with NumpyRNGContext(12345):
        data = np.random.normal(0.0, 0.05, (10, 10))
        data[2, 2] = 1.0e5
        sobj1 = SigmaClip(sigma=1, maxiters=2, cenfunc="mean")
        sobj2 = SigmaClip(sigma=1, maxiters=2, cenfunc=np.nanmean)
        assert_equal(sobj1(data), sobj2(data))
        assert_equal(sobj1(data, axis=0), sobj2(data, axis=0))


def test_sigma_clip_invalid_cenfunc_stdfunc():
    with pytest.raises(ValueError):
        SigmaClip(cenfunc="invalid")

    with pytest.raises(ValueError):
        SigmaClip(stdfunc="invalid")


def test_sigma_clipped_stats():
    """Test list data with input mask or mask_value (#3268)."""
    # test list data with mask
    data = [0, 1]
    mask = np.array([True, False])
    result = sigma_clipped_stats(data, mask=mask)
    # Check that the result of np.ma.median was converted to a scalar
    assert isinstance(result[1], float)
    assert result == (1.0, 1.0, 0.0)

    result2 = sigma_clipped_stats(data, mask=mask, axis=0)
    assert_equal(result, result2)

    # test list data with mask_value
    result = sigma_clipped_stats(data, mask_value=0.0)
    assert isinstance(result[1], float)
    assert result == (1.0, 1.0, 0.0)

    # test without mask
    data = [0, 2]
    result = sigma_clipped_stats(data)
    assert isinstance(result[1], float)
    assert result == (1.0, 1.0, 1.0)

    _data = np.arange(10)
    data = np.ma.MaskedArray([_data, _data, 10 * _data])
    mean = sigma_clip(data, axis=0, sigma=1).mean(axis=0)
    assert_equal(mean, _data)
    mean, median, stddev = sigma_clipped_stats(data, axis=0, sigma=1)
    assert_equal(mean, _data)
    assert_equal(median, _data)
    assert_equal(stddev, np.zeros_like(_data))


def test_sigma_clipped_stats_ddof():
    with NumpyRNGContext(12345):
        data = np.random.randn(10000)
        data[10] = 1.0e5
        mean1, median1, stddev1 = sigma_clipped_stats(data)
        mean2, median2, stddev2 = sigma_clipped_stats(data, std_ddof=1)
        assert mean1 == mean2
        assert median1 == median2
        assert_allclose(stddev1, 0.98156805711673156)
        assert_allclose(stddev2, 0.98161731654802831)


def test_sigma_clipped_stats_masked_col():
    # see https://github.com/astropy/astropy/issues/13281
    arr = np.ma.masked_array([1, 2, 3], mask=[False, True, False])
    sigma_clipped_stats(arr)

    col = MaskedColumn(data=arr)
    sigma_clipped_stats(col)


@pytest.mark.parametrize("units", [False, True])
@pytest.mark.parametrize("sigma", np.linspace(2, 5, 10))
def test_sigmaclippedstats_stats(units, sigma):
    """Test SigmaClippedStats class."""
    data = np.ones((101, 205))
    rng = np.random.default_rng(0)
    idx = rng.integers(0, 100, 10)
    data[idx, idx] = 1000

    if units:
        unit = u.m
        data <<= unit
    else:
        unit = 1

    stats = SigmaClippedStats(data, sigma=sigma)
    assert stats.min() == 1.0 * unit
    assert stats.max() == 1.0 * unit
    assert stats.sum() == np.sum(data) - (1000 * 10) * unit
    assert stats.mean() == 1.0 * unit
    assert stats.median() == 1.0 * unit
    assert stats.mode() == 1.0 * unit
    assert stats.mode(median_factor=3.5, mean_factor=2.5) == 1.0 * unit
    assert stats.std() == 0.0 * unit
    assert stats.var() == 0.0 * unit**2
    assert stats.mad_std() == 0.0 * unit
    assert stats.biweight_location() == 1.0 * unit
    assert stats.biweight_scale() == 0.0 * unit

    # test nanvar with Angle
    data = Angle([10, 20, 30, 40], "deg")
    stats = SigmaClippedStats(data, sigma=5)
    assert isinstance(stats.mean(), u.Quantity)
    assert isinstance(stats.var(), u.Quantity)
    assert not isinstance(stats.var(), Angle)
    assert stats.mean() == 25 * u.deg
    assert stats.var() == 125 * u.deg**2


def test_sigmaclippedstats_mask():
    data = np.ones((101, 205))
    rng = np.random.default_rng(0)
    idx = rng.integers(0, 100, 10)
    data[idx, idx] = 1000
    stats = SigmaClippedStats(data, sigma=3.2)
    assert stats.min() == 1.0
    assert stats.max() == 1.0
    assert stats.sum() == np.sum(data) - 1000 * 10
    assert stats.mean() == 1.0
    assert stats.median() == 1.0
    assert stats.mode() == 1.0
    assert stats.mode(median_factor=3.5, mean_factor=2.5) == 1.0
    assert stats.std() == 0.0
    assert stats.var() == 0.0
    assert stats.mad_std() == 0.0
    assert stats.biweight_location() == 1.0
    assert stats.biweight_scale() == 0.0

    mask = np.zeros_like(data, dtype=bool)
    mask[idx, idx] = True
    stats3 = SigmaClippedStats(data, sigma=3.2, mask=mask)
    assert stats3.mode() == stats.mode()

    mask = np.zeros_like(data, dtype=bool)
    stats4 = SigmaClippedStats(data, sigma=3.2, mask_value=1000)
    assert stats4.mode() == stats3.mode()

    data = np.ones((101, 205))
    mask = np.ones_like(data, dtype=bool)
    match = "input data is all masked"
    with pytest.raises(ValueError, match=match):
        SigmaClippedStats(data, sigma=3, mask=mask)


def test_sigmaclippedstats_axis():
    data = np.ones((101, 205))
    stats = SigmaClippedStats(data, sigma=2.8, axis=1)
    shape = (data.shape[0],)
    assert stats.mean().shape == shape
    assert stats.mode().shape == shape
    assert stats.var().shape == shape
    assert stats.mad_std().shape == shape
    assert stats.biweight_location().shape == shape


@pytest.mark.slow
@pytest.mark.skipif(
    not HAS_BOTTLENECK,
    reason="test a workaround for upstream bug in bottleneck",
)
@pytest.mark.parametrize("shape", [(1024, 1024), (6388, 9576)])
def test_sigma_clip_large_float32_arrays(shape):
    # see https://github.com/astropy/astropy/issues/17185
    rng = np.random.default_rng(0)

    expected = (0.5, 0.5, 0.288)  # mean, median, stddev

    arr = rng.random(size=shape, dtype="f4")
    for byteorder in (">", "<"):
        data = arr.astype(dtype=f"{byteorder}f4", copy=COPY_IF_NEEDED)
        res = sigma_clipped_stats(data, sigma=3, maxiters=5)
        assert_allclose(res, expected, rtol=3e-3)


def test_invalid_sigma_clip():
    """Test sigma_clip of data containing invalid values."""

    data = np.ones((5, 5))
    data[2, 2] = 1000
    data[3, 4] = np.nan
    data[1, 1] = np.inf

    data_ma = np.ma.MaskedArray(data)

    with pytest.warns(AstropyUserWarning, match=r"Input data contains invalid values"):
        result = sigma_clip(data)
    with pytest.warns(AstropyUserWarning, match=r"Input data contains invalid values"):
        result_ma = sigma_clip(data_ma)

    assert_equal(result.data, result_ma.data)
    assert_equal(result.mask, result_ma.mask)

    # Pre #4051 if data contains any NaN or infs sigma_clip returns the
    # mask containing `False` only or TypeError if data also contains a
    # masked value.
    assert result.mask[2, 2]
    assert result.mask[3, 4]
    assert result.mask[1, 1]

    with pytest.warns(AstropyUserWarning, match=r"Input data contains invalid values"):
        result2 = sigma_clip(data, axis=0)
    assert result2.mask[1, 1]
    assert result2.mask[3, 4]

    with pytest.warns(AstropyUserWarning, match=r"Input data contains invalid values"):
        result3 = sigma_clip(data, axis=0, copy=False)
    assert result3.mask[1, 1]
    assert result3.mask[3, 4]

    # stats along axis with all nans
    data[0, :] = np.nan  # row of all nans
    with pytest.warns(AstropyUserWarning, match=r"Input data contains invalid values"):
        _, minarr, maxarr = sigma_clip(data, axis=1, masked=False, return_bounds=True)
    assert np.isnan(minarr[0])
    assert np.isnan(maxarr[0])


def test_sigmaclip_negative_axis():
    """Test that dimensions are expanded correctly even if axis is negative."""
    data = np.ones((3, 4))
    # without correct expand_dims this would raise a ValueError
    sigma_clip(data, axis=-1)


def test_sigmaclip_fully_masked():
    """
    Make sure a fully masked array is returned when sigma clipping a
    fully masked array.
    """
    data = np.ma.MaskedArray(
        data=[[1.0, 0.0], [0.0, 1.0]], mask=[[True, True], [True, True]]
    )
    clipped_data = sigma_clip(data)
    assert np.ma.allequal(data, clipped_data)

    clipped_data = sigma_clip(data, masked=False)
    assert not isinstance(clipped_data, np.ma.MaskedArray)
    assert np.all(np.isnan(clipped_data))

    clipped_data, low, high = sigma_clip(data, return_bounds=True)
    assert np.ma.allequal(data, clipped_data)
    assert np.isnan(low)
    assert np.isnan(high)


def test_sigmaclip_empty_masked():
    """
    Make sure an empty masked array is returned when sigma clipping an
    empty masked array.
    """
    data = np.ma.MaskedArray(data=[], mask=[])
    clipped_data = sigma_clip(data)
    assert np.ma.allequal(data, clipped_data)

    clipped_data, low, high = sigma_clip(data, return_bounds=True)
    assert np.ma.allequal(data, clipped_data)
    assert np.isnan(low)
    assert np.isnan(high)


def test_sigmaclip_empty():
    """
    Make sure an empty array is returned when sigma clipping an empty
    array.
    """
    data = np.array([])
    clipped_data = sigma_clip(data)
    assert isinstance(clipped_data, np.ma.MaskedArray)
    assert_equal(data, clipped_data.data)

    clipped_data, low, high = sigma_clip(data, return_bounds=True)
    assert_equal(data, clipped_data)
    assert np.isnan(low)
    assert np.isnan(high)


def test_sigma_clip_axis_tuple_3D():
    """Test sigma clipping over a subset of axes (issue #7227)."""
    data = np.sin(0.78 * np.arange(27)).reshape(3, 3, 3)
    mask = np.zeros_like(data, dtype=np.bool_)

    data_t = np.rollaxis(data, 1, 0)
    mask_t = np.rollaxis(mask, 1, 0)

    # Loop over what was originally axis 1 and clip each plane directly:
    for data_plane, mask_plane in zip(data_t, mask_t):
        mean = data_plane.mean()
        maxdev = 1.5 * data_plane.std()
        mask_plane[:] = np.logical_or(
            data_plane < mean - maxdev, data_plane > mean + maxdev
        )

    # Do the equivalent thing using sigma_clip:
    result = sigma_clip(data, sigma=1.5, cenfunc=np.mean, maxiters=1, axis=(0, -1))

    assert_equal(result.mask, mask)


def test_sigmaclip_repr():
    sigclip = SigmaClip()

    sigclip_repr = (
        "SigmaClip(sigma=3.0, sigma_lower=3.0, sigma_upper=3.0,"
        " maxiters=5, cenfunc='median', stdfunc='std', "
        "grow=False)"
    )
    sigclip_str = (
        "<SigmaClip>\n    sigma: 3.0\n    sigma_lower: 3.0\n"
        "    sigma_upper: 3.0\n    maxiters: 5\n"
        "    cenfunc: 'median'\n    stdfunc: 'std'\n"
        "    grow: False"
    )

    assert repr(sigclip) == sigclip_repr
    assert str(sigclip) == sigclip_str


def test_sigma_clippped_stats_unit():
    data = np.array([1, 1]) * u.kpc
    result = sigma_clipped_stats(data)
    assert result == (1.0 * u.kpc, 1.0 * u.kpc, 0.0 * u.kpc)


def test_sigma_clippped_stats_all_masked():
    """
    Test sigma_clipped_stats when the input array is completely masked.
    """

    arr = np.ma.MaskedArray(np.arange(10), mask=True)
    result = sigma_clipped_stats(arr)
    assert result == (np.ma.masked, np.ma.masked, np.ma.masked)

    arr = np.ma.MaskedArray(np.zeros(10), mask=False)
    result = sigma_clipped_stats(arr, mask_value=0.0)
    assert result == (np.ma.masked, np.ma.masked, np.ma.masked)

    arr = np.ma.MaskedArray(np.arange(10), mask=False)
    mask = arr < 20
    result = sigma_clipped_stats(arr, mask=mask)
    assert result == (np.ma.masked, np.ma.masked, np.ma.masked)


def test_sigma_clip_masked_data_values():
    """
    Test that the data values & type returned by sigma_clip are the same as
    its input when using masked=True (rather than being upcast to float64 &
    containing NaNs as in issue #10605) and also that the input data get
    copied or referenced as appropriate.
    """

    data = np.array([-2, 5, -5, -6, 20, 14, 1])

    result = sigma_clip(data, sigma=1.5, maxiters=3, axis=None, masked=True, copy=True)

    assert result.dtype == data.dtype
    assert_equal(result.data, data)
    assert not np.shares_memory(result.data, data)

    result = sigma_clip(data, sigma=1.5, maxiters=3, axis=None, masked=True, copy=False)

    assert result.dtype == data.dtype
    assert_equal(result.data, data)
    assert np.shares_memory(result.data, data)
    # (The fact that the arrays share memory probably also means they're the
    # same, but doesn't strictly prove it, eg. one could be reversed.)

    result = sigma_clip(data, sigma=1.5, maxiters=3, axis=0, masked=True, copy=True)

    assert result.dtype == data.dtype
    assert_equal(result.data, data)
    assert not np.shares_memory(result.data, data)

    result = sigma_clip(data, sigma=1.5, maxiters=3, axis=0, masked=True, copy=False)

    assert result.dtype == data.dtype
    assert_equal(result.data, data)
    assert np.shares_memory(result.data, data)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_sigma_clip_grow():
    """
    Test sigma_clip with growth of masking to include the neighbours within a
    specified radius of deviant values.
    """

    # We could use a random seed here, but enumerating the data guarantees that
    # we test sigma_clip itself and not random number generation.

    # fmt: off
    data = np.array(
        [
            -0.2 ,  0.48, -0.52, -0.56,  1.97,  1.39,  0.09,  0.28,  0.77,  1.25,
             1.01, -1.3 ,  0.27,  0.23,  1.35,  0.89, -2.  , -0.37,  1.67, -0.44,
            -0.54,  0.48,  3.25, -1.02, -0.58,  0.12,  0.3 ,  0.52,  0.  ,  1.34,
            -0.71, -0.83, -2.37, -1.86, -0.86,  0.56, -1.27,  0.12, -1.06,  0.33,
            -2.36, -0.2 , -1.54, -0.97, -1.31,  0.29,  0.38, -0.75,  0.33,  1.35,
             0.07,  0.25, -0.01,  1.  ,  1.33, -0.92, -1.55,  0.02,  0.76, -0.66,
             0.86, -0.01,  0.05,  0.67,  0.85, -0.96, -0.02, -2.3 , -0.65, -1.22,
            -1.33,  1.07,  0.72,  0.69,  1.  , -0.5 , -0.62, -0.92, -0.73,  0.22,
             0.05, -1.16,  0.82,  0.43,  1.01,  1.82, -1.  ,  0.85, -0.13,  0.91,
             0.19,  2.17, -0.11,  2.  ,  0.03,  0.8 ,  0.12, -0.75,  0.58,  0.15,
        ]
    )
    # fmt: on

    # Test growth to immediate neighbours in simple 1D case:
    filtered_data = sigma_clip(data, sigma=2, maxiters=3, grow=1)

    # Indices of the 26/100 points expected to be masked:

    # fmt: off
    expected = np.array(
        [
             3,  4,  5, 15, 16, 17, 21, 22, 23, 31, 32, 33, 39,
            40, 41, 66, 67, 68, 84, 85, 86, 90, 91, 92, 93, 94,
        ]
    )
    # fmt: on

    assert np.array_equal(np.where(filtered_data.mask)[0], expected)

    # Test block growth in 2 of 3 dimensions (as in a 2D model set):
    data = data.reshape(4, 5, 5)
    filtered_data = sigma_clip(data, sigma=2.1, maxiters=1, grow=1.5, axis=(1, 2))

    # fmt: off
    expected = np.array(
        [
            [
                0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            ],
            [
                3, 3, 3, 4, 4, 4, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4,
                2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4,
            ],
            [
                1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 0, 1, 2, 3, 0, 1, 0, 1,
                1, 2, 3, 1, 2, 3, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0, 1, 2,
            ],
        ]
    )
    # fmt: on

    assert np.array_equal(np.where(filtered_data.mask), expected)

    # Test ~spherical growth (of a single very-deviant point) in 3D data:
    data[1, 2, 2] = 100.0
    filtered_data = sigma_clip(data, sigma=3.0, maxiters=1, grow=2.0)

    # fmt: off
    expected = np.array(
        [
            [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 3
            ],
            [
                1, 1, 1, 2, 2, 2, 3, 3, 3, 0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 1,
                1, 1, 2, 2, 2, 3, 3, 3, 2
            ],
            [
                1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 1, 2, 3, 0, 1, 2, 3, 4, 1, 2, 3, 2, 1,
                2, 3, 1, 2, 3, 1, 2, 3, 2
            ],
        ]
    )
    # fmt: on

    assert np.array_equal(np.where(filtered_data.mask), expected)


@pytest.mark.parametrize(
    ("axis", "bounds_shape"),
    [
        (0, (4, 5, 6, 7)),
        (1, (3, 5, 6, 7)),
        (-1, (3, 4, 5, 6)),
        ((1, 3), (3, 5, 7)),
        ((3, 1), (3, 5, 7)),
        ((1, 2, 4), (3, 6)),
    ],
)
def test_sigma_clip_axis_shapes(axis, bounds_shape):
    # Check the shapes of the output for different use cases

    with NumpyRNGContext(12345):
        array = np.random.random((3, 4, 5, 6, 7))

    result1 = sigma_clip(array, axis=axis)
    assert result1.shape == array.shape

    result2, bound1, bound2 = sigma_clip(array, axis=axis, return_bounds=True)
    assert result2.shape == array.shape
    assert bound1.shape == bounds_shape
    assert bound2.shape == bounds_shape


@pytest.mark.parametrize(
    "dtype", [">f2", "<f2", ">f4", "<f4", ">f8", "<f8", "<i4", ">i8"]
)
def test_sigma_clip_dtypes(dtype):
    # Compare the outputs for various dtypes
    with NumpyRNGContext(12345):
        array = np.random.randint(-5, 5, 1000).astype(float)
    array[30] = 100
    reference = sigma_clip(array, copy=True, masked=False)
    actual = sigma_clip(array.astype(dtype), copy=True, masked=False)
    assert_equal(reference, actual)

    # Check that the output dtype is the same as the input dtype
    arr = np.ones((10, 10), dtype=np.float32)
    arr2 = sigma_clip(arr, cenfunc=np.mean, axis=0, masked=False)
    assert arr2.dtype == arr.dtype

    # Check that int dtype is converted to float32, not float
    arr = np.ones((10, 10), dtype=int)
    arr2 = sigma_clip(arr, axis=0, masked=False)
    assert arr2.dtype == np.float32


def test_mad_std():
    # Check with a small array where we know how the result should differ from std

    # Choose an array with few elements and a high proportion of outliers since
    # in this case std and mad_std will be very different.
    array = np.array([1, 10000, 4, 3, 10000])

    # First check with regular std, which shouldn't remove any values
    result_std = sigma_clip(
        array, cenfunc="median", stdfunc="std", maxiters=1, sigma=5, masked=False
    )
    assert_equal(result_std, array)

    # Whereas using mad_std should result in the high values being removed
    result_mad_std = sigma_clip(
        array, cenfunc="median", stdfunc="mad_std", maxiters=1, sigma=5, masked=False
    )
    assert_equal(result_mad_std, [1, 4, 3])

    # We now check this again but with the axis= keyword set since at the time
    # of writing this test this relies on a fast C implementation in which we
    # have re-inplemented mad_std.

    result_std = sigma_clip(
        array,
        cenfunc="median",
        stdfunc="std",
        maxiters=1,
        sigma=5,
        masked=False,
        axis=0,
    )
    assert_equal(result_std, array)

    result_mad_std = sigma_clip(
        array,
        cenfunc="median",
        stdfunc="mad_std",
        maxiters=1,
        sigma=5,
        masked=False,
        axis=0,
    )
    assert_equal(result_mad_std, [1, np.nan, 4, 3, np.nan])


def test_mad_std_large():
    # And now test with a larger array and compare with Python mad_std function

    with NumpyRNGContext(12345):
        array = np.random.uniform(-1, 2, (30, 40))

    def nan_mad_std(data, axis=None):
        return mad_std(data, axis=axis, ignore_nan=True)

    result1 = sigma_clip(
        array, sigma=2, maxiters=None, stdfunc=nan_mad_std, axis=0, masked=False
    )

    result2 = sigma_clip(
        array, sigma=2, maxiters=None, stdfunc="mad_std", axis=0, masked=False
    )

    assert_allclose(result1, result2)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_propagation_of_mask():
    # Quick test to check that the mask is propagated correctly
    x = np.array([1, 1, 1, 1, 1, 1, 1, 1, 5, 5]).astype(float)
    y = np.ma.masked_where(x > 1, x)

    assert_allclose(sigma_clipped_stats(y, grow=1), (1, 1, 0))

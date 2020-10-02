# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import assert_equal, assert_allclose

try:
    from scipy import stats  # used in testing
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

from astropy import units as u
from astropy.stats.sigma_clipping import sigma_clip, SigmaClip, sigma_clipped_stats
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
        filtered_data2 = sigma_clip(randvar, sigma=1, maxiters=2,
                                    stdfunc=np.var)
        assert not np.all(filtered_data.mask == filtered_data2.mask)

        filtered_data3 = sigma_clip(randvar, sigma=1, maxiters=2,
                                    cenfunc=np.mean)
        assert not np.all(filtered_data.mask == filtered_data3.mask)

        # make sure the maxiters=None method works at all.
        filtered_data = sigma_clip(randvar, sigma=3, maxiters=None)

        # test copying
        assert filtered_data.data[0] == randvar[0]
        filtered_data.data[0] += 1.
        assert filtered_data.data[0] != randvar[0]

        filtered_data = sigma_clip(randvar, sigma=3, maxiters=None,
                                   copy=False)
        assert filtered_data.data[0] == randvar[0]
        filtered_data.data[0] += 1.
        assert filtered_data.data[0] == randvar[0]

        # test axis
        data = np.arange(5) + np.random.normal(0., 0.05, (5, 5)) + \
            np.diag(np.ones(5))
        filtered_data = sigma_clip(data, axis=0, sigma=2.3)
        assert filtered_data.count() == 20
        filtered_data = sigma_clip(data, axis=1, sigma=2.3)
        assert filtered_data.count() == 25


@pytest.mark.skipif('not HAS_SCIPY')
def test_compare_to_scipy_sigmaclip():
    # need to seed the numpy RNG to make sure we don't get some
    # amazingly flukey random number that breaks one of the tests

    with NumpyRNGContext(12345):

        randvar = np.random.randn(10000)

        astropyres = sigma_clip(randvar, sigma=3, maxiters=None,
                                cenfunc=np.mean)
        scipyres = stats.sigmaclip(randvar, 3, 3)[0]

        assert astropyres.count() == len(scipyres)
        assert_equal(astropyres[~astropyres.mask].data, scipyres)


def test_sigma_clip_scalar_mask():
    """Test that the returned mask is not a scalar."""
    data = np.arange(5)
    result = sigma_clip(data, sigma=100., maxiters=1)
    assert result.mask.shape != ()


def test_sigma_clip_class():
    with NumpyRNGContext(12345):
        data = np.random.randn(100)
        data[10] = 1.e5
        sobj = SigmaClip(sigma=1, maxiters=2)
        sfunc = sigma_clip(data, sigma=1, maxiters=2)
        assert_equal(sobj(data), sfunc)


def test_sigma_clip_mean():
    with NumpyRNGContext(12345):
        data = np.random.normal(0., 0.05, (10, 10))
        data[2, 2] = 1.e5
        sobj1 = SigmaClip(sigma=1, maxiters=2, cenfunc='mean')
        sobj2 = SigmaClip(sigma=1, maxiters=2, cenfunc=np.nanmean)
        assert_equal(sobj1(data), sobj2(data))
        assert_equal(sobj1(data, axis=0), sobj2(data, axis=0))


def test_sigma_clip_invalid_cenfunc_stdfunc():
    with pytest.raises(ValueError):
        SigmaClip(cenfunc='invalid')

    with pytest.raises(ValueError):
        SigmaClip(stdfunc='invalid')


def test_sigma_clipped_stats():
    """Test list data with input mask or mask_value (#3268)."""
    # test list data with mask
    data = [0, 1]
    mask = np.array([True, False])
    result = sigma_clipped_stats(data, mask=mask)
    # Check that the result of np.ma.median was converted to a scalar
    assert isinstance(result[1], float)
    assert result == (1., 1., 0.)

    result2 = sigma_clipped_stats(data, mask=mask, axis=0)
    assert_equal(result, result2)

    # test list data with mask_value
    result = sigma_clipped_stats(data, mask_value=0.)
    assert isinstance(result[1], float)
    assert result == (1., 1., 0.)

    # test without mask
    data = [0, 2]
    result = sigma_clipped_stats(data)
    assert isinstance(result[1], float)
    assert result == (1., 1., 1.)

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
        data[10] = 1.e5
        mean1, median1, stddev1 = sigma_clipped_stats(data)
        mean2, median2, stddev2 = sigma_clipped_stats(data, std_ddof=1)
        assert mean1 == mean2
        assert median1 == median2
        assert_allclose(stddev1, 0.98156805711673156)
        assert_allclose(stddev2, 0.98161731654802831)


def test_invalid_sigma_clip():
    """Test sigma_clip of data containing invalid values."""

    data = np.ones((5, 5))
    data[2, 2] = 1000
    data[3, 4] = np.nan
    data[1, 1] = np.inf

    with pytest.warns(AstropyUserWarning,
                      match=r'Input data contains invalid values'):
        result = sigma_clip(data)

    # Pre #4051 if data contains any NaN or infs sigma_clip returns the
    # mask containing `False` only or TypeError if data also contains a
    # masked value.
    assert result.mask[2, 2]
    assert result.mask[3, 4]
    assert result.mask[1, 1]

    with pytest.warns(AstropyUserWarning,
                      match=r'Input data contains invalid values'):
        result2 = sigma_clip(data, axis=0)
    assert result2.mask[1, 1]
    assert result2.mask[3, 4]

    with pytest.warns(AstropyUserWarning,
                      match=r'Input data contains invalid values'):
        result3 = sigma_clip(data, axis=0, copy=False)
    assert result3.mask[1, 1]
    assert result3.mask[3, 4]

    # stats along axis with all nans
    data[0, :] = np.nan     # row of all nans
    with pytest.warns(AstropyUserWarning,
                      match=r'Input data contains invalid values'):
        result4, minarr, maxarr = sigma_clip(data, axis=1, masked=False,
                                             return_bounds=True)
    assert np.isnan(minarr[0])
    assert np.isnan(maxarr[0])


def test_sigmaclip_negative_axis():
    """Test that dimensions are expanded correctly even if axis is negative."""
    data = np.ones((3, 4))
    # without correct expand_dims this would raise a ValueError
    sigma_clip(data, axis=-1)


def test_sigmaclip_fully_masked():
    """Make sure a fully masked array is returned when sigma clipping a fully
    masked array.
    """

    data = np.ma.MaskedArray(data=[[1., 0.], [0., 1.]],
                             mask=[[True, True], [True, True]])
    clipped_data = sigma_clip(data)
    np.ma.allequal(data, clipped_data)

    clipped_data = sigma_clip(data, masked=False)
    assert not isinstance(clipped_data, np.ma.MaskedArray)
    assert np.all(np.isnan(clipped_data))


def test_sigmaclip_empty_masked():
    """Make sure a empty masked array is returned when sigma clipping an empty
    masked array.
    """

    data = np.ma.MaskedArray(data=[], mask=[])
    clipped_data = sigma_clip(data)
    np.ma.allequal(data, clipped_data)


def test_sigmaclip_empty():
    """Make sure a empty array is returned when sigma clipping an empty array.
    """

    data = np.array([])
    clipped_data = sigma_clip(data)
    assert_equal(data, clipped_data)


def test_sigma_clip_axis_tuple_3D():
    """Test sigma clipping over a subset of axes (issue #7227).
    """

    data = np.sin(0.78 * np.arange(27)).reshape(3, 3, 3)
    mask = np.zeros_like(data, dtype=np.bool_)

    data_t = np.rollaxis(data, 1, 0)
    mask_t = np.rollaxis(mask, 1, 0)

    # Loop over what was originally axis 1 and clip each plane directly:
    for data_plane, mask_plane in zip(data_t, mask_t):

        mean = data_plane.mean()
        maxdev = 1.5 * data_plane.std()
        mask_plane[:] = np.logical_or(data_plane < mean - maxdev,
                                      data_plane > mean + maxdev)

    # Do the equivalent thing using sigma_clip:
    result = sigma_clip(data, sigma=1.5, cenfunc=np.mean, maxiters=1,
                        axis=(0, -1))

    assert_equal(result.mask, mask)


def test_sigmaclip_repr():

    sigclip = SigmaClip()
    median_str = str(sigclip._parse_cenfunc('median'))
    std_str = str(sigclip._parse_stdfunc('std'))

    sigclip_repr = ('SigmaClip(sigma=3.0, sigma_lower=3.0, sigma_upper=3.0,'
                    ' maxiters=5, cenfunc={}, stdfunc={}, '
                    'grow=False)'.format(median_str, std_str))
    sigclip_str = ('<SigmaClip>\n    sigma: 3.0\n    sigma_lower: 3.0\n'
                   '    sigma_upper: 3.0\n    maxiters: 5\n'
                   '    cenfunc: {}\n    stdfunc: {}\n'
                   '    grow: False'.format(median_str, std_str))

    assert repr(sigclip) == sigclip_repr
    assert str(sigclip) == sigclip_str


def test_sigma_clippped_stats_unit():
    data = np.array([1, 1]) * u.kpc
    result = sigma_clipped_stats(data)
    assert result == (1. * u.kpc, 1. * u.kpc, 0. * u.kpc)


def test_sigma_clippped_stats_all_masked():
    """
    Test sigma_clipped_stats when the input array is completely masked.
    """

    arr = np.ma.MaskedArray(np.arange(10), mask=True)
    result = sigma_clipped_stats(arr)
    assert result == (np.ma.masked, np.ma.masked, np.ma.masked)

    arr = np.ma.MaskedArray(np.zeros(10), mask=False)
    result = sigma_clipped_stats(arr, mask_value=0.)
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

    result = sigma_clip(data, sigma=1.5, maxiters=3, axis=None, masked=True,
                        copy=True)

    assert result.dtype == data.dtype
    assert_equal(result.data, data)
    assert not np.shares_memory(result.data, data)

    result = sigma_clip(data, sigma=1.5, maxiters=3, axis=None, masked=True,
                        copy=False)

    assert result.dtype == data.dtype
    assert_equal(result.data, data)
    assert np.shares_memory(result.data, data)
    # (The fact that the arrays share memory probably also means they're the
    # same, but doesn't strictly prove it, eg. one could be reversed.)

    result = sigma_clip(data, sigma=1.5, maxiters=3, axis=0, masked=True,
                        copy=True)

    assert result.dtype == data.dtype
    assert_equal(result.data, data)
    assert not np.shares_memory(result.data, data)

    result = sigma_clip(data, sigma=1.5, maxiters=3, axis=0, masked=True,
                        copy=False)

    assert result.dtype == data.dtype
    assert_equal(result.data, data)
    assert np.shares_memory(result.data, data)


@pytest.mark.skipif('not HAS_SCIPY')
def test_sigma_clip_grow():
    """
    Test sigma_clip with growth of masking to include the neighbours within a
    specified radius of deviant values.
    """

    # We could use a random seed here, but enumerating the data guarantees that
    # we test sigma_clip itself and not random number generation.
    data = np.array(
        [-0.2 ,  0.48, -0.52, -0.56,  1.97,  1.39,  0.09,  0.28,  0.77,  1.25,
          1.01, -1.3 ,  0.27,  0.23,  1.35,  0.89, -2.  , -0.37,  1.67, -0.44,
         -0.54,  0.48,  3.25, -1.02, -0.58,  0.12,  0.3 ,  0.52,  0.  ,  1.34,
         -0.71, -0.83, -2.37, -1.86, -0.86,  0.56, -1.27,  0.12, -1.06,  0.33,
         -2.36, -0.2 , -1.54, -0.97, -1.31,  0.29,  0.38, -0.75,  0.33,  1.35,
          0.07,  0.25, -0.01,  1.  ,  1.33, -0.92, -1.55,  0.02,  0.76, -0.66,
          0.86, -0.01,  0.05,  0.67,  0.85, -0.96, -0.02, -2.3 , -0.65, -1.22,
         -1.33,  1.07,  0.72,  0.69,  1.  , -0.5 , -0.62, -0.92, -0.73,  0.22,
          0.05, -1.16,  0.82,  0.43,  1.01,  1.82, -1.  ,  0.85, -0.13,  0.91,
          0.19,  2.17, -0.11,  2.  ,  0.03,  0.8 ,  0.12, -0.75,  0.58,  0.15]
    )

    # Test growth to immediate neighbours in simple 1D case:
    filtered_data = sigma_clip(data, sigma=2, maxiters=3, grow=1)

    # Indices of the 26/100 points expected to be masked:
    expected = np.array([3, 4, 5, 15, 16, 17, 21, 22, 23, 31, 32, 33, 39, 40,
                         41, 66, 67, 68, 84, 85, 86, 90, 91, 92, 93, 94])

    assert np.array_equal(np.where(filtered_data.mask)[0], expected)

    # Test block growth in 2 of 3 dimensions (as in a 2D model set):
    data = data.reshape(4, 5, 5)
    filtered_data = sigma_clip(data, sigma=2.1, maxiters=1, grow=1.5,
                               axis=(1, 2))
    expected = np.array(
        [[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3],
         [3, 3, 3, 4, 4, 4, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4,
          2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, 4, 4],
         [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 0, 1, 2, 3, 0, 1, 0, 1,
          1, 2, 3, 1, 2, 3, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0, 1, 2]]
    )

    assert np.array_equal(np.where(filtered_data.mask), expected)

    # Test ~spherical growth (of a single very-deviant point) in 3D data:
    data[1, 2, 2] = 100.
    filtered_data = sigma_clip(data, sigma=3., maxiters=1, grow=2.)

    expected = np.array(
        [[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
          2, 2, 2, 2, 2, 2, 2, 2, 3],
         [1, 1, 1, 2, 2, 2, 3, 3, 3, 0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 1,
          1, 1, 2, 2, 2, 3, 3, 3, 2],
         [1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 1, 2, 3, 0, 1, 2, 3, 4, 1, 2, 3, 2, 1,
          2, 3, 1, 2, 3, 1, 2, 3, 2]]
    )

    assert np.array_equal(np.where(filtered_data.mask), expected)

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

from ..sigma_clipping import sigma_clip, SigmaClip, sigma_clipped_stats
from ...utils.misc import NumpyRNGContext


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

    result = sigma_clip(data)

    # Pre #4051 if data contains any NaN or infs sigma_clip returns the
    # mask containing `False` only or TypeError if data also contains a
    # masked value.
    assert result.mask[2, 2]
    assert result.mask[3, 4]
    assert result.mask[1, 1]

    result2 = sigma_clip(data, axis=0)
    assert result2.mask[1, 1]
    assert result2.mask[3, 4]

    result3 = sigma_clip(data, axis=0, copy=False)
    assert result3.mask[1, 1]
    assert result3.mask[3, 4]

    # stats along axis with all nans
    data[0, :] = np.nan     # row of all nans
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
    mask = np.zeros_like(data, dtype=np.bool)

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
    sigclip_repr = ('SigmaClip(sigma=3.0, sigma_lower=3.0, sigma_upper=3.0,'
                    ' maxiters=5, cenfunc=')
    sigclip_str = ('<SigmaClip>\n    sigma: 3.0\n    sigma_lower: 3.0\n'
                   '    sigma_upper: 3.0\n    maxiters: 5\n    cenfunc: ')

    assert repr(sigclip).startswith(sigclip_repr)
    assert str(sigclip).startswith(sigclip_str)

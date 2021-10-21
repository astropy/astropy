# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys

import pytest
import numpy as np
from numpy.testing import assert_equal

from astropy import units as u
from astropy.time import Time
from astropy.utils.exceptions import AstropyUserWarning

from astropy.timeseries.sampled import TimeSeries
from astropy.timeseries.downsample import aggregate_downsample, reduceat


INPUT_TIME = Time(['2016-03-22T12:30:31', '2016-03-22T12:30:32',
                   '2016-03-22T12:30:33', '2016-03-22T12:30:34',
                   '2016-03-22T12:30:35'])


def test_reduceat():
    add_output = np.add.reduceat(np.arange(8),[0, 4, 1, 5, 2, 6, 3, 7])
    # Similar to np.add for an array input.
    sum_output = reduceat(np.arange(8), [0, 4, 1, 5, 2, 6, 3, 7], np.sum)
    assert_equal(sum_output, add_output)

    mean_output = reduceat(np.arange(8), np.arange(8)[::2], np.mean)
    assert_equal(mean_output, np.array([0.5, 2.5, 4.5, 6.5]))
    nanmean_output = reduceat(np.arange(8), [0, 4, 1, 5, 2, 6, 3, 7], np.mean)
    assert_equal(nanmean_output, np.array([1.5, 4, 2.5, 5, 3.5, 6, 4.5, 7.]))
    assert_equal(reduceat(np.arange(8), np.arange(8)[::2], np.mean),
                 reduceat(np.arange(8), np.arange(8)[::2], np.nanmean))


def test_timeseries_invalid():
    with pytest.raises(TypeError) as exc:
        aggregate_downsample(None)
    assert exc.value.args[0] == ("time_series should be a TimeSeries")


def test_time_bin_invalid():

    # Make sure to raise the right exception when time_bin_* is passed incorrectly.

    with pytest.raises(TypeError, match=r"'time_bin_size' should be a Quantity or a TimeDelta"):
        aggregate_downsample(TimeSeries(), time_bin_size=1)


def test_binning_arg_invalid():
    ts = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4, 5]], names=['a'])
    with pytest.raises(TypeError, match=r"With single 'time_bin_start' either 'n_bins', "
                                         "'time_bin_size' or time_bin_end' must be provided"):
        aggregate_downsample(ts)


def test_time_bin_conversion():

    ts = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4, 5]], names=['a'])

    # Make sure time_bin_start and time_bin_end are properly converted to Time
    down_start = aggregate_downsample(ts, time_bin_start=['2016-03-22T12:30:31'],
                                      time_bin_size=[1]*u.s)
    assert_equal(down_start.time_bin_start.isot, ['2016-03-22T12:30:31.000'])

    down_end = aggregate_downsample(ts, time_bin_start=['2016-03-22T12:30:31', '2016-03-22T12:30:33'],
                                    time_bin_end='2016-03-22T12:30:34')
    assert_equal(down_end.time_bin_end.isot, ['2016-03-22T12:30:33.000', '2016-03-22T12:30:34.000'])


def test_time_bin_end_auto():

    ts = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4, 5]], names=['a'])

    # Interpret `time_bin_end` as the end of timeseries when `time_bin_start` is
    # an array and `time_bin_size` is not provided
    down_auto_end = aggregate_downsample(ts, time_bin_start=['2016-03-22T12:30:31', '2016-03-22T12:30:33'])
    assert_equal(down_auto_end.time_bin_end.isot, ['2016-03-22T12:30:33.000', '2016-03-22T12:30:35.000'])


def test_time_bin_start_array():

    ts = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4, 5]], names=['a'])

    # When `time_bin_end` is an array and `time_bin_start` is not provided, `time_bin_start` is converted
    # to an array with its first element set to the start of the timeseries and rest populated using
    #`time_bin_end`. This case is separately tested since `BinnedTimeSeries` allows `time_bin_end` to
    # be an array only if `time_bin_start` is an array.
    down_start_array = aggregate_downsample(ts, time_bin_end=['2016-03-22T12:30:33', '2016-03-22T12:30:35'])
    assert_equal(down_start_array.time_bin_start.isot, ['2016-03-22T12:30:31.000', '2016-03-22T12:30:33.000'])


def test_nbins():

    ts = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4, 5]], names=['a'])

    # n_bins should default to the number needed to fit all the original points
    down_nbins = aggregate_downsample(ts, n_bins=2)
    assert_equal(down_nbins.time_bin_start.isot, ['2016-03-22T12:30:31.000', '2016-03-22T12:30:33.000'])


def test_downsample():
    ts = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4, 5]], names=['a'])
    ts_units = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4, 5] * u.count], names=['a'])

    # Avoid precision problems with floating-point comparisons on 32bit
    if sys.maxsize > 2**32:
        # 64 bit
        time_bin_incr = 1 * u.s
        time_bin_start = None
    else:
        # 32 bit
        time_bin_incr = (1 - 1e-6) * u.s
        time_bin_start = ts.time[0] - 1 * u.ns

    down_1 = aggregate_downsample(ts, time_bin_size=time_bin_incr, time_bin_start=time_bin_start)
    u.isclose(down_1.time_bin_size, [1, 1, 1, 1, 1]*time_bin_incr)
    assert_equal(down_1.time_bin_start.isot, Time(['2016-03-22T12:30:31.000', '2016-03-22T12:30:32.000',
                                                   '2016-03-22T12:30:33.000', '2016-03-22T12:30:34.000',
                                                   '2016-03-22T12:30:35.000']))
    assert_equal(down_1["a"].data.data, np.array([1, 2, 3, 4, 5]))

    down_2 = aggregate_downsample(ts, time_bin_size=2*time_bin_incr, time_bin_start=time_bin_start)
    u.isclose(down_2.time_bin_size, [2, 2, 2]*time_bin_incr)
    assert_equal(down_2.time_bin_start.isot, Time(['2016-03-22T12:30:31.000', '2016-03-22T12:30:33.000',
                                                   '2016-03-22T12:30:35.000']))
    assert_equal(down_2["a"].data.data, np.array([1, 3, 5]))

    down_3 = aggregate_downsample(ts, time_bin_size=3*time_bin_incr, time_bin_start=time_bin_start)
    u.isclose(down_3.time_bin_size, [3, 3]*time_bin_incr)
    assert_equal(down_3.time_bin_start.isot, Time(['2016-03-22T12:30:31.000', '2016-03-22T12:30:34.000']))
    assert_equal(down_3["a"].data.data, np.array([2, 4]))

    down_4 = aggregate_downsample(ts, time_bin_size=4*time_bin_incr, time_bin_start=time_bin_start)
    u.isclose(down_4.time_bin_size, [4, 4]*time_bin_incr)
    assert_equal(down_4.time_bin_start.isot, Time(['2016-03-22T12:30:31.000', '2016-03-22T12:30:35.000']))
    assert_equal(down_4["a"].data.data, np.array([2, 5]))

    down_units = aggregate_downsample(ts_units, time_bin_size=4*time_bin_incr, time_bin_start=time_bin_start)
    u.isclose(down_units.time_bin_size, [4, 4]*time_bin_incr)
    assert_equal(down_units.time_bin_start.isot, Time(['2016-03-22T12:30:31.000', '2016-03-22T12:30:35.000']))
    assert down_units["a"].unit.name == 'ct'
    assert_equal(down_units["a"].data, np.array([2.5, 5.0]))

    # Contiguous bins with uneven bin sizes: `time_bin_size` is an array
    down_uneven_bins = aggregate_downsample(ts, time_bin_size=[2, 1, 1]*time_bin_incr,
                                            time_bin_start=time_bin_start)
    u.isclose(down_uneven_bins.time_bin_size, [2, 1, 1]*time_bin_incr)
    assert_equal(down_uneven_bins.time_bin_start.isot, Time(['2016-03-22T12:30:31.000',
                                                             '2016-03-22T12:30:33.000',
                                                             '2016-03-22T12:30:34.000']))
    assert_equal(down_uneven_bins["a"].data.data, np.array([1, 3, 4]))

    # Uncontiguous bins with even bin sizes: `time_bin_start` and `time_bin_end` are both arrays
    down_time_array = aggregate_downsample(ts, time_bin_start=Time(['2016-03-22T12:30:31.000',
                                                                    '2016-03-22T12:30:34.000']),
                                           time_bin_end=Time(['2016-03-22T12:30:32.000',
                                                              '2016-03-22T12:30:35.000']))
    u.isclose(down_time_array.time_bin_size, [1, 1]*u.second)
    assert_equal(down_time_array.time_bin_start.isot, Time(['2016-03-22T12:30:31.000',
                                                            '2016-03-22T12:30:34.000']))
    assert_equal(down_time_array["a"].data.data, np.array([1, 4]))

    # Overlapping bins
    with pytest.warns(AstropyUserWarning, match="Overlapping bins should be avoided since they "
                                                "can lead to double-counting of data during binning."):
        down_overlap_bins = aggregate_downsample(ts, time_bin_start=Time(['2016-03-22T12:30:31.000',
                                                                          '2016-03-22T12:30:33.000']),
                                                 time_bin_end=Time(['2016-03-22T12:30:34',
                                                              '2016-03-22T12:30:36.000']))
        assert_equal(down_overlap_bins["a"].data, np.array([2, 5]))

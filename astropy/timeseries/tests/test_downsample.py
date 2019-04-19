# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import assert_equal

from astropy import units as u
from astropy.time import Time

from astropy.timeseries.sampled import TimeSeries
from astropy.timeseries.downsample import aggregate_downsample, reduceat

INPUT_TIME = Time(['2016-03-22T12:30:31', '2016-03-22T12:30:32',
                   '2016-03-22T12:30:33', '2016-03-22T12:30:34'])
ts = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4]], names=['a'])
ts_units = TimeSeries(time=INPUT_TIME, data=[[1, 2, 3, 4] * u.count], names=['a'])


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

    with pytest.raises(TypeError) as exc:
        aggregate_downsample(TimeSeries())
    assert exc.value.args[0] == ("time_bin_size should be a astropy.unit quantity")


def test_downsample():
    down_1 = aggregate_downsample(ts, time_bin_size=1*u.second)
    u.isclose(down_1.time_bin_size, [1, 1, 1, 1]*u.second)
    assert_equal(down_1.time_bin_start.isot, Time(['2016-03-22T12:30:31.000', '2016-03-22T12:30:32.000',
                                                   '2016-03-22T12:30:33.000', '2016-03-22T12:30:34.000']))
    assert_equal(down_1["a"].data, np.array([1, 2, 3, 4]))

    down_2 = aggregate_downsample(ts, time_bin_size=2*u.second)
    u.isclose(down_2.time_bin_size, [2, 2]*u.second)
    assert_equal(down_2.time_bin_start.isot, Time(['2016-03-22T12:30:31.000', '2016-03-22T12:30:33.000']))
    assert_equal(down_2["a"].data, np.array([1, 3]))

    down_3 = aggregate_downsample(ts, time_bin_size=3*u.second)
    u.isclose(down_3.time_bin_size, [3, 3]*u.second)
    assert_equal(down_3.time_bin_start.isot, Time(['2016-03-22T12:30:31.000', '2016-03-22T12:30:34.000']))
    assert_equal(down_3["a"].data, np.array([2, 4]))

    down_4 = aggregate_downsample(ts, time_bin_size=4*u.second)
    u.isclose(down_4.time_bin_size, [4]*u.second)
    assert_equal(down_4.time_bin_start.isot, Time(['2016-03-22T12:30:31.000']))
    assert_equal(down_4["a"].data, np.array([2]))

    down_units = aggregate_downsample(ts_units, time_bin_size=4*u.second)
    u.isclose(down_units.time_bin_size, [4]*u.second)
    assert_equal(down_units.time_bin_start.isot, Time(['2016-03-22T12:30:31.000']))
    assert down_units["a"].unit.name == 'ct'
    assert_equal(down_units["a"].data, np.array([2.5]))

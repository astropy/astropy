# Licensed under a 3-clause BSD style license - see LICENSE.rst

from datetime import datetime

import pytest
from numpy.testing import assert_equal

from ... import units as u
from ...table import Table, QTable, Column, hstack, vstack
from ...time import Time, TimeDelta
from ..timeseries import SampledTimeSeries, BinnedTimeSeries


INPUT_TIME = Time(['2016-03-22T12:30:31', '2015-01-21T12:30:32', '2016-03-22T12:30:40'])
PLAIN_TABLE = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a', 'b', 'c'])


# Tests of SampledTimeSeries class

def test_sampled_empty_initialization():
    with pytest.raises(TypeError) as exc:
        SampledTimeSeries()
    assert exc.value.args[0] == "'time' has not been specified"


def test_sampled_initialize_only_time():
    ts = SampledTimeSeries(time=INPUT_TIME)
    assert ts['time'] is ts.time
    # NOTE: the object in the table is a copy
    assert_equal(ts.time.isot, INPUT_TIME.isot)


def test_sampled_initialization_with_data():
    ts = SampledTimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert_equal(ts.time.isot, INPUT_TIME.isot)
    assert_equal(ts['a'], [10, 2, 3])
    assert_equal(ts['b'], [4, 5, 6])


def test_sampled_initialization_with_table():
    ts = SampledTimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    assert ts.colnames == ['time', 'a', 'b', 'c']


def test_sampled_initialization_with_deltatime():
    ts = SampledTimeSeries(time=datetime(2018, 7, 1, 10, 10, 10),
                           time_delta=TimeDelta(3, format='sec'),
                           data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert_equal(ts.time.isot, ['2018-07-01T10:10:10.000', '2018-07-01T10:10:13.000', '2018-07-01T10:10:16.000'])


def test_sampled_initialization_with_time_in_data():

    data = PLAIN_TABLE.copy()
    data['time'] = INPUT_TIME

    ts1 = SampledTimeSeries(data=data)

    assert set(ts1.colnames) == set(['time', 'a', 'b', 'c'])
    assert all(ts1.time == INPUT_TIME)

    ts2 = SampledTimeSeries(data=[[10, 2, 3], INPUT_TIME], names=['a', 'time'])
    assert set(ts2.colnames) == set(['time', 'a'])
    assert all(ts2.time == INPUT_TIME)

    with pytest.raises(TypeError) as exc:
        # Don't allow ambiguous cases of passing multiple 'time' columns
        SampledTimeSeries(data=data, time=INPUT_TIME)
    assert exc.value.args[0] == "'time' has been given both in the table and as a keyword argument"

    with pytest.raises(TypeError) as exc:
        # 'time' is a protected name, don't allow ambiguous cases
        SampledTimeSeries(time=INPUT_TIME, data=[[10, 2, 3], INPUT_TIME], names=['a', 'time'])
    assert exc.value.args[0] == "'time' has been given both in the table and as a keyword argument"

# Tests of BinnedTimeSeries class

def test_binned_empty_initialization():
    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries()
    assert exc.value.args[0] == "'start_time' has not been specified"


def test_binned_even_contiguous():
    # Initialize a ``BinnedTimeSeries`` with even contiguous bins by specifying the bin width:
    ts = BinnedTimeSeries(start_time='2016-03-22T12:30:31', bin_size=3 * u.s, data=[[1, 4, 3]])
    assert_equal(ts.start_time.isot, ['2016-03-22T12:30:31.000', '2016-03-22T12:30:34.000', '2016-03-22T12:30:37.000'])
    assert_equal(ts.end_time.isot, ['2016-03-22T12:30:34.000', '2016-03-22T12:30:37.000', '2016-03-22T12:30:40.000'])


def test_binned_uneven_contiguous():
    # Initialize a ``BinnedTimeSeries`` with uneven contiguous bins by giving an end time:
    ts = BinnedTimeSeries(start_time=['2016-03-22T12:30:31', '2016-03-22T12:30:32', '2016-03-22T12:30:40'],
                          end_time='2016-03-22T12:30:55',
                          data=[[1, 4, 3]])
    assert_equal(ts.start_time.isot, ['2016-03-22T12:30:31.000', '2016-03-22T12:30:32.000', '2016-03-22T12:30:40.000'])
    assert_equal(ts.end_time.isot, ['2016-03-22T12:30:32.000', '2016-03-22T12:30:40.000', '2016-03-22T12:30:55.000'])


def test_binned_uneven_non_contiguous():
    # Initialize a ``BinnedTimeSeries`` with uneven non-contiguous bins with lists of start times, bin sizes and data:
    ts = BinnedTimeSeries(start_time=['2016-03-22T12:30:31', '2016-03-22T12:30:38', '2016-03-22T12:34:40'],
                          bin_size=[5, 100, 2]*u.s,
                          data=[[1, 4, 3]])
    assert_equal(ts.start_time.isot, ['2016-03-22T12:30:31.000', '2016-03-22T12:30:38.000', '2016-03-22T12:34:40.000'])
    assert_equal(ts.end_time.isot, ['2016-03-22T12:30:36.000', '2016-03-22T12:32:18.000', '2016-03-22T12:34:42.000'])


def test_binned_uneven_non_contiguous_full():
    # Initialize a ``BinnedTimeSeries`` with uneven non-contiguous bins by specifying the start and end times for the bins:
    ts = BinnedTimeSeries(start_time=['2016-03-22T12:30:31', '2016-03-22T12:30:33', '2016-03-22T12:30:40'],
                          end_time=['2016-03-22T12:30:32', '2016-03-22T12:30:35', '2016-03-22T12:30:41'],
                          data=[[1, 4, 3]])
    assert_equal(ts.start_time.isot, ['2016-03-22T12:30:31.000', '2016-03-22T12:30:33.000', '2016-03-22T12:30:40.000'])
    assert_equal(ts.end_time.isot, ['2016-03-22T12:30:32.000', '2016-03-22T12:30:35.000', '2016-03-22T12:30:41.000'])


# Tests that apply to both

class CommonTimeSeriesTests:

    @pytest.mark.xfail
    def test_stacking(self):
        ts = vstack([self.series, self.series])
        assert isinstance(ts, self.series.__class__)

    @pytest.mark.xfail
    def test_row_slicing(self):
        ts = self.series[:2]
        assert isinstance(ts, self.series.__class__)

    def test_row_indexing(self):
        self.series[0][self.time_attr] == Time('2015-01-21T12:30:32')
        self.series[self.time_attr][0] == Time('2015-01-21T12:30:32')

    @pytest.mark.xfail
    def test_column_slicing(self):
        ts = self.series[self.time_attr, 'a']
        assert isinstance(ts, self.series.__class__)

    def test_column_indexing(self):
        assert_equal(self.series['a'], [1, 2, 11])

    @pytest.mark.xfail
    def test_column_slicing_notime(self):
        tab = self.series['a', 'b']
        assert not isinstance(tab, self.series.__class__)
        assert isinstance(tab, QTable)

    def test_add_column(self):
        self.series['d'] = [1, 2, 3]


class TestSampledTimeSeries(CommonTimeSeriesTests):
    def setup_method(self, method):
        self.series = SampledTimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
        self.time_attr = 'time'


class TestBinnedTimeSeries(CommonTimeSeriesTests):
    def setup_method(self, method):
        self.series = BinnedTimeSeries(start_time=INPUT_TIME, bin_size=3 * u.s, data=PLAIN_TABLE)
        self.time_attr = 'start_time'


# Other tests that need to be tidied up/fixed


@pytest.mark.xfail
def test_sampled_initialization_with_timeseries():
    ts = SampledTimeSeries(data=PLAIN_TABLE, time=INPUT_TIME)
    table = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a2', 'b2', 'c2'])

    # TODO: decide which initializations to allow
    ts2 = SampledTimeSeries([ts, table])



# TODO remove redundancy for input table initialization, or just simply structure better


@pytest.mark.xfail
def test_sampled_sorting():
    # TODO: decide which sorting we would like to allow

    sorted_input_time = INPUT_TIME.sort()
    input_data = Table([[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    ts = TimeSeries(time=INPUT_TIME, data=input_data)
    ts_sort = TimeSeries(time=sorted_input_time, data=input_data)

    # There is automated sorting by time
    assert ts['a'][0] != input_data['a'][0]
    assert all(ts == ts_sort)

    ts.sort()
    ts.sort('a')
    ts.sort('time')


def test_sampled_adding_index_column():
    # TODO: decide which sorting we would like to allow, if it's only by time,
    # the indices should already been set 'time'

    ts = SampledTimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    ts.add_index('a')
    ts.add_index('time')

    # TODO: add asserts

# Licensed under a 3-clause BSD style license - see LICENSE.rst

from datetime import datetime

import pytest
from numpy.testing import assert_equal

from ... import units as u
from ...table import Table, QTable, Column
from ...time import Time, TimeDelta
from ..timeseries import TimeSeries, BinnedTimeSeries


INPUT_TIME = Time(['2016-03-22T12:30:31', '2015-01-21T12:30:32', '2016-03-22T12:30:40'])
PLAIN_TABLE = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a', 'b', 'c'])


# Tests of SampledTimeSeries class

def test_empty_initialization():
    with pytest.raises(ValueError) as exc:
        TimeSeries()
    assert exc.value.args[0] == "Expected a 'time' column in the time series"


def test_initialize_only_time():
    ts = TimeSeries(time=INPUT_TIME)
    assert ts['time'] is ts.time
    assert ts['time'] is INPUT_TIME


def test_initialization_with_data():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert ts['time'] is ts.time
    assert ts['time'] is INPUT_TIME


def test_initialization_with_table():
    ts = TimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    assert ts.colnames == ['time', 'a', 'b', 'c']


def test_initialization_with_deltatime():
    date = datetime(2018, 7, 1, 10, 10, 10)

    input_time = Time(date) + TimeDelta(360, format='sec') * range(3)
    ts = TimeSeries(time=input_time, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])

    # TODO: we probably want to allow to avoid constructing the time array above and
    #  initialize with the start time and step size


def test_initialization_with_time_in_data():
    data = PLAIN_TABLE.copy()
    data['time'] = INPUT_TIME

    ts = TimeSeries(data=data)

    assert isinstance(ts, TimeSeries)
    assert isinstance(ts, QTable)
    assert set(ts.colnames) == set(['a', 'b', 'c', 'time'])
    assert all(ts.time == INPUT_TIME)

    ts2 = TimeSeries(data=[[10, 2, 3], INPUT_TIME], names=['a', 'time'])

    with pytest.raises(ValueError):
        # Don't allow ambiguous cases of passing multiple 'time' columns
        TimeSeries(data=data, time=INPUT_TIME)

    with pytest.raises(ValueError):
        # 'time' is a protected name, don't allow ambiguous cases
        TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], INPUT_TIME], names=['a', 'time'])


@pytest.mark.xfail
def test_initialization_with_timeseries():
    ts = TimeSeries(data=PLAIN_TABLE, time=INPUT_TIME)
    table = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a2', 'b2', 'c2'])

    # TODO: decide which initializations to allow
    ts2 = TimeSeries([ts, table])



# TODO remove redundancy for input table initialization, or just simply structure better


def test_adding_more_columns():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    ts['c'] = Column([5, 4, 3])

    ts.add_column(Column([6, 5, 4], name='d'))


def test_access_column():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])

    # TODO update this to TimeSeriesColumn if we end up adding extra features to them
    assert isinstance(ts['a'], Column)
    # assert all(ts['a'].time == INPUT_TIME)  # not sure whether we need this

    assert ts['a'].name == 'a'


def test_access_time():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])

    assert all(ts['time'] == INPUT_TIME)


@pytest.mark.xfail
def test_access_row():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 20, 3], [4, 5, 6]], names=['a', 'b'])
    assert isinstance(ts[0], TimeSeriesRow)
    # TODO more content checking


@pytest.mark.xfail
def test_access_col_value():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 20, 3], [4, 5, 6]], names=['a', 'b'])
    # TODO: we need to get the sorting by time to get this right
    assert ts['a'][0] == 20


@pytest.mark.xfail
def test_access_time_value():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 20, 3], [4, 5, 6]], names=['a', 'b'])

    # assert ts['a'].time[0] == ts['time'][0]  # not sure whether we need this

    # TODO: we need to get the sorting by time to get this right
    assert ts['time'][0] == Time('2015-01-21T12:30:32')


def test_normal_Columns():
    ts = TimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    assert all(ts.columns['a'] == PLAIN_TABLE['a'])
    assert all(ts['a'][()] == PLAIN_TABLE['a'])


@pytest.mark.xfail
def test_operations():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    new_column = Column([5, 4, 3], name='c')

    tstack = hstack([ts, new_column])

    assert isinstance(tstack, TimeSeries)

    # TODO: do join and vstack, too (but need to decide desired Row behaviour first


@pytest.mark.xfail
def test_sorting():
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


@pytest.mark.xfail
def test_adding_index_column():
    # TODO: decide which sorting we would like to allow, if it's only by time,
    # the indices should already been set 'time'

    ts = TimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    ts.add_index('a')
    ts.add_index('time')

    # TODO: add asserts


@pytest.mark.xfail
def test_access_multiple_columns():
    ts = TimeSeries(time=INPUT_TIME, data=[[1, 20, 3], [4, 5, 6], [4, 3, 2]], names=['a', 'b', 'c'])
    t_out = QTable(data=[[1, 20, 3], [4, 5, 6]], names=['a', 'b'])
    ts_out = TimeSeries(time=INPUT_TIME, data=[[4, 5, 6], ], names=['b', ])

    # TODO: For this we need to be able to do TimeSeries([TimeSeries, TimeSeries, Time])
    # type initialization, or override TimeSeries.__getitem__

    t = ts['a', 'b']
    assert not isinstance(t, TimeSeries)
    assert isinstance(t, QTable)
    assert all(t == t_out)

    t = ts['time', 'b']
    assert isinstance(t, TimeSeries)
    assert all(t == ts_out)


# Tests of BinnedTimeSeries class


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

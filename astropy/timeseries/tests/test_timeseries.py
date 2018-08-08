# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
from datetime import datetime

from ...table import Table, QTable, Column
from ..timeseries import TimeSeries, TimeSeriesColumn, TimeSeriesRow
from ...time import Time, TimeDelta


INPUT_TIME = Time(['2016-03-22T12:30:31', '2015-01-21T12:30:32', '2016-03-22T12:30:40'])
INPUT_DATA = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a', 'b', 'c'])


def test_empty_initialization():
    ts = TimeSeries()
    assert ts.time is None

    assert isinstance(ts, TimeSeries)
    assert isinstance(ts, QTable)


def test_initialize_only_time():
    ts = TimeSeries(time=INPUT_TIME)

    assert isinstance(ts, TimeSeries)
    assert isinstance(ts, QTable)


def test_initialization():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])

    assert isinstance(ts, TimeSeries)
    assert isinstance(ts, QTable)

    assert all(ts.time == INPUT_TIME)

    ts = TimeSeries(time=INPUT_TIME, data=INPUT_DATA)


def test_initialization_with_deltatime():
    date = datetime(2018, 7, 1, 10, 10, 10)

    input_time = Time(date) + TimeDelta(360, format='sec') * range(3)
    ts = TimeSeries(time=input_time, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])

    # TODO: we probably want to allow to avoid constructing the time array above and
    #  initialize with the start time and step size


@pytest.mark.xfail
def test_initialization_with_time_in_data():
    data = INPUT_DATA.copy()
    data['time'] = INPUT_TIME

    ts = TimeSeries(data=data)

    assert isinstance(ts, TimeSeries)
    assert isinstance(ts, QTable)
    assert set(ts.colnames) == set(['a', 'b', 'c', 'time'])
    assert all(ts.time == INPUT_TIME)

    # TODO: add more valid cases when a time is passed in data & names
    ts2 = TimeSeries(data=[[10, 2, 3], INPUT_TIME], names=['a', 'time'])


    with pytest.raises.ValueError:
        # Don't allow ambiguous cases of passing multiple 'time' columns
        TimeSeries(data=data, time=INPUT_TIME)

    with pytest.raises.ValueError:
        # 'time' is a protected name, don't allow ambiguous cases
        TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], INPUT_TIME], names=['a', 'time'])



# TODO remove redundancy for input table initialization, or just simply structure better


def test_adding_more_columns():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    ts['c'] = Column([5, 4, 3])

    ts.add_column(Column([6, 5, 4], name='d'))


@pytest.mark.xfail
def test_access_column():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert isinstance(ts['a'], TimeSeriesColumn)
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
    ts = TimeSeries(time=INPUT_TIME, data=INPUT_DATA)

    assert all(ts.columns['a'] == INPUT_DATA['a'])

    assert all(ts['a'][()] == INPUT_DATA['a'])


@pytest.mark.xfail
def test_operations():
    # TODO: hstack, vstack, join
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    new_column = Column([5, 4, 3], name='c')

    hstack([ts, new_column])


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

    ts = TimeSeries(time=INPUT_TIME, data=INPUT_DATA)
    ts.add_index('a')
    ts.add_index('time')

    # TODO: add asserts


@pytest.mark.xfail
def test_access_multiple_columns():
    ts = TimeSeries(time=INPUT_TIME, data=[[1, 20, 3], [4, 5, 6], [4, 3, 2]], names=['a', 'b', 'c'])
    ts_out = TimeSeries(time=INPUT_TIME, data=[[1, 20, 3], [4, 5, 6]], names=['a', 'b'])

    # TODO: For this we need to be able to do TimeSeries([TimeSeries, TimeSeries, Time])
    # type initialization, or override TimeSeries.__getitem__

    t = ts['a', 'b']
    assert isinstance(t, TimeSeries)
    assert t == ts_out

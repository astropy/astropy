# Licensed under a 3-clause BSD style license - see LICENSE.rst

from datetime import datetime

import pytest
from numpy.testing import assert_equal

from ...table import Table
from ...time import Time, TimeDelta
from ..sampled import SampledTimeSeries


INPUT_TIME = Time(['2016-03-22T12:30:31', '2015-01-21T12:30:32', '2016-03-22T12:30:40'])
PLAIN_TABLE = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a', 'b', 'c'])


# Tests of SampledTimeSeries class

def test_empty_initialization():
    ts = SampledTimeSeries()
    ts['time'] = Time([1, 2, 3], format='mjd')


def test_empty_initialization_invalid():
    # Make sure things crash when the first column added is not a time column
    ts = SampledTimeSeries()
    with pytest.raises(ValueError) as exc:
        ts['flux'] = [1, 2, 3]
    assert exc.value.args[0] == ("SampledTimeSeries requires a column called "
                                 "'time' to be set before data can be added")


def test_initialize_only_time():
    ts = SampledTimeSeries(time=INPUT_TIME)
    assert ts['time'] is ts.time
    # NOTE: the object in the table is a copy
    assert_equal(ts.time.isot, INPUT_TIME.isot)


def test_initialization_with_data():
    ts = SampledTimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert_equal(ts.time.isot, INPUT_TIME.isot)
    assert_equal(ts['a'], [10, 2, 3])
    assert_equal(ts['b'], [4, 5, 6])


def test_initialization_with_table():
    ts = SampledTimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    assert ts.colnames == ['time', 'a', 'b', 'c']


def test_initialization_with_deltatime():
    ts = SampledTimeSeries(time=datetime(2018, 7, 1, 10, 10, 10),
                           time_delta=TimeDelta(3, format='sec'),
                           data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert_equal(ts.time.isot, ['2018-07-01T10:10:10.000', '2018-07-01T10:10:13.000', '2018-07-01T10:10:16.000'])


def test_initialization_with_time_in_data():

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

# Tests that apply to both


# Other tests that need to be tidied up/fixed


@pytest.mark.xfail
def test_initialization_with_timeseries():
    ts = SampledTimeSeries(data=PLAIN_TABLE, time=INPUT_TIME)
    table = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a2', 'b2', 'c2'])

    # TODO: decide which initializations to allow
    ts2 = SampledTimeSeries([ts, table])



# TODO remove redundancy for input table initialization, or just simply structure better


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


def test_adding_index_column():
    # TODO: decide which sorting we would like to allow, if it's only by time,
    # the indices should already been set 'time'

    ts = SampledTimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    ts.add_index('a')
    ts.add_index('time')

    # TODO: add asserts

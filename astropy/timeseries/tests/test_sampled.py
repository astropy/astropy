# Licensed under a 3-clause BSD style license - see LICENSE.rst

from datetime import datetime

import pytest

from numpy.testing import assert_equal, assert_allclose

from astropy.table import Table
from astropy.time import Time, TimeDelta
from astropy import units as u

from ..sampled import TimeSeries

INPUT_TIME = Time(['2016-03-22T12:30:31',
                   '2015-01-21T12:30:32',
                   '2016-03-22T12:30:40'])
PLAIN_TABLE = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a', 'b', 'c'])


def test_empty_initialization():
    ts = TimeSeries()
    ts['time'] = Time([1, 2, 3], format='mjd')


def test_empty_initialization_invalid():
    # Make sure things crash when the first column added is not a time column
    ts = TimeSeries()
    with pytest.raises(ValueError) as exc:
        ts['flux'] = [1, 2, 3]
    assert exc.value.args[0] == ("TimeSeries requires a column called "
                                 "'time' to be set before data can be added")


def test_initialize_only_time():
    ts = TimeSeries(time=INPUT_TIME)
    assert ts['time'] is ts.time
    # NOTE: the object in the table is a copy
    assert_equal(ts.time.isot, INPUT_TIME.isot)


def test_initialization_with_data():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert_equal(ts.time.isot, INPUT_TIME.isot)
    assert_equal(ts['a'], [10, 2, 3])
    assert_equal(ts['b'], [4, 5, 6])


def test_initialize_only_data():
    with pytest.raises(TypeError) as exc:
        TimeSeries(data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert exc.value.args[0] == "'time' has not been specified"


def test_initialization_with_table():
    ts = TimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    assert ts.colnames == ['time', 'a', 'b', 'c']


def test_initialization_with_time_delta():
    ts = TimeSeries(time=datetime(2018, 7, 1, 10, 10, 10),
                    time_delta=TimeDelta(3, format='sec'),
                    data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert_equal(ts.time.isot, ['2018-07-01T10:10:10.000',
                                '2018-07-01T10:10:13.000',
                                '2018-07-01T10:10:16.000'])


def test_initialization_missing_time_delta():
    with pytest.raises(TypeError) as exc:
        TimeSeries(time=datetime(2018, 7, 1, 10, 10, 10),
                   data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert exc.value.args[0] == "'time' is scalar, so 'time_delta' is required"


def test_initialization_invalid_time_delta():
    with pytest.raises(TypeError) as exc:
        TimeSeries(time=datetime(2018, 7, 1, 10, 10, 10),
                   time_delta=[1, 4, 3],
                   data=[[10, 2, 3], [4, 5, 6]], names=['a', 'b'])
    assert exc.value.args[0] == "'time_delta' should be a Quantity or a TimeDelta"


def test_initialization_with_time_in_data():

    data = PLAIN_TABLE.copy()
    data['time'] = INPUT_TIME

    ts1 = TimeSeries(data=data)

    assert set(ts1.colnames) == set(['time', 'a', 'b', 'c'])
    assert all(ts1.time == INPUT_TIME)

    ts2 = TimeSeries(data=[[10, 2, 3], INPUT_TIME], names=['a', 'time'])
    assert set(ts2.colnames) == set(['time', 'a'])
    assert all(ts2.time == INPUT_TIME)

    with pytest.raises(TypeError) as exc:
        # Don't allow ambiguous cases of passing multiple 'time' columns
        TimeSeries(data=data, time=INPUT_TIME)
    assert exc.value.args[0] == "'time' has been given both in the table and as a keyword argument"

    with pytest.raises(TypeError) as exc:
        # 'time' is a protected name, don't allow ambiguous cases
        TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], INPUT_TIME], names=['a', 'time'])
    assert exc.value.args[0] == "'time' has been given both in the table and as a keyword argument"


def test_initialization_n_samples():

    # Make sure things crash with incorrect n_samples

    with pytest.raises(TypeError) as exc:
        TimeSeries(time=INPUT_TIME, data=PLAIN_TABLE, n_samples=1000)
    assert exc.value.args[0] == ("'n_samples' has been given both and it is not the "
                                 "same length as the input data.")


def test_initialization_length_mismatch():
    with pytest.raises(ValueError) as exc:
        TimeSeries(time=INPUT_TIME, data=[[10, 2], [4, 5]], names=['a', 'b'])
    assert exc.value.args[0] == "Length of 'time' (3) should match data length (2)"


def test_initialization_invalid_both_time_and_time_delta():
    with pytest.raises(TypeError) as exc:
        TimeSeries(time=INPUT_TIME, time_delta=TimeDelta(3, format='sec'))
    assert exc.value.args[0] == ("'time_delta' should not be specified since "
                                 "'time' is an array")


def test_fold():

    times = Time([1, 2, 3, 8, 9, 12], format='unix')

    ts = TimeSeries(time=times)
    ts['flux'] = [1, 4, 4, 3, 2, 3]

    # Try without midpoint epoch, as it should default to the first time
    tsf = ts.fold(period=3 * u.s)
    assert isinstance(tsf.time, TimeDelta)
    assert_allclose(tsf.time.sec, [0, 1, -1, 1, -1, -1], rtol=1e-6)

    # Try with midpoint epoch
    tsf = ts.fold(period=4 * u.s, midpoint_epoch=Time(2.5, format='unix'))
    assert isinstance(tsf.time, TimeDelta)
    assert_allclose(tsf.time.sec, [-1.5, -0.5, 0.5, 1.5, -1.5, 1.5], rtol=1e-6)


def test_pandas():
    pandas = pytest.importorskip("pandas")

    df1 = pandas.DataFrame()
    df1['a'] = [1, 2, 3]
    df1.set_index(pandas.DatetimeIndex(INPUT_TIME.datetime64), inplace=True)

    ts = TimeSeries.from_pandas(df1)
    assert_equal(ts.time.isot, INPUT_TIME.isot)
    assert ts.colnames == ['time', 'a']
    assert len(ts.indices) == 1
    assert (ts.indices['time'].columns[0] == INPUT_TIME).all()

    ts_tcb = TimeSeries.from_pandas(df1, time_scale='tcb')
    assert ts_tcb.time.scale == 'tcb'

    df2 = ts.to_pandas()
    assert (df2.index.values == pandas.Index(INPUT_TIME.datetime64).values).all()
    assert df2.columns == pandas.Index(['a'])
    assert (df1['a'] == df2['a']).all()

    with pytest.raises(TypeError) as exc:
        TimeSeries.from_pandas(None)
    assert exc.value.args[0] == 'Input should be a pandas DataFrame'

    df4 = pandas.DataFrame()
    df4['a'] = [1, 2, 3]

    with pytest.raises(TypeError) as exc:
        TimeSeries.from_pandas(df4)
    assert exc.value.args[0] == 'DataFrame does not have a DatetimeIndex'

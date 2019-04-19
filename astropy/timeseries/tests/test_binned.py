# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
from numpy.testing import assert_equal, assert_allclose

from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.utils.data import get_pkg_data_filename

from astropy.timeseries.periodograms import BoxLeastSquares, LombScargle
from astropy.timeseries.binned import BinnedTimeSeries
from astropy.tests.helper import assert_quantity_allclose


CSV_FILE = get_pkg_data_filename('data/binned.csv')


def test_empty_initialization():
    ts = BinnedTimeSeries()
    ts['time_bin_start'] = Time([1, 2, 3], format='mjd')


def test_empty_initialization_invalid():

    # Make sure things crash when the first column added is not a time column

    ts = BinnedTimeSeries()
    with pytest.raises(ValueError) as exc:
        ts['flux'] = [1, 2, 3]
    assert exc.value.args[0] == ("BinnedTimeSeries object is invalid - expected "
                                 "'time_bin_start' as the first column but found 'flux'")


def test_initialization_time_bin_invalid():

    # Make sure things crash when time_bin_* is passed incorrectly.

    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries(data=[[1, 4, 3]])
    assert exc.value.args[0] == ("'time_bin_start' has not been specified")

    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries(time_bin_start='2016-03-22T12:30:31', data=[[1, 4, 3]])
    assert exc.value.args[0] == ("Either 'time_bin_size' or 'time_bin_end' should be specified")


def test_initialization_time_bin_both():

    # Make sure things crash when time_bin_* is passed twice.

    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries(data={"time_bin_start": ["2016-03-22T12:30:31"]},
                         time_bin_start="2016-03-22T12:30:31")
    assert exc.value.args[0] == ("'time_bin_start' has been given both in the table "
                                 "and as a keyword argument")

    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries(data={"time_bin_size": ["2016-03-22T12:30:31"]},
                         time_bin_size=[1]*u.s)
    assert exc.value.args[0] == ("'time_bin_size' has been given both in the table "
                                 "and as a keyword argument")


def test_initialization_time_bin_size():

    # Make sure things crash when time_bin_size has no units

    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries(data={"time": ["2016-03-22T12:30:31"]},
                         time_bin_start="2016-03-22T12:30:31",
                         time_bin_size=1)
    assert exc.value.args[0] == ("'time_bin_size' should be a Quantity or a TimeDelta")

    # TimeDelta for time_bin_size
    ts = BinnedTimeSeries(data={"time": ["2016-03-22T12:30:31"]},
                          time_bin_start="2016-03-22T12:30:31",
                          time_bin_size=TimeDelta(1))
    assert isinstance(ts.time_bin_size, u.quantity.Quantity)


def test_initialization_time_bin_start_scalar():

    # Make sure things crash when time_bin_start is a scalar with no time_bin_size

    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries(data={"time": ["2016-03-22T12:30:31"]},
                         time_bin_start=Time(1, format='mjd'),
                         time_bin_end=Time(1, format='mjd'))
    assert exc.value.args[0] == ("'time_bin_start' is scalar, so 'time_bin_size' is required")


def test_initialization_n_bins():

    # Make sure things crash with incorrect n_bins

    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries(data={"time": ["2016-03-22T12:30:31"]},
                         time_bin_start=Time(1, format='mjd'),
                         time_bin_size=1*u.s,
                         time_bin_end=Time(1, format='mjd'),
                         n_bins=10)
    assert exc.value.args[0] == ("'n_bins' has been given and it is not the "
                                 "same length as the input data.")


def test_initialization_non_scalar_time():

    # Make sure things crash with incorrect size of time_bin_start

    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries(data={"time": ["2016-03-22T12:30:31"]},
                         time_bin_start=["2016-03-22T12:30:31", "2016-03-22T12:30:32"],
                         time_bin_size=1*u.s,
                         time_bin_end=Time(1, format='mjd'))
    assert exc.value.args[0] == ("Length of 'time_bin_start' (2) should match table length (1)")

    with pytest.raises(TypeError) as exc:
        BinnedTimeSeries(data={"time": ["2016-03-22T12:30:31"]},
                         time_bin_start=["2016-03-22T12:30:31"],
                         time_bin_size=None,
                         time_bin_end=None)
    assert exc.value.args[0] == ("Either 'time_bin_size' or 'time_bin_end' should be specified")


def test_even_contiguous():

    # Initialize a ``BinnedTimeSeries`` with even contiguous bins by specifying
    # the bin width:

    ts = BinnedTimeSeries(time_bin_start='2016-03-22T12:30:31',
                          time_bin_size=3 * u.s, data=[[1, 4, 3]])

    assert_equal(ts.time_bin_start.isot, ['2016-03-22T12:30:31.000',
                                          '2016-03-22T12:30:34.000',
                                          '2016-03-22T12:30:37.000'])

    assert_equal(ts.time_bin_center.isot, ['2016-03-22T12:30:32.500',
                                           '2016-03-22T12:30:35.500',
                                           '2016-03-22T12:30:38.500'])

    assert_equal(ts.time_bin_end.isot, ['2016-03-22T12:30:34.000',
                                        '2016-03-22T12:30:37.000',
                                        '2016-03-22T12:30:40.000'])


def test_uneven_contiguous():

    # Initialize a ``BinnedTimeSeries`` with uneven contiguous bins by giving an
    # end time:

    ts = BinnedTimeSeries(time_bin_start=['2016-03-22T12:30:31',
                                          '2016-03-22T12:30:32',
                                          '2016-03-22T12:30:40'],
                          time_bin_end='2016-03-22T12:30:55',
                          data=[[1, 4, 3]])

    assert_equal(ts.time_bin_start.isot, ['2016-03-22T12:30:31.000',
                                          '2016-03-22T12:30:32.000',
                                          '2016-03-22T12:30:40.000'])

    assert_equal(ts.time_bin_center.isot, ['2016-03-22T12:30:31.500',
                                           '2016-03-22T12:30:36.000',
                                           '2016-03-22T12:30:47.500'])

    assert_equal(ts.time_bin_end.isot, ['2016-03-22T12:30:32.000',
                                        '2016-03-22T12:30:40.000',
                                        '2016-03-22T12:30:55.000'])


def test_uneven_non_contiguous():

    # Initialize a ``BinnedTimeSeries`` with uneven non-contiguous bins with
    # lists of start times, bin sizes and data:

    ts = BinnedTimeSeries(time_bin_start=['2016-03-22T12:30:31',
                                          '2016-03-22T12:30:38',
                                          '2016-03-22T12:34:40'],
                          time_bin_size=[5, 100, 2]*u.s,
                          data=[[1, 4, 3]])

    assert_equal(ts.time_bin_start.isot, ['2016-03-22T12:30:31.000',
                                          '2016-03-22T12:30:38.000',
                                          '2016-03-22T12:34:40.000'])

    assert_equal(ts.time_bin_center.isot, ['2016-03-22T12:30:33.500',
                                           '2016-03-22T12:31:28.000',
                                           '2016-03-22T12:34:41.000'])

    assert_equal(ts.time_bin_end.isot, ['2016-03-22T12:30:36.000',
                                        '2016-03-22T12:32:18.000',
                                        '2016-03-22T12:34:42.000'])


def test_uneven_non_contiguous_full():

    # Initialize a ``BinnedTimeSeries`` with uneven non-contiguous bins by
    # specifying the start and end times for the bins:

    ts = BinnedTimeSeries(time_bin_start=['2016-03-22T12:30:31',
                                          '2016-03-22T12:30:33',
                                          '2016-03-22T12:30:40'],
                          time_bin_end=['2016-03-22T12:30:32',
                                        '2016-03-22T12:30:35',
                                        '2016-03-22T12:30:41'],
                          data=[[1, 4, 3]])

    assert_equal(ts.time_bin_start.isot, ['2016-03-22T12:30:31.000',
                                          '2016-03-22T12:30:33.000',
                                          '2016-03-22T12:30:40.000'])

    assert_equal(ts.time_bin_center.isot, ['2016-03-22T12:30:31.500',
                                           '2016-03-22T12:30:34.000',
                                           '2016-03-22T12:30:40.500'])

    assert_equal(ts.time_bin_end.isot, ['2016-03-22T12:30:32.000',
                                        '2016-03-22T12:30:35.000',
                                        '2016-03-22T12:30:41.000'])


def test_read_empty():
    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries.read(CSV_FILE, format='csv')
    assert exc.value.args[0] == '``time_bin_start_column`` should be provided since the default Table readers are being used.'


def test_read_no_size_end():
    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='time_start', format='csv')
    assert exc.value.args[0] == 'Either `time_bin_end_column` or `time_bin_size_column` should be provided.'


def test_read_both_extra_bins():
    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='time_start', time_bin_end_column='END', time_bin_size_column='bin_size', format='csv')
    assert exc.value.args[0] == "Cannot specify both `time_bin_end_column` and `time_bin_size_column`."


def test_read_size_no_unit():
    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='time_start', time_bin_size_column='bin_size', format='csv')
    assert exc.value.args[0] == "The bin size unit should be specified as an astropy Unit using ``time_bin_size_unit``."


def test_read_start_time_missing():
    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='abc', time_bin_size_column='bin_size', time_bin_size_unit=u.second, format='csv')
    assert exc.value.args[0] == "Bin start time column 'abc' not found in the input data."


def test_read_end_time_missing():
    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='time_start', time_bin_end_column="missing", format='csv')
    assert exc.value.args[0] == "Bin end time column 'missing' not found in the input data."


def test_read_size_missing():
    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='time_start', time_bin_size_column="missing", time_bin_size_unit=u.second, format='csv')
    assert exc.value.args[0] == "Bin size column 'missing' not found in the input data."


def test_read_time_unit_missing():
    with pytest.raises(ValueError) as exc:
        BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='time_start', time_bin_size_column="bin_size", format='csv')
    assert exc.value.args[0] == "The bin size unit should be specified as an astropy Unit using ``time_bin_size_unit``."


def test_read():

    timeseries = BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='time_start',
                                       time_bin_end_column='time_end', format='csv')
    assert timeseries.colnames == ['time_bin_start', 'time_bin_size', 'bin_size', 'A', 'B', 'C', 'D', 'E', 'F']
    assert len(timeseries) == 10
    assert timeseries['B'].sum() == 1151.54

    timeseries = BinnedTimeSeries.read(CSV_FILE, time_bin_start_column='time_start',
                                       time_bin_size_column='bin_size',
                                       time_bin_size_unit=u.second, format='csv')
    assert timeseries.colnames == ['time_bin_start', 'time_bin_size', 'time_end', 'A', 'B', 'C', 'D', 'E', 'F']
    assert len(timeseries) == 10
    assert timeseries['B'].sum() == 1151.54


@pytest.mark.parametrize('cls', [BoxLeastSquares, LombScargle])
def test_periodogram(cls):

    # Note that we don't need to check the actual results from the periodogram
    # classes here since these are tested extensively in
    # astropy.timeseries.periodograms.

    ts = BinnedTimeSeries(time_bin_start='2016-03-22T12:30:31',
                          time_bin_size=3 * u.s, data=[[1, 4, 3], [3, 4, 3]], names=['a', 'b'])

    p1 = cls.from_timeseries(ts, 'a')
    assert isinstance(p1, cls)
    assert_allclose(p1.t.jd, ts.time_bin_center.jd)
    assert_equal(p1.y, ts['a'])
    assert p1.dy is None

    p2 = cls.from_timeseries(ts, 'a', uncertainty='b')
    assert_quantity_allclose(p2.dy, ts['b'])

    p3 = cls.from_timeseries(ts, 'a', uncertainty=0.1)
    assert_allclose(p3.dy, 0.1)

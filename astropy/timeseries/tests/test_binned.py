import pytest
from numpy.testing import assert_equal

from astropy import units as u
from astropy.table import Table
from astropy.time import Time, TimeDelta

from ..binned import BinnedTimeSeries


def test_empty_initialization():
    ts = BinnedTimeSeries()
    ts['time_bin_start'] = Time([1, 2, 3], format='mjd')


def test_empty_initialization_invalid():

    # Make sure things crash when the first column added is not a time column

    ts = BinnedTimeSeries()
    with pytest.raises(ValueError) as exc:
        ts['flux'] = [1, 2, 3]
    assert exc.value.args[0] == ("BinnedTimeSeries requires a column called "
                                 "'time_bin_start' to be set before data can be added")


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

import pytest
from numpy.testing import assert_equal

from ... import units as u
from ...table import Table
from ...time import Time
from ..binned import BinnedTimeSeries


INPUT_TIME = Time(['2016-03-22T12:30:31',
                   '2015-01-21T12:30:32',
                   '2016-03-22T12:30:40'])
PLAIN_TABLE = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=['a', 'b', 'c'])


def test_empty_initialization():
    ts = BinnedTimeSeries()
    ts['start_time'] = Time([1, 2, 3], format='mjd')


def test_empty_initialization_invalid():

    # Make sure things crash when the first column added is not a time column

    ts = BinnedTimeSeries()
    with pytest.raises(ValueError) as exc:
        ts['flux'] = [1, 2, 3]
    assert exc.value.args[0] == ("BinnedTimeSeries requires a column called "
                                 "'start_time' to be set before data can be added")


def test_even_contiguous():

    # Initialize a ``BinnedTimeSeries`` with even contiguous bins by specifying
    # the bin width:

    ts = BinnedTimeSeries(start_time='2016-03-22T12:30:31',
                          bin_size=3 * u.s, data=[[1, 4, 3]])
    assert_equal(ts.start_time.isot, ['2016-03-22T12:30:31.000',
                                      '2016-03-22T12:30:34.000',
                                      '2016-03-22T12:30:37.000'])
    assert_equal(ts.end_time.isot, ['2016-03-22T12:30:34.000',
                                    '2016-03-22T12:30:37.000',
                                    '2016-03-22T12:30:40.000'])


def test_uneven_contiguous():

    # Initialize a ``BinnedTimeSeries`` with uneven contiguous bins by giving an
    # end time:

    ts = BinnedTimeSeries(start_time=['2016-03-22T12:30:31',
                                      '2016-03-22T12:30:32',
                                      '2016-03-22T12:30:40'],
                          end_time='2016-03-22T12:30:55',
                          data=[[1, 4, 3]])
    assert_equal(ts.start_time.isot, ['2016-03-22T12:30:31.000',
                                      '2016-03-22T12:30:32.000',
                                      '2016-03-22T12:30:40.000'])
    assert_equal(ts.end_time.isot, ['2016-03-22T12:30:32.000',
                                    '2016-03-22T12:30:40.000',
                                    '2016-03-22T12:30:55.000'])


def test_uneven_non_contiguous():

    # Initialize a ``BinnedTimeSeries`` with uneven non-contiguous bins with
    # lists of start times, bin sizes and data:

    ts = BinnedTimeSeries(start_time=['2016-03-22T12:30:31',
                                      '2016-03-22T12:30:38',
                                      '2016-03-22T12:34:40'],
                          bin_size=[5, 100, 2]*u.s,
                          data=[[1, 4, 3]])
    assert_equal(ts.start_time.isot, ['2016-03-22T12:30:31.000',
                                      '2016-03-22T12:30:38.000',
                                      '2016-03-22T12:34:40.000'])
    assert_equal(ts.end_time.isot, ['2016-03-22T12:30:36.000',
                                    '2016-03-22T12:32:18.000',
                                    '2016-03-22T12:34:42.000'])


def test_uneven_non_contiguous_full():

    # Initialize a ``BinnedTimeSeries`` with uneven non-contiguous bins by
    # specifying the start and end times for the bins:

    ts = BinnedTimeSeries(start_time=['2016-03-22T12:30:31',
                                      '2016-03-22T12:30:33',
                                      '2016-03-22T12:30:40'],
                          end_time=['2016-03-22T12:30:32',
                                    '2016-03-22T12:30:35',
                                    '2016-03-22T12:30:41'],
                          data=[[1, 4, 3]])

    assert_equal(ts.start_time.isot, ['2016-03-22T12:30:31.000',
                                      '2016-03-22T12:30:33.000',
                                      '2016-03-22T12:30:40.000'])

    assert_equal(ts.end_time.isot, ['2016-03-22T12:30:32.000',
                                    '2016-03-22T12:30:35.000',
                                    '2016-03-22T12:30:41.000'])

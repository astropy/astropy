# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from numpy.testing import assert_equal

from astropy import units as u
from astropy.table import Table, QTable, vstack, join
from astropy.time import Time

from astropy.timeseries.sampled import TimeSeries
from astropy.timeseries.binned import BinnedTimeSeries


INPUT_TIME = Time(['2016-03-22T12:30:31', '2015-01-21T12:30:32', '2016-03-22T12:30:40'])
PLAIN_TABLE = Table([[1., 2., 11.], [3, 4, 1], ['x', 'y', 'z']], names=['a', 'b', 'c'])


class CommonTimeSeriesTests:

    def test_stacking(self):
        ts = vstack([self.series, self.series])
        assert isinstance(ts, self.series.__class__)

    def test_row_slicing(self):
        ts = self.series[:2]
        assert isinstance(ts, self.series.__class__)

    def test_row_indexing(self):
        self.series[0][self.time_attr] == Time('2015-01-21T12:30:32')
        self.series[self.time_attr][0] == Time('2015-01-21T12:30:32')

    def test_column_indexing(self):
        assert_equal(self.series['a'], [1, 2, 11])

    def test_column_slicing_notime(self):
        tab = self.series['a', 'b']
        assert not isinstance(tab, self.series.__class__)
        assert isinstance(tab, QTable)

    def test_add_column(self):
        self.series['d'] = [1, 2, 3]

    def test_add_row(self):
        self.series.add_row(self._row)

    def test_set_unit(self):
        self.series['d'] = [1, 2, 3]
        self.series['d'].unit = 's'

    def test_replace_column(self):
        self.series.replace_column('c', [1, 3, 4])

    def test_required_after_stacking(self):
        # When stacking, we have to temporarily relax the checking of the
        # columns in the time series, but we need to make sure that the
        # checking works again afterwards
        ts = vstack([self.series, self.series])
        with pytest.raises(ValueError) as exc:
            ts.remove_columns(ts.colnames)
        assert 'TimeSeries object is invalid' in exc.value.args[0]

    def test_join(self):
        ts_other = self.series.copy()
        ts_other.add_row(self._row)
        ts_other['d'] = [11, 22, 33, 44]
        ts_other.remove_columns(['a', 'b'])
        ts = join(self.series, ts_other)
        assert len(ts) == len(self.series)
        ts = join(self.series, ts_other, join_type='outer')
        assert len(ts) == len(ts_other)


class TestTimeSeries(CommonTimeSeriesTests):

    _row = {'time': '2016-03-23T12:30:40', 'a': 1., 'b': 2, 'c': 'a'}

    def setup_method(self, method):
        self.series = TimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
        self.time_attr = 'time'

    def test_column_slicing(self):
        ts = self.series['time', 'a']
        assert isinstance(ts, TimeSeries)


class TestBinnedTimeSeries(CommonTimeSeriesTests):

    _row = {'time_bin_start': '2016-03-23T12:30:40',
            'time_bin_size': 2 * u.s, 'a': 1., 'b': 2, 'c': 'a'}

    def setup_method(self, method):
        self.series = BinnedTimeSeries(time_bin_start=INPUT_TIME,
                                       time_bin_size=3 * u.s,
                                       data=PLAIN_TABLE)
        self.time_attr = 'time_bin_start'

    def test_column_slicing(self):
        ts = self.series['time_bin_start', 'time_bin_size', 'a']
        assert isinstance(ts, BinnedTimeSeries)

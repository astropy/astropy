# Licensed under a 3-clause BSD style license - see LICENSE.rst
import urllib.request

import pytest
import numpy as np

from astropy.time import Time
from astropy.utils.iers import iers
from astropy.utils.data import get_pkg_data_filename


# Test leap_seconds.list in test/data.
LEAP_SECOND_LIST = get_pkg_data_filename('data/leap-seconds.list')


class TestReading:
    """Basic tests that leap seconds can be read."""

    def verify_day_month_year(self, ls):
        assert np.all(ls['day'] == 1)
        assert np.all((ls['month'] == 1) | (ls['month'] == 7) |
                      (ls['year'] < 1970))
        assert np.all(ls['year'] >= 1960)
        t = Time({'year': ls['year'], 'month': ls['month'], 'day': ls['day']},
                 format='ymdhms')
        assert np.all(t == Time(ls['mjd'], format='mjd'))

    def test_read_leap_second_dat(self):
        ls = iers.LeapSeconds.from_iers_leap_seconds(
            iers.IERS_LEAP_SECOND_FILE)
        # Below, >= to take into account we might ship and updated file.
        assert ls.expires >= Time('2020-06-28')
        assert ls['mjd'][0] == 41317
        assert ls['tai_utc'][0] == 10
        assert ls['mjd'][-1] >= 57754
        assert ls['tai_utc'][-1] >= 37
        self.verify_day_month_year(ls)

    @pytest.mark.parametrize('file', (
        LEAP_SECOND_LIST,
        "file:" + urllib.request.pathname2url(LEAP_SECOND_LIST)))
    def test_read_leap_seconds_list(self, file):
        ls = iers.LeapSeconds.from_leap_seconds_list(file)
        assert ls.expires == Time('2020-06-28')
        assert ls['mjd'][0] == 41317
        assert ls['tai_utc'][0] == 10
        assert ls['mjd'][-1] == 57754
        assert ls['tai_utc'][-1] == 37
        self.verify_day_month_year(ls)

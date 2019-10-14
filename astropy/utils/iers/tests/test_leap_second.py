# Licensed under a 3-clause BSD style license - see LICENSE.rst
import urllib.request

import pytest
import numpy as np
from numpy.testing import assert_array_equal

from astropy import _erfa as erfa
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


class ERFALeapSecondsSafe:
    def setup(self):
        # Keep current leap-second table and expiration.
        self.erfa_ls = self._erfa_ls = erfa.leap_seconds.get()
        self._expires = erfa.leap_seconds._expires

    def teardown(self):
        # Restore leap-second table and expiration.
        erfa.leap_seconds.set(self.erfa_ls)
        erfa.leap_seconds._expires = self._expires


class TestFromERFA(ERFALeapSecondsSafe):
    def test_get_erfa_ls(self):
        ls = iers.LeapSeconds.from_erfa()
        assert ls.colnames == ['year', 'month', 'tai_utc']
        assert ls.expires == erfa.leap_seconds.expires
        ls_array = np.array(ls['year', 'month', 'tai_utc'])
        assert np.all(ls_array == self.erfa_ls)

    def test_get_modified_erfa_ls(self):
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        ls = iers.LeapSeconds.from_erfa()
        ls_array = np.array(ls['year', 'month', 'tai_utc'])
        assert np.all(ls_array == self.erfa_ls[:-2])


class TestUpdateLeapSeconds(ERFALeapSecondsSafe):
    def setup(self):
        super().setup()
        # Read default leap second table.
        self.ls = iers.LeapSeconds.from_iers_leap_seconds()
        # For tests, reset ERFA table to built-in default.
        erfa.leap_seconds.set()
        self.erfa_ls = erfa.leap_seconds.get()

    def test_built_in_up_to_date(self):
        """Leap second should match between built-in and ERFA."""
        erfa_since_1970 = self.erfa_ls[self.erfa_ls['year'] > 1970]
        assert len(self.ls) >= len(erfa_since_1970), \
            "built-in leap seconds out of date"
        assert len(self.ls) <= len(erfa_since_1970), \
            "ERFA leap seconds out of date"
        overlap = np.array(self.ls['year', 'month', 'tai_utc'])
        assert np.all(overlap == erfa_since_1970.astype(overlap.dtype))

    def test_update_with_built_in(self):
        """An update with built-in should not do anything."""
        n_update = self.ls.update_erfa_leap_seconds()
        assert n_update == 0
        new_erfa_ls = erfa.leap_seconds.get()
        assert np.all(new_erfa_ls == self.erfa_ls)

    @pytest.mark.parametrize('n_short', (1, 3))
    def test_update(self, n_short):
        """Check whether we can recover removed leap seconds."""
        erfa.leap_seconds.set(self.erfa_ls[:-n_short])
        n_update = self.ls.update_erfa_leap_seconds()
        assert n_update == n_short
        new_erfa_ls = erfa.leap_seconds.get()
        assert_array_equal(new_erfa_ls, self.erfa_ls)
        # Check that a second update does not do anything.
        n_update2 = self.ls.update_erfa_leap_seconds()
        assert n_update2 == 0
        new_erfa_ls2 = erfa.leap_seconds.get()
        assert_array_equal(new_erfa_ls2, self.erfa_ls)

    def test_update_initialize_erfa(self):
        # With pre-initialization, update does nothing.
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        n_update = self.ls.update_erfa_leap_seconds(initialize_erfa=True)
        assert n_update == 0
        new_erfa_ls = erfa.leap_seconds.get()
        assert_array_equal(new_erfa_ls, self.erfa_ls)

    def test_bad_jump(self):
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        bad = self.ls.copy()
        bad['tai_utc'][-1] = 5
        with pytest.raises(ValueError, match='jump'):
            bad.update_erfa_leap_seconds()
        # With an error the ERFA table should not change.
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls[:-2])

        # Unless we initialized it beforehand.
        with pytest.raises(ValueError, match='jump'):
            bad.update_erfa_leap_seconds(initialize_erfa=True)
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls)

        # Of course, we get no errors if we initialize only.
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        n_update = bad.update_erfa_leap_seconds(initialize_erfa='only')
        assert n_update == 0
        new_erfa_ls = erfa.leap_seconds.get()
        assert_array_equal(new_erfa_ls, self.erfa_ls)

    def test_bad_day(self):
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        bad = self.ls.copy()
        bad['day'][-1] = 5
        with pytest.raises(ValueError, match='not on 1st'):
            bad.update_erfa_leap_seconds()

    def test_bad_month(self):
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        bad = self.ls.copy()
        bad['month'][-1] = 5
        with pytest.raises(ValueError, match='January'):
            bad.update_erfa_leap_seconds()
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls[:-2])

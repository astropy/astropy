# Licensed under a 3-clause BSD style license - see LICENSE.rst
import urllib.request
import os
import locale
import platform

import pytest
import numpy as np
from numpy.testing import assert_array_equal
import erfa

from astropy.time import Time, TimeDelta
from astropy.utils.iers import iers
from astropy.utils.data import get_pkg_data_filename
from astropy.tests.tests.test_imports import test_imports

# Import every top-level astropy module as a test that the ERFA leap second
# table is not updated for normal imports.
test_imports()

# Now test that the erfa leap_seconds table has not been updated. This must be
# done at the module level, which unfortunately will abort the entire test run
# if if fails. Running within a normal pytest test will not work because the
# other tests will end up updating this attribute by virtue of doing Time UTC
# transformations.
assert erfa.leap_seconds._expires is None

# Tests in this module assume that the erfa.leap_seconds attribute has been
# updated from the `erfa` package built-in table to the astropy built-in
# leap-second table. That has the effect of ensuring that the
# `erfa.leap_seconds.expires` property is sufficiently in the future.
iers_table = iers.LeapSeconds.auto_open()
erfa.leap_seconds.update(iers_table)
assert erfa.leap_seconds._expires is not None

SYSTEM_FILE = '/usr/share/zoneinfo/leap-seconds.list'

# Test leap_seconds.list in test/data.
LEAP_SECOND_LIST = get_pkg_data_filename('data/leap-seconds.list')


def test_configuration():
    # This test just ensures things stay consistent.
    # Adjust if changes are made.
    assert iers.conf.iers_leap_second_auto_url == iers.IERS_LEAP_SECOND_URL
    assert iers.conf.ietf_leap_second_auto_url == iers.IETF_LEAP_SECOND_URL


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
        assert ls.expires >= Time('2020-06-28', scale='tai')
        assert ls['mjd'][0] == 41317
        assert ls['tai_utc'][0] == 10
        assert ls['mjd'][-1] >= 57754
        assert ls['tai_utc'][-1] >= 37
        self.verify_day_month_year(ls)

    def test_read_leap_second_dat_locale(self):
        current = locale.setlocale(locale.LC_ALL)
        try:
            if platform.system() == 'Darwin':
                locale.setlocale(locale.LC_ALL, 'fr_FR')
            else:
                locale.setlocale(locale.LC_ALL, 'fr_FR.utf8')

            ls = iers.LeapSeconds.from_iers_leap_seconds(
                iers.IERS_LEAP_SECOND_FILE)
        except locale.Error as e:
            pytest.skip(f'Locale error: {e}')
        finally:
            locale.setlocale(locale.LC_ALL, current)

        # Below, >= to take into account we might ship and updated file.
        assert ls.expires >= Time('2020-06-28', scale='tai')

    def test_open_leap_second_dat(self):
        ls = iers.LeapSeconds.from_iers_leap_seconds(
            iers.IERS_LEAP_SECOND_FILE)
        ls2 = iers.LeapSeconds.open(iers.IERS_LEAP_SECOND_FILE)
        assert np.all(ls == ls2)

    @pytest.mark.parametrize('file', (
        LEAP_SECOND_LIST,
        "file:" + urllib.request.pathname2url(LEAP_SECOND_LIST)))
    def test_read_leap_seconds_list(self, file):
        ls = iers.LeapSeconds.from_leap_seconds_list(file)
        assert ls.expires == Time('2020-06-28', scale='tai')
        assert ls['mjd'][0] == 41317
        assert ls['tai_utc'][0] == 10
        assert ls['mjd'][-1] == 57754
        assert ls['tai_utc'][-1] == 37
        self.verify_day_month_year(ls)

    @pytest.mark.parametrize('file', (
        LEAP_SECOND_LIST,
        "file:" + urllib.request.pathname2url(LEAP_SECOND_LIST)))
    def test_open_leap_seconds_list(self, file):
        ls = iers.LeapSeconds.from_leap_seconds_list(file)
        ls2 = iers.LeapSeconds.open(file)
        assert np.all(ls == ls2)

    @pytest.mark.skipif(not os.path.isfile(SYSTEM_FILE),
                        reason=f'system does not have {SYSTEM_FILE}')
    def test_open_system_file(self):
        ls = iers.LeapSeconds.open(SYSTEM_FILE)
        expired = ls.expires < Time.now()
        if expired:
            pytest.skip("System leap second file is expired.")
        assert not expired


def make_fake_file(expiration, tmpdir):
    """copy the built-in IERS file but set a different expiration date."""
    ls = iers.LeapSeconds.from_iers_leap_seconds()
    fake_file = str(tmpdir.join('fake_leap_seconds.dat'))
    with open(fake_file, 'w') as fh:
        fh.write('\n'.join([f'#  File expires on {expiration}']
                           + str(ls).split('\n')[2:-1]))
        return fake_file


def test_fake_file(tmpdir):
    fake_file = make_fake_file('28 June 2345', tmpdir)
    fake = iers.LeapSeconds.from_iers_leap_seconds(fake_file)
    assert fake.expires == Time('2345-06-28', scale='tai')


class TestAutoOpenExplicitLists:
    # For this set of tests, leap-seconds are allowed to be expired
    # except as explicitly tested.
    @pytest.mark.filterwarnings(iers.IERSStaleWarning)
    def test_auto_open_simple(self):
        ls = iers.LeapSeconds.auto_open([iers.IERS_LEAP_SECOND_FILE])
        assert ls.meta['data_url'] == iers.IERS_LEAP_SECOND_FILE

    @pytest.mark.filterwarnings(iers.IERSStaleWarning)
    def test_auto_open_erfa(self):
        ls = iers.LeapSeconds.auto_open(['erfa', iers.IERS_LEAP_SECOND_FILE])
        assert ls.meta['data_url'] in ['erfa', iers.IERS_LEAP_SECOND_FILE]

    @pytest.mark.filterwarnings(iers.IERSStaleWarning)
    def test_fake_future_file(self, tmpdir):
        fake_file = make_fake_file('28 June 2345', tmpdir)
        # Try as system file for auto_open, setting auto_max_age such
        # that any ERFA or system files are guaranteed to be expired,
        # while the fake file is guaranteed to be OK.
        with iers.conf.set_temp('auto_max_age', -100000):
            ls = iers.LeapSeconds.auto_open([
                'erfa', iers.IERS_LEAP_SECOND_FILE, fake_file])
            assert ls.expires == Time('2345-06-28', scale='tai')
            assert ls.meta['data_url'] == str(fake_file)
            # And as URL
            fake_url = "file:" + urllib.request.pathname2url(fake_file)
            ls2 = iers.LeapSeconds.auto_open([
                'erfa', iers.IERS_LEAP_SECOND_FILE, fake_url])
            assert ls2.expires == Time('2345-06-28', scale='tai')
            assert ls2.meta['data_url'] == str(fake_url)

    def test_fake_expired_file(self, tmpdir):
        fake_file1 = make_fake_file('28 June 2010', tmpdir)
        fake_file2 = make_fake_file('27 June 2012', tmpdir)
        # Between these and the built-in one, the built-in file is best.
        ls = iers.LeapSeconds.auto_open([fake_file1, fake_file2,
                                         iers.IERS_LEAP_SECOND_FILE])
        assert ls.meta['data_url'] == iers.IERS_LEAP_SECOND_FILE

        # But if we remove the built-in one, the least expired one will be
        # used and we get a warning that it is stale.
        with pytest.warns(iers.IERSStaleWarning):
            ls2 = iers.LeapSeconds.auto_open([fake_file1, fake_file2])
        assert ls2.meta['data_url'] == fake_file2
        assert ls2.expires == Time('2012-06-27', scale='tai')

        # Use the fake files to make sure auto_max_age is safe.
        # Should have no warning in either example.
        with iers.conf.set_temp('auto_max_age', None):
            ls3 = iers.LeapSeconds.auto_open([fake_file1,
                                              iers.IERS_LEAP_SECOND_FILE])
        assert ls3.meta['data_url'] == iers.IERS_LEAP_SECOND_FILE
        with iers.conf.set_temp('auto_max_age', None):
            ls4 = iers.LeapSeconds.auto_open([fake_file1, fake_file2])
        assert ls4.meta['data_url'] == fake_file2


@pytest.mark.remote_data
class TestRemoteURLs:
    def setup_class(cls):
        # Need auto_download so that IERS_B won't be loaded and cause tests to
        # fail.
        iers.conf.auto_download = True

    def teardown_class(cls):
        # This setting is to be consistent with astropy/conftest.py
        iers.conf.auto_download = False

    # In these tests, the results may be cached.
    # This is fine - no need to download again.
    def test_iers_url(self):
        ls = iers.LeapSeconds.auto_open([iers.IERS_LEAP_SECOND_URL])
        assert ls.expires > Time.now()

    def test_ietf_url(self):
        ls = iers.LeapSeconds.auto_open([iers.IETF_LEAP_SECOND_URL])
        assert ls.expires > Time.now()


class TestDefaultAutoOpen:
    """Test auto_open with different _auto_open_files."""
    def setup(self):
        # Identical to what is used in LeapSeconds.auto_open().
        self.good_enough = (iers.LeapSeconds._today()
                            + TimeDelta(180 - iers._none_to_float(iers.conf.auto_max_age),
                                        format='jd'))
        self._auto_open_files = iers.LeapSeconds._auto_open_files.copy()

    def teardown(self):
        iers.LeapSeconds._auto_open_files = self._auto_open_files

    def remove_auto_open_files(self, *files):
        """Remove some files from the auto-opener.

        The default set is restored in teardown.
        """
        for f in files:
            iers.LeapSeconds._auto_open_files.remove(f)

    def test_erfa_found(self):
        # Set huge maximum age such that whatever ERFA has is OK.
        # Since it is checked first, it should thus be found.
        with iers.conf.set_temp('auto_max_age', 100000):
            ls = iers.LeapSeconds.open()
        assert ls.meta['data_url'] == 'erfa'

    def test_builtin_found(self):
        # Set huge maximum age such that built-in file is always OK.
        # If we remove 'erfa', it should thus be found.
        self.remove_auto_open_files('erfa')
        with iers.conf.set_temp('auto_max_age', 100000):
            ls = iers.LeapSeconds.open()
        assert ls.meta['data_url'] == iers.IERS_LEAP_SECOND_FILE

    # The test below is marked remote_data only to ensure it runs
    # as an allowed-fail job on CI: i.e., we will notice it (eventually)
    # but will not be misled in thinking that a PR is bad.
    @pytest.mark.remote_data
    def test_builtin_not_expired(self):
        # TODO: would be nice to have automatic PRs for this!
        ls = iers.LeapSeconds.open(iers.IERS_LEAP_SECOND_FILE)
        assert ls.expires > self.good_enough, (
            "The leap second file built in to astropy is expired. Fix with:\n"
            "cd astropy/utils/iers/data/; . update_builtin_iers.sh\n"
            "and commit as a PR (for details, see release procedure).")

    def test_fake_future_file(self, tmpdir):
        fake_file = make_fake_file('28 June 2345', tmpdir)
        # Try as system file for auto_open, setting auto_max_age such
        # that any ERFA or system files are guaranteed to be expired.
        with iers.conf.set_temp('auto_max_age', -100000), \
                iers.conf.set_temp('system_leap_second_file', fake_file):
            ls = iers.LeapSeconds.open()
        assert ls.expires == Time('2345-06-28', scale='tai')
        assert ls.meta['data_url'] == str(fake_file)
        # And as URL
        fake_url = "file:" + urllib.request.pathname2url(fake_file)
        with iers.conf.set_temp('auto_max_age', -100000), \
                iers.conf.set_temp('iers_leap_second_auto_url', fake_url):
            ls2 = iers.LeapSeconds.open()
        assert ls2.expires == Time('2345-06-28', scale='tai')
        assert ls2.meta['data_url'] == str(fake_url)

    def test_fake_expired_file(self, tmpdir):
        self.remove_auto_open_files('erfa', 'iers_leap_second_auto_url',
                                    'ietf_leap_second_auto_url')
        fake_file = make_fake_file('28 June 2010', tmpdir)
        with iers.conf.set_temp('system_leap_second_file', fake_file):
            # If we try this directly, the built-in file will be found.
            ls = iers.LeapSeconds.open()
            assert ls.meta['data_url'] == iers.IERS_LEAP_SECOND_FILE

            # But if we remove the built-in one, the expired one will be
            # used and we get a warning that it is stale.
            self.remove_auto_open_files(iers.IERS_LEAP_SECOND_FILE)
            with pytest.warns(iers.IERSStaleWarning):
                ls2 = iers.LeapSeconds.open()
            assert ls2.meta['data_url'] == fake_file
            assert ls2.expires == Time('2010-06-28', scale='tai')

    @pytest.mark.skipif(not os.path.isfile(SYSTEM_FILE),
                        reason=f'system does not have {SYSTEM_FILE}')
    def test_system_file_used_if_not_expired(self, tmpdir):
        # We skip the test if the system file is on a CI and is expired -
        # we should not depend on CI keeping it up to date, but if it is,
        # we should check that it is used if possible.
        if (iers.LeapSeconds.open(SYSTEM_FILE).expires <= self.good_enough):
            pytest.skip("System leap second file is expired.")

        self.remove_auto_open_files('erfa')
        with iers.conf.set_temp('system_leap_second_file', SYSTEM_FILE):
            ls = iers.LeapSeconds.open()
            assert ls.expires > self.good_enough
            assert ls.meta['data_url'] in (iers.IERS_LEAP_SECOND_FILE,
                                           SYSTEM_FILE)

            # Also check with a "built-in" file that is expired
            fake_file = make_fake_file('28 June 2017', tmpdir)
            iers.LeapSeconds._auto_open_files[0] = fake_file
            ls2 = iers.LeapSeconds.open()
            assert ls2.expires > Time.now()
            assert ls2.meta['data_url'] == SYSTEM_FILE

    @pytest.mark.remote_data
    def test_auto_open_urls_always_good_enough(self):
        # Avoid using the erfa, built-in and system files, as they might
        # be good enough already.
        try:
            # Need auto_download so that IERS_B won't be loaded and
            # cause tests to fail.
            iers.conf.auto_download = True

            self.remove_auto_open_files('erfa', iers.IERS_LEAP_SECOND_FILE,
                                        'system_leap_second_file')
            ls = iers.LeapSeconds.open()
            assert ls.expires > self.good_enough
            assert ls.meta['data_url'].startswith('http')
        finally:
            # This setting is to be consistent with astropy/conftest.py
            iers.conf.auto_download = False


class ERFALeapSecondsSafe:
    """Base class for tests that change the ERFA leap-second tables.

    It ensures the original state is restored.
    """
    def setup(self):
        # Keep current leap-second table and expiration.
        self.erfa_ls = self._erfa_ls = erfa.leap_seconds.get()
        self.erfa_expires = self._expires = erfa.leap_seconds._expires

    def teardown(self):
        # Restore leap-second table and expiration.
        erfa.leap_seconds.set(self.erfa_ls)
        erfa.leap_seconds._expires = self._expires


class TestFromERFA(ERFALeapSecondsSafe):
    def test_get_erfa_ls(self):
        ls = iers.LeapSeconds.from_erfa()
        assert ls.colnames == ['year', 'month', 'tai_utc']
        assert isinstance(ls.expires, Time)
        assert ls.expires == self.erfa_expires
        ls_array = np.array(ls['year', 'month', 'tai_utc'])
        assert np.all(ls_array == self.erfa_ls)

    def test_get_built_in_erfa_ls(self):
        ls = iers.LeapSeconds.from_erfa(built_in=True)
        assert ls.colnames == ['year', 'month', 'tai_utc']
        assert isinstance(ls.expires, Time)
        ls_array = np.array(ls['year', 'month', 'tai_utc'])
        assert np.all(ls_array == self.erfa_ls[:len(ls_array)])

    def test_get_modified_erfa_ls(self):
        erfa.leap_seconds.set(self.erfa_ls[:-10])
        ls = iers.LeapSeconds.from_erfa()
        assert len(ls) == len(self.erfa_ls)-10
        ls_array = np.array(ls['year', 'month', 'tai_utc'])
        assert np.all(ls_array == self.erfa_ls[:-10])
        ls2 = iers.LeapSeconds.from_erfa(built_in=True)
        assert len(ls2) > len(ls)
        erfa.leap_seconds.set(None)
        erfa_built_in = erfa.leap_seconds.get()
        assert len(ls2) == len(erfa_built_in)
        ls2_array = np.array(ls2['year', 'month', 'tai_utc'])
        assert np.all(ls2_array == erfa_built_in)

    def test_open(self):
        ls = iers.LeapSeconds.open('erfa')
        assert isinstance(ls.expires, Time)
        assert ls.expires == self.erfa_expires
        ls_array = np.array(ls['year', 'month', 'tai_utc'])
        assert np.all(ls_array == self.erfa_ls)


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

    def test_update_overwrite(self):
        n_update = self.ls.update_erfa_leap_seconds(initialize_erfa='empty')
        assert n_update == len(self.ls)
        new_erfa_ls = erfa.leap_seconds.get()
        assert new_erfa_ls['year'].min() > 1970
        n_update2 = self.ls.update_erfa_leap_seconds()
        assert n_update2 == 0
        new_erfa_ls2 = erfa.leap_seconds.get()
        assert_array_equal(new_erfa_ls2, new_erfa_ls)
        n_update3 = self.ls.update_erfa_leap_seconds(initialize_erfa=True)
        assert n_update3 == 0
        new_erfa_ls3 = erfa.leap_seconds.get()
        assert_array_equal(new_erfa_ls3, self.erfa_ls)

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

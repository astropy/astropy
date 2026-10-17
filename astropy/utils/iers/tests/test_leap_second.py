# Licensed under a 3-clause BSD style license - see LICENSE.rst
import locale
import os
import platform
import urllib.request
from datetime import datetime, timedelta

import erfa
import numpy as np
import pytest
from astropy_iers_data import IERS_LEAP_SECOND_FILE
from numpy.testing import assert_array_equal

from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.iers import iers

SYSTEM_FILE = "/usr/share/zoneinfo/leap-seconds.list"

# Test leap_seconds.list in test/data.
LEAP_SECOND_LIST = get_pkg_data_filename("data/leap-seconds.list")


# Reading leap-second files is deprecated, but still tested until removal.
@pytest.mark.filterwarnings(
    "ignore::astropy.utils.exceptions.AstropyDeprecationWarning"
)
class TestReading:
    """Basic tests that leap seconds can be read."""

    def verify_day_month_year(self, ls):
        assert np.all(ls["day"] == 1)
        assert np.all((ls["month"] == 1) | (ls["month"] == 7) | (ls["year"] < 1970))
        assert np.all(ls["year"] >= 1960)
        t = Time(
            {"year": ls["year"], "month": ls["month"], "day": ls["day"]},
            format="ymdhms",
        )
        assert np.all(t == Time(ls["mjd"], format="mjd"))

    def test_read_leap_second_dat(self):
        ls = iers.LeapSeconds.from_iers_leap_seconds(IERS_LEAP_SECOND_FILE)
        # Below, >= to take into account we might ship and updated file.
        assert ls.expires >= Time("2020-06-28", scale="tai")
        assert ls["mjd"][0] == 41317
        assert ls["tai_utc"][0] == 10
        assert ls["mjd"][-1] >= 57754
        assert ls["tai_utc"][-1] >= 37
        self.verify_day_month_year(ls)

    def test_read_leap_second_dat_locale(self):
        current = locale.setlocale(locale.LC_ALL)
        try:
            if platform.system() == "Darwin":
                locale.setlocale(locale.LC_ALL, "fr_FR")
            else:
                locale.setlocale(locale.LC_ALL, "fr_FR.utf8")

            ls = iers.LeapSeconds.from_iers_leap_seconds(IERS_LEAP_SECOND_FILE)
        except locale.Error as e:
            pytest.skip(f"Locale error: {e}")
        finally:
            locale.setlocale(locale.LC_ALL, current)

        # Below, >= to take into account we might ship and updated file.
        assert ls.expires >= Time("2020-06-28", scale="tai")

    def test_open_leap_second_dat(self):
        ls = iers.LeapSeconds.from_iers_leap_seconds(IERS_LEAP_SECOND_FILE)
        ls2 = iers.LeapSeconds.open(IERS_LEAP_SECOND_FILE)
        assert np.all(ls == ls2)

    @pytest.mark.parametrize(
        "file",
        (LEAP_SECOND_LIST, "file:" + urllib.request.pathname2url(LEAP_SECOND_LIST)),
    )
    def test_read_leap_seconds_list(self, file):
        ls = iers.LeapSeconds.from_leap_seconds_list(file)
        assert ls.expires == Time("2020-06-28", scale="tai")
        assert ls["mjd"][0] == 41317
        assert ls["tai_utc"][0] == 10
        assert ls["mjd"][-1] == 57754
        assert ls["tai_utc"][-1] == 37
        self.verify_day_month_year(ls)

    @pytest.mark.parametrize(
        "file",
        (LEAP_SECOND_LIST, "file:" + urllib.request.pathname2url(LEAP_SECOND_LIST)),
    )
    def test_open_leap_seconds_list(self, file):
        ls = iers.LeapSeconds.from_leap_seconds_list(file)
        ls2 = iers.LeapSeconds.open(file)
        assert np.all(ls == ls2)

    @pytest.mark.skipif(
        not os.path.isfile(SYSTEM_FILE), reason=f"system does not have {SYSTEM_FILE}"
    )
    def test_open_system_file(self):
        ls = iers.LeapSeconds.open(SYSTEM_FILE)
        assert ls["mjd"][-1] == 57754
        assert ls["tai_utc"][-1] == 37


def make_fake_file(expiration, tmp_path):
    """copy the built-in IERS file but set a different expiration date."""
    ls = iers.LeapSeconds._from_iers_leap_seconds()
    fake_file = str(tmp_path / "fake_leap_seconds.dat")
    with open(fake_file, "w") as fh:
        fh.write(
            "\n".join([f"#  File expires on {expiration}"] + str(ls).split("\n")[2:-1])
        )
        return fake_file


@pytest.mark.filterwarnings(
    "ignore::astropy.utils.exceptions.AstropyDeprecationWarning"
)
def test_fake_file(tmp_path):
    fake_file = make_fake_file("28 June 2345", tmp_path)
    fake = iers.LeapSeconds.from_iers_leap_seconds(fake_file)
    assert fake.expires == Time("2345-06-28", scale="tai")


class TestDeprecatedAutoOpen:
    """Leap seconds were discontinued, so finding a newer table is deprecated."""

    def test_auto_open_default(self):
        with pytest.warns(AstropyDeprecationWarning, match="auto_open is deprecated"):
            ls = iers.LeapSeconds.auto_open()
        assert ls.meta["data_url"] == IERS_LEAP_SECOND_FILE
        assert np.all(ls == iers.LeapSeconds._from_iers_leap_seconds())

    def test_open_none(self):
        with pytest.warns(AstropyDeprecationWarning, match="other than 'erfa'"):
            ls = iers.LeapSeconds.open(None)
        assert ls.meta["data_url"] == IERS_LEAP_SECOND_FILE

    def test_auto_open_picks_latest_expiration(self, tmp_path):
        fake_file = make_fake_file("28 June 2345", tmp_path)
        with pytest.warns(AstropyDeprecationWarning, match="auto_open is deprecated"):
            ls = iers.LeapSeconds.auto_open(
                ["nonsense", IERS_LEAP_SECOND_FILE, fake_file]
            )
        assert ls.meta["data_url"] == fake_file
        assert ls.expires == Time("2345-06-28", scale="tai")

    def test_auto_open_no_readable_files(self):
        with (
            pytest.warns(AstropyDeprecationWarning, match="auto_open is deprecated"),
            pytest.raises(ValueError, match="none of the files could be read"),
        ):
            iers.LeapSeconds.auto_open(["nonsense"])


class ERFALeapSecondsSafe:
    """Base class for tests that change the ERFA leap-second tables.

    It ensures the original state is restored.
    """

    def setup_method(self):
        # Keep current leap-second table and expiration.
        self.erfa_ls = self._erfa_ls = erfa.leap_seconds.get()
        self._expires = erfa.leap_seconds._expires
        # Expected LeapSeconds.from_erfa().expires: leap seconds were
        # discontinued, so an expired ERFA table is treated as valid for
        # another 180 days.  Compare dates, since pyerfa may give a time too.
        if erfa.leap_seconds.expired:
            self.erfa_expires = datetime.now().date() + timedelta(180)
        else:
            self.erfa_expires = erfa.leap_seconds.expires.date()

    def teardown_method(self):
        # Restore leap-second table and expiration.
        erfa.leap_seconds.set(self.erfa_ls)
        erfa.leap_seconds._expires = self._expires


class TestFromERFA(ERFALeapSecondsSafe):
    def test_get_erfa_ls(self):
        ls = iers.LeapSeconds.from_erfa()
        assert ls.colnames == ["year", "month", "tai_utc", "day", "mjd"]
        assert isinstance(ls.expires, Time)
        assert ls.expires.datetime.date() == self.erfa_expires
        ls_array = np.array(ls["year", "month", "tai_utc"])
        assert np.all(ls_array == self.erfa_ls)

    def test_get_built_in_erfa_ls(self):
        ls = iers.LeapSeconds.from_erfa(built_in=True)
        assert ls.colnames == ["year", "month", "tai_utc", "day", "mjd"]
        assert isinstance(ls.expires, Time)
        ls_array = np.array(ls["year", "month", "tai_utc"])
        assert np.all(ls_array == self.erfa_ls[: len(ls_array)])

    def test_get_modified_erfa_ls(self):
        erfa.leap_seconds.set(self.erfa_ls[:-10])
        ls = iers.LeapSeconds.from_erfa()
        assert len(ls) == len(self.erfa_ls) - 10
        ls_array = np.array(ls["year", "month", "tai_utc"])
        assert np.all(ls_array == self.erfa_ls[:-10])
        ls2 = iers.LeapSeconds.from_erfa(built_in=True)
        assert len(ls2) > len(ls)
        erfa.leap_seconds.set(None)
        erfa_built_in = erfa.leap_seconds.get()
        assert len(ls2) == len(erfa_built_in)
        ls2_array = np.array(ls2["year", "month", "tai_utc"])
        assert np.all(ls2_array == erfa_built_in)

    def test_open(self):
        ls = iers.LeapSeconds.open("erfa")
        assert isinstance(ls.expires, Time)
        assert ls.expires.datetime.date() == self.erfa_expires
        ls_array = np.array(ls["year", "month", "tai_utc"])
        assert np.all(ls_array == self.erfa_ls)

    def test_day_mjd_columns(self):
        ls = iers.LeapSeconds.from_erfa()
        assert np.all(ls["day"] == 1)
        t = Time(
            {"year": ls["year"], "month": ls["month"], "day": ls["day"]},
            format="ymdhms",
            scale="tai",
        )
        assert np.all(ls["mjd"] == t.mjd)
        # From 1972 on, the table matches that in the bundled IERS file.
        ls_file = iers.LeapSeconds._from_iers_leap_seconds()
        ls_1972 = ls[ls["year"] >= 1972]
        for name in ("year", "month", "day", "mjd", "tai_utc"):
            assert np.all(ls_1972[name] == ls_file[name])

    def test_open_erfa_is_default(self):
        ls = iers.LeapSeconds.open()
        assert np.all(ls == iers.LeapSeconds.from_erfa())


class TestDeprecations(ERFALeapSecondsSafe):
    """Leap seconds were discontinued, so only the ERFA table remains."""

    @pytest.mark.parametrize(
        "method", ["from_iers_leap_seconds", "from_leap_seconds_list", "open"]
    )
    def test_reading_files(self, method):
        file = LEAP_SECOND_LIST if method == "from_leap_seconds_list" else None
        file = file or IERS_LEAP_SECOND_FILE
        with pytest.warns(AstropyDeprecationWarning, match="from_erfa"):
            ls = getattr(iers.LeapSeconds, method)(file)
        assert ls["tai_utc"][-1] == 37

    def test_update_erfa_leap_seconds(self):
        ls = iers.LeapSeconds.from_erfa()
        with pytest.warns(AstropyDeprecationWarning, match="erfa.leap_seconds.set"):
            assert ls.update_erfa_leap_seconds() == 0

    @pytest.mark.parametrize(
        "name",
        ["IERS_LEAP_SECOND_FILE", "IERS_LEAP_SECOND_URL", "IETF_LEAP_SECOND_URL"],
    )
    def test_constants(self, name):
        import astropy_iers_data

        from astropy.utils import iers as iers_package

        expected = getattr(
            astropy_iers_data,
            "IERS_LEAP_SECOND_URL_MIRROR" if name.startswith("IETF") else name,
        )
        for module in (iers, iers_package):
            with pytest.warns(AstropyDeprecationWarning, match=f"{name} is deprecated"):
                assert getattr(module, name) == expected
        assert name not in iers.__all__

    def test_missing_attribute(self):
        from astropy.utils import iers as iers_package

        for module in (iers, iers_package):
            with pytest.raises(AttributeError, match="no attribute 'nonsense'"):
                module.nonsense


# Updating ERFA from a LeapSeconds table is deprecated, but still tested.
@pytest.mark.filterwarnings(
    "ignore::astropy.utils.exceptions.AstropyDeprecationWarning"
)
class TestUpdateLeapSeconds(ERFALeapSecondsSafe):
    def setup_method(self):
        super().setup_method()
        # Read default leap second table.
        self.ls = iers.LeapSeconds.from_iers_leap_seconds()
        # For tests, reset ERFA table to built-in default.
        erfa.leap_seconds.set()
        self.erfa_ls = erfa.leap_seconds.get()

    def test_built_in_up_to_date(self):
        """Leap second should match between built-in and ERFA."""
        erfa_since_1970 = self.erfa_ls[self.erfa_ls["year"] > 1970]
        assert len(self.ls) >= len(erfa_since_1970), "built-in leap seconds out of date"
        assert len(self.ls) <= len(erfa_since_1970), "ERFA leap seconds out of date"
        overlap = np.array(self.ls["year", "month", "tai_utc"])
        assert np.all(overlap == erfa_since_1970.astype(overlap.dtype))

    def test_update_with_built_in(self):
        """An update with built-in should not do anything."""
        n_update = self.ls.update_erfa_leap_seconds()
        assert n_update == 0
        new_erfa_ls = erfa.leap_seconds.get()
        assert np.all(new_erfa_ls == self.erfa_ls)

    @pytest.mark.parametrize("n_short", (1, 3))
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
        n_update = self.ls.update_erfa_leap_seconds(initialize_erfa="empty")
        assert n_update == len(self.ls)
        new_erfa_ls = erfa.leap_seconds.get()
        assert new_erfa_ls["year"].min() > 1970
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
        bad["tai_utc"][-1] = 5
        with pytest.raises(ValueError, match="jump"):
            bad.update_erfa_leap_seconds()
        # With an error the ERFA table should not change.
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls[:-2])

        # Unless we initialized it beforehand.
        with pytest.raises(ValueError, match="jump"):
            bad.update_erfa_leap_seconds(initialize_erfa=True)
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls)

        # Of course, we get no errors if we initialize only.
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        n_update = bad.update_erfa_leap_seconds(initialize_erfa="only")
        assert n_update == 0
        new_erfa_ls = erfa.leap_seconds.get()
        assert_array_equal(new_erfa_ls, self.erfa_ls)

    def test_bad_day(self):
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        bad = self.ls.copy()
        bad["day"][-1] = 5
        with pytest.raises(ValueError, match="not on 1st"):
            bad.update_erfa_leap_seconds()

    def test_bad_month(self):
        erfa.leap_seconds.set(self.erfa_ls[:-2])
        bad = self.ls.copy()
        bad["month"][-1] = 5
        with pytest.raises(ValueError, match="January"):
            bad.update_erfa_leap_seconds()
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls[:-2])

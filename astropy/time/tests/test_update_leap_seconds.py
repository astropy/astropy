# Licensed under a 3-clause BSD style license - see LICENSE.rst
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta

import erfa
import pytest

import astropy.time.core
from astropy.time import Time, update_leap_seconds
from astropy.utils import iers
from astropy.utils.exceptions import AstropyWarning


class TestUpdateLeapSeconds:
    def setup_method(self):
        self.built_in = iers.LeapSeconds.from_iers_leap_seconds()
        self.erfa_ls = iers.LeapSeconds.from_erfa()
        now = datetime.now()
        self.good_enough = now + timedelta(150)

    def teardown_method(self):
        self.erfa_ls.update_erfa_leap_seconds(initialize_erfa=True)

    def test_auto_update_leap_seconds(self):
        # Sanity check.
        assert erfa.dat(2018, 1, 1, 0.0) == 37.0
        # Set expired leap seconds
        expired = self.erfa_ls[self.erfa_ls["year"] < 2017]
        expired.update_erfa_leap_seconds(initialize_erfa="empty")
        # Check the 2017 leap second is indeed missing.
        assert erfa.dat(2018, 1, 1, 0.0) == 36.0

        # Update with missing leap seconds.
        n_update = update_leap_seconds([iers.IERS_LEAP_SECOND_FILE])
        assert n_update >= 1
        assert erfa.leap_seconds.expires == self.built_in.expires
        assert erfa.dat(2018, 1, 1, 0.0) == 37.0

        # Doing it again does not change anything
        n_update2 = update_leap_seconds([iers.IERS_LEAP_SECOND_FILE])
        assert n_update2 == 0
        assert erfa.dat(2018, 1, 1, 0.0) == 37.0

    @pytest.mark.remote_data
    def test_never_expired_if_connected(self):
        assert self.erfa_ls.expires > datetime.now()
        assert self.erfa_ls.expires >= self.good_enough

    @pytest.mark.remote_data
    def test_auto_update_always_good(self):
        self.erfa_ls.update_erfa_leap_seconds(initialize_erfa="only")
        update_leap_seconds()
        assert not erfa.leap_seconds.expired
        assert erfa.leap_seconds.expires > self.good_enough

    def test_auto_update_bad_file(self):
        with pytest.warns(AstropyWarning, match="FileNotFound"):
            update_leap_seconds(["nonsense"])

    def test_auto_update_corrupt_file(self, tmp_path):
        bad_file = str(tmp_path / "no_expiration")
        with open(iers.IERS_LEAP_SECOND_FILE) as fh:
            lines = fh.readlines()
        with open(bad_file, "w") as fh:
            fh.write("\n".join([line for line in lines if not line.startswith("#")]))

        with pytest.warns(AstropyWarning, match="ValueError.*did not find expiration"):
            update_leap_seconds([bad_file])

    def test_auto_update_expired_file(self, tmp_path):
        # Set up expired ERFA leap seconds.
        expired = self.erfa_ls[self.erfa_ls["year"] < 2017]
        expired.update_erfa_leap_seconds(initialize_erfa="empty")
        # Create similarly expired file.
        expired_file = str(tmp_path / "expired.dat")
        with open(expired_file, "w") as fh:
            fh.write(
                "\n".join(
                    ["# File expires on 28 June 2010"] + [str(item) for item in expired]
                )
            )

        with pytest.warns(iers.IERSStaleWarning):
            update_leap_seconds(["erfa", expired_file])

    def test_init_thread_safety(self, monkeypatch):
        # Set up expired ERFA leap seconds.
        expired = self.erfa_ls[self.erfa_ls["year"] < 2017]
        expired.update_erfa_leap_seconds(initialize_erfa="empty")
        # Force re-initialization, even if another test already did it
        monkeypatch.setattr(
            astropy.time.core,
            "_LEAP_SECONDS_CHECK",
            astropy.time.core._LeapSecondsCheck.NOT_STARTED,
        )
        workers = 4
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(lambda: str(Time("2019-01-01 00:00:00.000").tai))
                for i in range(workers)
            ]
            results = [future.result() for future in futures]
            assert results == ["2019-01-01 00:00:37.000"] * workers

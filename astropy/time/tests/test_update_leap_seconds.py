# Licensed under a 3-clause BSD style license - see LICENSE.rst

import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path
from textwrap import dedent

import erfa
import pytest

from astropy.time import update_leap_seconds
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
        bad_file = tmp_path / "no_expiration"
        lines = Path(iers.IERS_LEAP_SECOND_FILE).read_text().splitlines()
        bad_file.write_text(
            "\n".join(line for line in lines if not line.startswith("#"))
        )

        with pytest.warns(AstropyWarning, match="ValueError.*did not find expiration"):
            update_leap_seconds([bad_file])

    def test_auto_update_expired_file(self, tmp_path):
        # Copy the built-in file, but with an expiration date in the past.
        # Leap seconds were discontinued, so the ERFA table itself is never
        # treated as expired; only a file can be.
        expired_file = tmp_path / "expired.dat"
        lines = Path(iers.IERS_LEAP_SECOND_FILE).read_text().splitlines()
        expired_file.write_text(
            "\n".join(
                "#  File expires on 28 June 2010" if "File expires" in line else line
                for line in lines
            )
        )
        with pytest.warns(iers.IERSStaleWarning):
            update_leap_seconds([expired_file])


@pytest.mark.skipif(
    sys.platform == "emscripten", reason="subprocesses not supported on emscripten"
)
def test_utc_conversion_does_not_update_erfa():
    """Leap seconds were discontinued, so Time no longer updates ERFA's table.

    Getting the table with LeapSeconds.from_erfa() ensures that the ERFA
    table is not considered expired, also for pyerfa versions that do not
    know leap seconds were discontinued.  This needs a fresh process, since
    the test process has already done UTC conversions and other tests change
    the ERFA table.
    """
    cmd = dedent("""
    from datetime import datetime, timedelta

    import erfa
    import numpy as np

    built_in = erfa.leap_seconds.get().copy()
    expires = erfa.leap_seconds._expires

    from astropy.time import Time
    from astropy.utils.iers import LeapSeconds

    assert Time("2020-01-01", scale="utc").tai.iso == "2020-01-01 00:00:37.000"
    assert Time("2020-01-01", scale="tai").utc.iso == "2019-12-31 23:59:23.000"
    assert erfa.leap_seconds._expires is expires
    assert np.all(erfa.leap_seconds.get() == built_in)

    ls = LeapSeconds.from_erfa()
    assert not erfa.leap_seconds.expired
    assert erfa.leap_seconds.expires > datetime.now() + timedelta(179)
    assert ls.expires.datetime.date() == erfa.leap_seconds.expires.date()
    assert np.all(erfa.leap_seconds.get() == built_in)
    """)
    subprocess.check_call([sys.executable, "-c", cmd])

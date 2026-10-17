# Licensed under a 3-clause BSD style license - see LICENSE.rst

import subprocess
import sys
from textwrap import dedent

import erfa
import pytest

from astropy.time import update_leap_seconds
from astropy.utils.exceptions import AstropyDeprecationWarning


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


class TestUpdateLeapSeconds:
    def setup_method(self):
        self.erfa_ls = erfa.leap_seconds.get()
        self.erfa_expires = erfa.leap_seconds._expires

    def teardown_method(self):
        erfa.leap_seconds.set(self.erfa_ls)
        erfa.leap_seconds._expires = self.erfa_expires

    def test_deprecated_but_still_updates(self):
        # Remove the 2017 leap second from ERFA.
        expired = self.erfa_ls[self.erfa_ls["year"] < 2017]
        erfa.leap_seconds.set(expired)
        assert erfa.dat(2018, 1, 1, 0.0) == 36.0

        with pytest.warns(AstropyDeprecationWarning, match="discontinued"):
            n_update = update_leap_seconds()
        assert n_update == 1
        assert erfa.dat(2018, 1, 1, 0.0) == 37.0

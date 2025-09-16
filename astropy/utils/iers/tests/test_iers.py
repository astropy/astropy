# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import re
import warnings
from pathlib import Path

import numpy as np
import pytest
from astropy_iers_data import (
    IERS_A_README,
    IERS_B_FILE,
    IERS_B_README,
    IERS_LEAP_SECOND_FILE,
)

from astropy import units as u
from astropy.config import set_temp_cache
from astropy.table import QTable
from astropy.tests.helper import CI, assert_quantity_allclose
from astropy.time import Time, TimeDelta
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.iers import iers

FILE_NOT_FOUND_ERROR = getattr(__builtins__, "FileNotFoundError", OSError)

IERS_A_EXCERPT = get_pkg_data_filename(os.path.join("data", "iers_a_excerpt"))


class TestBasic:
    """Basic tests that IERS_B returns correct values"""

    @pytest.mark.parametrize("iers_cls", (iers.IERS_B, iers.IERS))
    def test_simple(self, iers_cls):
        """Test the default behaviour for IERS_B and IERS."""
        # Arguably, IERS itself should not be used at all, but it used to
        # provide IERS_B by default so we check that it continues to do so.
        # Eventually, IERS should probably be deprecated.
        iers_cls.close()
        assert iers_cls.iers_table is None
        iers_tab = iers_cls.open()
        assert iers_cls.iers_table is not None
        assert iers_cls.iers_table is iers_tab
        assert isinstance(iers_tab, QTable)
        assert isinstance(iers_tab, iers.IERS_B)
        assert (iers_tab["UT1_UTC"].unit / u.second).is_unity()
        assert (iers_tab["PM_x"].unit / u.arcsecond).is_unity()
        assert (iers_tab["PM_y"].unit / u.arcsecond).is_unity()
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5, 2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0.0, 0.5])
        ut1_utc = iers_tab.ut1_utc(jd1, jd2)
        assert isinstance(ut1_utc, u.Quantity)
        assert (ut1_utc.unit / u.second).is_unity()
        # IERS files change at the 0.1 ms level; see gh-6981
        assert_quantity_allclose(
            ut1_utc,
            [-0.5868211, -0.5868184, -0.5868184, 0.4131816, 0.41328895] * u.s,
            atol=0.1 * u.ms,
        )
        # should be future-proof; surely we've moved to another planet by then
        with pytest.raises(IndexError):
            ut1_utc2, status2 = iers_tab.ut1_utc(1e11, 0.0)
        # also check it returns the right status
        ut1_utc2, status2 = iers_tab.ut1_utc(jd1, jd2, return_status=True)
        assert np.all(status2 == iers.FROM_IERS_B)
        ut1_utc4, status4 = iers_tab.ut1_utc(1e11, 0.0, return_status=True)
        assert status4 == iers.TIME_BEYOND_IERS_RANGE

        # check it works via Time too
        t = Time(jd1, jd2, format="jd", scale="utc")
        ut1_utc3 = iers_tab.ut1_utc(t)
        assert_quantity_allclose(
            ut1_utc3,
            [-0.5868211, -0.5868184, -0.5868184, 0.4131816, 0.41328895] * u.s,
            atol=0.1 * u.ms,
        )

        # Table behaves properly as a table (e.g. can be sliced)
        assert len(iers_tab[:2]) == 2

    def test_open_filename(self):
        iers.IERS_B.close()
        iers.IERS_B.open(iers.IERS_B_FILE)
        assert iers.IERS_B.iers_table is not None
        assert isinstance(iers.IERS_B.iers_table, QTable)
        iers.IERS_B.close()
        with pytest.raises(FILE_NOT_FOUND_ERROR):
            iers.IERS_B.open("surely this does not exist")

    def test_open_network_url(self):
        iers.IERS_A.close()
        iers.IERS_A.open(Path(IERS_A_EXCERPT).as_uri())
        assert iers.IERS_A.iers_table is not None
        assert isinstance(iers.IERS_A.iers_table, QTable)
        iers.IERS_A.close()


@pytest.mark.parametrize("path_transform", [os.fspath, Path])
def test_IERS_B_old_style_excerpt(path_transform):
    """Check that the instructions given in `IERS_B.read` actually work."""
    # If this test is changed, be sure to also adjust the instructions.
    #
    # TODO: this test and the note can probably be removed after
    # enough time has passed that old-style IERS_B files are simply
    # not around any more, say in 2025.  If so, also remove the excerpt
    # and the ReadMe.eopc04_IAU2000 file.
    old_style_file = path_transform(
        get_pkg_data_filename(os.path.join("data", "iers_b_old_style_excerpt"))
    )
    excerpt = iers.IERS_B.read(
        old_style_file,
        readme=get_pkg_data_filename(
            "data/ReadMe.eopc04_IAU2000", package="astropy.utils.iers"
        ),
        data_start=14,
    )
    assert isinstance(excerpt, QTable)
    assert "PM_x_dot" not in excerpt.colnames


class TestIERS_AExcerpt:
    @classmethod
    def teardown_class(cls):
        iers.IERS_A.close()

    def test_simple(self):
        # Test the IERS A reader. It is also a regression tests that ensures
        # values do not get overridden by IERS B; see #4933.
        iers_tab = iers.IERS_A.open(IERS_A_EXCERPT)

        assert (iers_tab["UT1_UTC"].unit / u.second).is_unity()
        assert "P" in iers_tab["UT1Flag"]
        assert "I" in iers_tab["UT1Flag"]
        assert "B" in iers_tab["UT1Flag"]
        assert np.all(
            (iers_tab["UT1Flag"] == "I")
            | (iers_tab["UT1Flag"] == "P")
            | (iers_tab["UT1Flag"] == "B")
        )

        assert (iers_tab["dX_2000A"].unit / u.marcsec).is_unity()
        assert (iers_tab["dY_2000A"].unit / u.marcsec).is_unity()
        assert "P" in iers_tab["NutFlag"]
        assert "I" in iers_tab["NutFlag"]
        assert "B" in iers_tab["NutFlag"]
        assert np.all(
            (iers_tab["NutFlag"] == "P")
            | (iers_tab["NutFlag"] == "I")
            | (iers_tab["NutFlag"] == "B")
        )

        assert (iers_tab["PM_x"].unit / u.arcsecond).is_unity()
        assert (iers_tab["PM_y"].unit / u.arcsecond).is_unity()
        assert "P" in iers_tab["PolPMFlag"]
        assert "I" in iers_tab["PolPMFlag"]
        assert "B" in iers_tab["PolPMFlag"]
        assert np.all(
            (iers_tab["PolPMFlag"] == "P")
            | (iers_tab["PolPMFlag"] == "I")
            | (iers_tab["PolPMFlag"] == "B")
        )

        t = Time([57053.0, 57054.0, 57055.0], format="mjd")
        ut1_utc, status = iers_tab.ut1_utc(t, return_status=True)
        assert status[0] == iers.FROM_IERS_B
        assert np.all(status[1:] == iers.FROM_IERS_A)
        # These values are *exactly* as given in the table, so they should
        # match to double precision accuracy.
        assert_quantity_allclose(
            ut1_utc, [-0.4916557, -0.4925323, -0.4934373] * u.s, atol=0.1 * u.ms
        )

        dcip_x, dcip_y, status = iers_tab.dcip_xy(t, return_status=True)
        assert status[0] == iers.FROM_IERS_B
        assert np.all(status[1:] == iers.FROM_IERS_A)
        # These values are *exactly* as given in the table, so they should
        # match to double precision accuracy.
        print(dcip_x)
        print(dcip_y)
        assert_quantity_allclose(
            dcip_x, [-0.086, -0.093, -0.087] * u.marcsec, atol=1.0 * u.narcsec
        )
        assert_quantity_allclose(
            dcip_y, [0.094, 0.081, 0.072] * u.marcsec, atol=1 * u.narcsec
        )

        pm_x, pm_y, status = iers_tab.pm_xy(t, return_status=True)
        assert status[0] == iers.FROM_IERS_B
        assert np.all(status[1:] == iers.FROM_IERS_A)
        assert_quantity_allclose(
            pm_x, [0.003734, 0.004581, 0.004623] * u.arcsec, atol=0.1 * u.marcsec
        )
        assert_quantity_allclose(
            pm_y, [0.310824, 0.313150, 0.315517] * u.arcsec, atol=0.1 * u.marcsec
        )

        # Table behaves properly as a table (e.g. can be sliced)
        assert len(iers_tab[:2]) == 2


class TestIERS_A:
    @classmethod
    def teardown_class(cls):
        iers.IERS_A.close()

    def test_simple(self):
        """Test that open() by default reads a 'finals2000A.all' file."""
        # Ensure we remove any cached table (gh-5131).
        iers.IERS_A.close()
        iers_tab = iers.IERS_A.open()
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5, 2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0.0, 0.5])
        ut1_utc, status = iers_tab.ut1_utc(jd1, jd2, return_status=True)
        assert np.all(status == iers.FROM_IERS_B)
        assert_quantity_allclose(
            ut1_utc,
            [-0.5868211, -0.5868184, -0.5868184, 0.4131816, 0.41328895] * u.s,
            atol=0.1 * u.ms,
        )
        ut1_utc2, status2 = iers_tab.ut1_utc(1e11, 0.0, return_status=True)
        assert status2 == iers.TIME_BEYOND_IERS_RANGE

        tnow = Time.now()

        ut1_utc3, status3 = iers_tab.ut1_utc(tnow, return_status=True)
        assert status3 == iers.FROM_IERS_A_PREDICTION
        assert ut1_utc3 != 0.0


class TestIERS_Auto:
    def setup_class(self):
        """Set up useful data for the tests."""
        self.N = 40
        self.ame = 30.0
        self.iers_a_file_1 = get_pkg_data_filename(
            os.path.join("data", "finals2000A-2016-02-30-test")
        )
        self.iers_a_file_2 = get_pkg_data_filename(
            os.path.join("data", "finals2000A-2016-04-30-test")
        )
        self.iers_a_url_1 = Path(self.iers_a_file_1).as_uri()
        self.iers_a_url_2 = Path(self.iers_a_file_2).as_uri()
        self.t = Time.now() + TimeDelta(10, format="jd") * np.arange(self.N)

        # This group of tests requires auto downloading to be on
        self._auto_download = iers.conf.auto_download
        iers.conf.auto_download = True

        # auto_download = False is tested in test_IERS_B_parameters_loading_into_IERS_Auto()

    def teardown_class(self):
        # Restore the auto downloading setting
        iers.conf.auto_download = self._auto_download

    def teardown_method(self, method):
        """Run this after every test."""
        iers.IERS_Auto.close()

    def test_interpolate_error_formatting(self):
        """Regression test: make sure the error message in
        IERS_Auto._check_interpolate_indices() is formatted correctly.
        """
        with iers.conf.set_temp("iers_auto_url", self.iers_a_url_1):
            with iers.conf.set_temp("iers_auto_url_mirror", self.iers_a_url_1):
                with iers.conf.set_temp("auto_max_age", self.ame):
                    with pytest.raises(
                        ValueError,
                        match=re.escape(iers.INTERPOLATE_ERROR.format(self.ame)),
                    ):
                        iers_table = iers.IERS_Auto.open()
                        with warnings.catch_warnings():
                            # Ignoring this if it comes up -- IERS_Auto predictive
                            # values are older than 30.0 days but downloading the
                            # latest table did not find newer values
                            warnings.simplefilter("ignore", iers.IERSStaleWarning)
                            iers_table.ut1_utc(self.t.jd1, self.t.jd2)

    def test_auto_max_age_none(self):
        """Make sure that iers.INTERPOLATE_ERROR's advice about setting
        auto_max_age = None actually works.
        """
        with iers.conf.set_temp("iers_auto_url", self.iers_a_url_1):
            with iers.conf.set_temp("auto_max_age", None):
                iers_table = iers.IERS_Auto.open()
                delta = iers_table.ut1_utc(self.t.jd1, self.t.jd2)
        assert isinstance(delta, np.ndarray)
        assert delta.shape == (self.N,)
        assert_quantity_allclose(delta, np.array([-0.2246227] * self.N) * u.s)

    def test_auto_max_age_minimum(self):
        """Check that the minimum auto_max_age is enforced."""
        with iers.conf.set_temp("iers_auto_url", self.iers_a_url_1):
            with iers.conf.set_temp("auto_max_age", 5.0):
                with pytest.raises(
                    ValueError,
                    match=(
                        r"IERS auto_max_age configuration value must be larger than 10"
                        r" days"
                    ),
                ):
                    iers_table = iers.IERS_Auto.open()
                    _ = iers_table.ut1_utc(self.t.jd1, self.t.jd2)

    def test_simple(self):
        with iers.conf.set_temp("iers_auto_url", self.iers_a_url_1):
            dat = iers.IERS_Auto.open()
            assert dat["MJD"][0] == 57359.0 * u.d
            assert dat["MJD"][-1] == 57539.0 * u.d

            # Pretend we are accessing at a time 7 days after start of predictive data
            predictive_mjd = dat.meta["predictive_mjd"]
            dat._time_now = Time(predictive_mjd, format="mjd") + 7 * u.d

            # Look at times before and after the test file begins.  0.1292934 is
            # the IERS-B value from MJD=57359.  The value in
            # finals2000A-2016-02-30-test has been replaced at this point.
            assert np.allclose(
                dat.ut1_utc(Time(50000, format="mjd").jd).value, 0.1292934
            )
            assert np.allclose(
                dat.ut1_utc(Time(60000, format="mjd").jd).value, -0.2246227
            )

            # Now pretend we are accessing at time 60 days after start of predictive data.
            # There will be a warning when downloading the file doesn't give new data
            # and an exception when extrapolating into the future with insufficient data.
            dat._time_now = Time(predictive_mjd, format="mjd") + 60 * u.d
            assert np.allclose(
                dat.ut1_utc(Time(50000, format="mjd").jd).value, 0.1292934
            )
            with (
                pytest.warns(
                    iers.IERSStaleWarning, match="IERS_Auto predictive values are older"
                ) as warns,
                pytest.raises(
                    ValueError,
                    match="interpolating from IERS_Auto using predictive values",
                ),
            ):
                dat.ut1_utc(Time(60000, format="mjd").jd)
            assert len(warns) == 1

            # Confirm that disabling the download means no warning because there is no
            # refresh to even fail, but there will still be the interpolation error
            with (
                iers.conf.set_temp("auto_download", False),
                pytest.raises(
                    ValueError,
                    match="interpolating from IERS_Auto using predictive values that are more",
                ),
            ):
                dat.ut1_utc(Time(60000, format="mjd").jd)

            # Warning only (i.e., no exception) if we are getting return status
            with pytest.warns(
                iers.IERSStaleWarning, match="IERS_Auto predictive values are older"
            ):
                dat.ut1_utc(Time(60000, format="mjd").jd, return_status=True)

            # Now set auto_max_age = None which says that we don't care how old the
            # available IERS-A file is.  There should be no warnings or exceptions.
            with iers.conf.set_temp("auto_max_age", None):
                dat.ut1_utc(Time(60000, format="mjd").jd)

        # Now point to a later file with same values but MJD increased by
        # 60 days and see that things work.  dat._time_now is still the same value
        # as before, i.e. right around the start of predictive values for the new file.
        # (In other words this is like downloading the latest file online right now).
        with iers.conf.set_temp("iers_auto_url", self.iers_a_url_2):
            # Look at times before and after the test file begins.  This forces a new download.
            assert np.allclose(
                dat.ut1_utc(Time(50000, format="mjd").jd).value, 0.1292934
            )
            assert np.allclose(dat.ut1_utc(Time(60000, format="mjd").jd).value, -0.3)

            # Now the time range should be different.
            assert dat["MJD"][0] == 57359.0 * u.d
            assert dat["MJD"][-1] == (57539.0 + 60) * u.d


@pytest.mark.parametrize("query", ["ut1_utc", "pm_xy"])
@pytest.mark.parametrize("jd", [np.array([]), Time([], format="mjd")])
@pytest.mark.parametrize("return_status", [False, True])
def test_empty_mjd(query, jd, return_status):
    # Regression test for gh-17008
    iers_table = iers.IERS_Auto.open()
    result = getattr(iers_table, query)(jd, return_status=return_status)
    n_exp = (1 if query == "ut1_utc" else 2) + (1 if return_status else 0)
    if n_exp == 1:
        assert isinstance(result, np.ndarray)
        assert result.size == 0
    else:
        assert len(result) == n_exp
        assert all(r.size == 0 for r in result)


def test_IERS_B_parameters_loading_into_IERS_Auto():
    # Make sure that auto downloading is off
    with iers.conf.set_temp("auto_download", False):
        A = iers.IERS_Auto.open()
    B = iers.IERS_B.open()

    ok_A = A["MJD"] <= B["MJD"][-1]
    assert not np.all(ok_A), "IERS B covers all of IERS A: should not happen"

    # We only overwrite IERS_B values in the IERS_A table that were already
    # there in the first place.  Better take that into account.
    ok_A &= np.isfinite(A["UT1_UTC_B"])

    i_B = np.searchsorted(B["MJD"], A["MJD"][ok_A])

    assert np.all(np.diff(i_B) == 1), "Valid region not contiguous"
    assert np.all(A["MJD"][ok_A] == B["MJD"][i_B])
    # Check that values are copied correctly.  Since units are not
    # necessarily the same, we use allclose with very strict tolerance.
    for name in ("UT1_UTC", "PM_x", "PM_y", "dX_2000A", "dY_2000A"):
        assert_quantity_allclose(
            A[name][ok_A],
            B[name][i_B],
            rtol=1e-15,
            err_msg=(
                f"Bug #9206 IERS B parameter {name} not copied over "
                "correctly to IERS Auto"
            ),
        )


# Issue with FTP, rework test into previous one when it's fixed
@pytest.mark.skipif(CI, reason="Flaky on CI")
@pytest.mark.remote_data
def test_iers_a_dl():
    iersa_tab = iers.IERS_A.open(iers.IERS_A_URL, cache=False)
    try:
        # some basic checks to ensure the format makes sense
        assert len(iersa_tab) > 0
        assert "UT1_UTC_A" in iersa_tab.colnames
    finally:
        iers.IERS_A.close()


@pytest.mark.remote_data
def test_iers_a_dl_mirror():
    iersa_tab = iers.IERS_A.open(iers.IERS_A_URL_MIRROR, cache=False)
    try:
        # some basic checks to ensure the format makes sense
        assert len(iersa_tab) > 0
        assert "UT1_UTC_A" in iersa_tab.colnames
    finally:
        iers.IERS_A.close()


@pytest.mark.remote_data
def test_iers_b_dl():
    iersb_tab = iers.IERS_B.open(iers.IERS_B_URL, cache=False)
    try:
        # some basic checks to ensure the format makes sense
        assert len(iersb_tab) > 0
        assert "UT1_UTC" in iersb_tab.colnames
    finally:
        iers.IERS_B.close()


def test_iers_b_out_of_range_handling():
    # The following error/warning applies only to IERS_B, not to the default IERS_Auto
    with iers.earth_orientation_table.set(iers.IERS_B.open()):
        now = Time.now()

        # Should be fine with bundled IERS-B
        (now - 300 * u.day).ut1

        # Default is to raise an error
        match = r"\(some\) times are outside of range covered by IERS table"
        with pytest.raises(iers.IERSRangeError, match=match):
            (now + 100 * u.day).ut1

        with iers.conf.set_temp("iers_degraded_accuracy", "warn"):
            with pytest.warns(iers.IERSDegradedAccuracyWarning, match=match):
                (now + 100 * u.day).ut1

        with iers.conf.set_temp("iers_degraded_accuracy", "ignore"):
            (now + 100 * u.day).ut1


@pytest.mark.remote_data
def test_iers_download_error_handling(tmp_path):
    # Make sure an IERS-A table isn't already loaded
    with set_temp_cache(tmp_path), iers.conf.set_temp("auto_download", True):
        iers.IERS_A.close()
        iers.IERS_Auto.close()
        iers.IERS.close()
        now = Time.now()

        # bad site name
        with iers.conf.set_temp("iers_auto_url", "FAIL FAIL"):
            # site that exists but doesn't have IERS data
            with iers.conf.set_temp("iers_auto_url_mirror", "https://google.com"):
                with pytest.warns(iers.IERSWarning) as record:
                    with iers.conf.set_temp("iers_degraded_accuracy", "ignore"):
                        (now + 400 * u.day).ut1

                assert len(record) == 3
                assert str(record[0].message).startswith(
                    "failed to download FAIL FAIL: Malformed URL"
                )
                assert str(record[1].message).startswith(
                    "malformed IERS table from https://google.com"
                )
                assert str(record[2].message).startswith(
                    "unable to download valid IERS file, using bundled IERS-A"
                )


OLD_DATA_FILES = {
    "Leap_Second.dat": IERS_LEAP_SECOND_FILE,
    "ReadMe.finals2000A": IERS_A_README,
    "ReadMe.eopc04": IERS_B_README,
    "eopc04.1962-now": IERS_B_FILE,
}


@pytest.mark.parametrize("data_file", sorted(OLD_DATA_FILES))
def test_get_pkg_data_filename_backcompat(data_file):
    # Check that get_pkg_data_filename continues to work without breakage
    # if users use it to access IERS tables and READMEs that used to be in
    # astropy/utils/iers/data.

    with pytest.warns(
        AstropyDeprecationWarning,
        match=f"Accessing {data_file} in this way is deprecated",
    ):
        filename = get_pkg_data_filename(
            "data/" + data_file, package="astropy.utils.iers"
        )

    assert filename == OLD_DATA_FILES[data_file]

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import numpy as np
import pytest

import astropy.units as u
from astropy.time import Time
from astropy.utils.iers import conf as iers_conf
from astropy.utils.iers import iers  # used in testing

allclose_jd = functools.partial(np.allclose, rtol=0, atol=1e-9)
allclose_sec = functools.partial(np.allclose, rtol=1e-15, atol=1e-4)
# 0.1 ms atol; IERS-B files change at that level.

try:
    iers.IERS_A.open()  # check if IERS_A is available
except OSError:
    HAS_IERS_A = False
else:
    HAS_IERS_A = True


def do_ut1_prediction_tst(iers_type):
    tnow = Time.now()
    iers_tab = iers_type.open()
    tnow.delta_ut1_utc, status = iers_tab.ut1_utc(tnow, return_status=True)
    assert status == iers.FROM_IERS_A_PREDICTION
    tnow_ut1_jd = tnow.ut1.jd
    assert allclose_jd(tnow_ut1_jd - tnow.jd, tnow.delta_ut1_utc / 86400)

    delta_ut1_utc = tnow.delta_ut1_utc
    with iers.earth_orientation_table.set(iers_type.open()):
        delta2, status2 = tnow.get_delta_ut1_utc(return_status=True)
        assert status2 == status
        assert delta2.to_value("s") == delta_ut1_utc

        tnow_ut1 = tnow.ut1
        assert tnow_ut1._delta_ut1_utc == delta_ut1_utc
        assert allclose_jd(tnow_ut1.jd - tnow.jd, tnow.delta_ut1_utc / 86400)


@pytest.mark.remote_data
class TestTimeUT1Remote:
    def setup_class(cls):
        # Need auto_download so that IERS_B won't be loaded and cause tests to
        # fail.
        iers_conf.auto_download = True

    def teardown_class(cls):
        # This setting is to be consistent with astropy/conftest.py
        iers_conf.auto_download = False

    def test_utc_to_ut1(self):
        "Test conversion of UTC to UT1, making sure to include a leap second"
        t = Time(
            [
                "2012-06-30 12:00:00",
                "2012-06-30 23:59:59",
                "2012-06-30 23:59:60",
                "2012-07-01 00:00:00",
                "2012-07-01 12:00:00",
            ],
            scale="utc",
        )
        t_ut1_jd = t.ut1.jd
        t_comp = np.array(
            [
                2456108.9999932079,
                2456109.4999816339,
                2456109.4999932083,
                2456109.5000047823,
                2456110.0000047833,
            ]
        )
        assert allclose_jd(t_ut1_jd, t_comp)
        t_back = t.ut1.utc
        assert allclose_jd(t.jd, t_back.jd)

        tnow = Time.now()

        tnow.ut1

    def test_ut1_iers_auto(self):
        do_ut1_prediction_tst(iers.IERS_Auto)


class TestTimeUT1:
    """Test Time.ut1 using IERS tables"""

    def test_ut1_to_utc(self):
        """Also test the reverse, around the leap second
        (round-trip test closes #2077)"""
        with iers_conf.set_temp("auto_download", False):
            t = Time(
                [
                    "2012-06-30 12:00:00",
                    "2012-06-30 23:59:59",
                    "2012-07-01 00:00:00",
                    "2012-07-01 00:00:01",
                    "2012-07-01 12:00:00",
                ],
                scale="ut1",
            )
            t_utc_jd = t.utc.jd
            t_comp = np.array(
                [
                    2456109.0000010049,
                    2456109.4999836441,
                    2456109.4999952177,
                    2456109.5000067917,
                    2456109.9999952167,
                ]
            )
            assert allclose_jd(t_utc_jd, t_comp)
            t_back = t.utc.ut1
            assert allclose_jd(t.jd, t_back.jd)

    def test_empty_ut1(self):
        """Testing for a zero-length Time object from UTC to UT1
        when an empty array is passed"""
        with iers_conf.set_temp("auto_download", False):
            t = Time(["2012-06-30 12:00:00"]) + np.arange(24) * u.hour
            t_empty = t[[]].ut1
            assert isinstance(t_empty, Time)
            assert t_empty.scale == "ut1"
            assert t_empty.size == 0

    @pytest.mark.parametrize(
        "jd1",
        [2441498.5, 2456108.5],  # Leap seconds of 1972-07-01 and 2012-07-01
    )
    @pytest.mark.parametrize("scale", ["tai", "utc", "tt", "tcb"])
    def test_ut1_conversion_at_leap_second(self, jd1, scale):
        """Regression test for gh-13517: right at a leap second, conversion
        from UT1 (then via UTC) could be off by a second due to rounding."""
        jd2 = np.array([0.9999999999999999, 1.0, 1.0000000000000002])
        # Use the bundled IERS-B table, which has data for 1972 as well.
        with iers.earth_orientation_table.set(iers.IERS_B.open()):
            t = Time(jd1, jd2, format="jd", scale="ut1")
            t2 = getattr(t, scale)
            assert t2[0] < t2[1] < t2[2]
            assert (t2[2] - t2[0]).sec < 1e-9
            # Round trip back to UT1.
            assert np.all(np.abs((t2.ut1 - t).sec) < 1e-9)
            # And the same for the other direction, i.e., starting with the
            # other scale, going to UT1 and back (the way back uses the
            # UT1 - UTC stored in the intermediate UT1 time, which needs care
            # within the leap second).
            t3 = Time(jd1, jd2, format="jd", scale=scale)
            t3_ut1 = t3.ut1
            assert t3_ut1[0] < t3_ut1[1] < t3_ut1[2]
            assert (t3_ut1[2] - t3_ut1[0]).sec < 1e-9
            assert np.all(np.abs((getattr(t3_ut1, scale) - t3).sec) < 1e-9)

    def test_ut1_tai_continuous_across_leap_second(self):
        """UT1 - TAI is continuous through a leap second, unlike UT1 - UTC
        and TAI - UTC (gh-13517)."""
        with iers_conf.set_temp("auto_download", False):
            t = Time("2012-07-01 00:00:00", scale="ut1") + np.arange(-20, 21) * (
                0.1 * u.s
            )
            tai = t.tai
            assert np.allclose((tai[1:] - tai[:-1]).sec, 0.1, rtol=0, atol=1e-8)
            # The leap second is present in UTC (UT1 - UTC is about -0.59 s
            # before it and +0.41 s after), i.e., UT1 - UTC jumps by 1 s.
            utc = t.utc
            assert utc.iso[20] == "2012-06-30 23:59:60.587"
            assert utc.iso[24] == "2012-06-30 23:59:60.987"
            assert utc.iso[25] == "2012-07-01 00:00:00.087"
            assert np.all(utc.delta_ut1_utc[:25] < 0)
            assert np.all(utc.delta_ut1_utc[25:] > 0)
            assert np.allclose(
                utc.delta_ut1_utc[25:] - utc.delta_ut1_utc[24], 1.0, atol=1e-6
            )
            assert np.all(np.abs((utc.ut1 - t).sec) < 1e-9)

    def test_ut1_to_utc_explicit_delta_at_leap_second(self):
        """With an explicitly set delta_ut1_utc, its sign is used to determine
        whether it is the value from before or after a nearby leap second
        (like in erfa.ut1utc): a positive leap second increases UT1 - UTC by
        one second, so UT1 - UTC is negative before and positive after it."""
        t = Time(
            ["2012-06-30 23:59:59", "2012-07-01 00:00:00", "2012-07-01 00:00:01"],
            scale="ut1",
        )
        # A non-negative value is the value from after the leap second, i.e.,
        # UT1 - UTC = -1 before it.
        t.delta_ut1_utc = 0.0
        assert np.all(
            t.utc.iso
            == [
                "2012-06-30 23:59:60.000",
                "2012-07-01 00:00:00.000",
                "2012-07-01 00:00:01.000",
            ]
        )
        # UT1 - UTC = -0.5 before the leap second and thus +0.5 after it, or,
        # equivalently, +0.5 after and -0.5 before, so both give the same.
        expected = [
            "2012-06-30 23:59:59.500",
            "2012-06-30 23:59:60.500",
            "2012-07-01 00:00:00.500",
        ]
        t.delta_ut1_utc = -0.5
        assert np.all(t.utc.iso == expected)
        t.delta_ut1_utc = 0.5
        assert np.all(t.utc.iso == expected)
        # Going to other scales gives consistent results.
        assert np.all(t.utc.tai.iso == t.tai.iso)
        assert np.all(t.tai.ut1.iso == t.iso)
        # UTC in and around the leap second, going to UT1 and back.
        t_utc = Time(
            ["2012-06-30 23:59:59.5", "2012-06-30 23:59:60.5", "2012-07-01 00:00:00.5"],
            scale="utc",
        )
        t_utc.delta_ut1_utc = -0.25
        assert np.all(
            t_utc.ut1.iso
            == [
                "2012-06-30 23:59:59.250",
                "2012-07-01 00:00:00.250",
                "2012-07-01 00:00:01.250",
            ]
        )
        assert np.all(t_utc.ut1.utc.iso == t_utc.iso)
        t_utc.delta_ut1_utc = 0.75  # Same as -0.25 before the leap second.
        assert np.all(t_utc.ut1.utc.iso == t_utc.iso)

    def test_delta_ut1_utc(self):
        """Accessing delta_ut1_utc should try to get it from IERS
        (closes #1924 partially)"""
        with iers_conf.set_temp("auto_download", False):
            t = Time("2012-06-30 12:00:00", scale="utc")
            assert not hasattr(t, "_delta_ut1_utc")
            # accessing delta_ut1_utc calculates it
            assert allclose_sec(t.delta_ut1_utc, -0.58682110003124965)
            # and keeps it around
            assert allclose_sec(t._delta_ut1_utc, -0.58682110003124965)


class TestTimeUT1SpecificIERSTable:
    @pytest.mark.skipif(not HAS_IERS_A, reason="requires IERS_A")
    def test_ut1_iers_A(self):
        do_ut1_prediction_tst(iers.IERS_A)

    def test_ut1_iers_B(self):
        tnow = Time.now()
        iers_b = iers.IERS_B.open()
        delta1, status1 = tnow.get_delta_ut1_utc(iers_b, return_status=True)
        assert status1 == iers.TIME_BEYOND_IERS_RANGE

        with iers.earth_orientation_table.set(iers.IERS_B.open()):
            delta2, status2 = tnow.get_delta_ut1_utc(return_status=True)
            assert status2 == status1

            with pytest.raises(iers.IERSRangeError):
                tnow.ut1

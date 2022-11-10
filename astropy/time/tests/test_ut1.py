# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import numpy as np
import pytest

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
    assert tnow_ut1_jd != tnow.jd

    delta_ut1_utc = tnow.delta_ut1_utc
    with iers.earth_orientation_table.set(iers_type.open()):
        delta2, status2 = tnow.get_delta_ut1_utc(return_status=True)
        assert status2 == status
        assert delta2.to_value("s") == delta_ut1_utc

        tnow_ut1 = tnow.ut1
        assert tnow_ut1._delta_ut1_utc == delta_ut1_utc
        assert tnow_ut1.jd != tnow.jd


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
        from astropy import units as u

        with iers_conf.set_temp("auto_download", False):
            t = Time(["2012-06-30 12:00:00"]) + np.arange(24) * u.hour
            t_empty = t[[]].ut1
            assert isinstance(t_empty, Time)
            assert t_empty.scale == "ut1"
            assert t_empty.size == 0

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

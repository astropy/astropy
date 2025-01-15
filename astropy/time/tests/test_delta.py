# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools
import itertools
import operator
from datetime import timedelta
from decimal import Decimal

import numpy as np
import pytest

from astropy import units as u
from astropy.table import Table
from astropy.time import (
    STANDARD_TIME_SCALES,
    TIME_DELTA_SCALES,
    TIME_SCALES,
    OperandTypeError,
    ScaleValueError,
    Time,
    TimeDelta,
    TimeDeltaMissingUnitWarning,
)
from astropy.utils import iers

allclose_jd = functools.partial(np.allclose, rtol=2.0**-52, atol=0)
allclose_jd2 = functools.partial(
    np.allclose, rtol=2.0**-52, atol=2.0**-52
)  # 20 ps atol
allclose_sec = functools.partial(
    np.allclose, rtol=2.0**-52, atol=2.0**-52 * 24 * 3600
)  # 20 ps atol
orig_auto_download = iers.conf.auto_download


def setup_module(module):
    """Use offline IERS table only."""
    iers.conf.auto_download = False


def teardown_module(module):
    """Restore original setting."""
    iers.conf.auto_download = orig_auto_download


class TestTimeDelta:
    """Test TimeDelta class"""

    def setup_method(self):
        self.t = Time("2010-01-01", scale="utc")
        self.t2 = Time("2010-01-02 00:00:01", scale="utc")
        self.t3 = Time(
            "2010-01-03 01:02:03",
            scale="utc",
            precision=9,
            in_subfmt="date_hms",
            out_subfmt="date_hm",
            location=(-75.0 * u.degree, 30.0 * u.degree, 500 * u.m),
        )
        self.t4 = Time("2010-01-01", scale="local")
        self.dt = TimeDelta(100.0, format="sec")
        self.dt_array = TimeDelta(np.arange(100, 1000, 100), format="sec")

    def test_sub(self):
        # time - time
        dt = self.t2 - self.t
        assert repr(dt).startswith(
            "<TimeDelta object: scale='tai' format='jd' value=1.00001157407"
        )
        assert allclose_jd(dt.jd, 86401.0 / 86400.0)
        assert allclose_sec(dt.sec, 86401.0)

        # time - delta_time
        t = self.t2 - dt
        assert t.iso == self.t.iso

        # delta_time - delta_time
        dt2 = dt - self.dt
        assert allclose_sec(dt2.sec, 86301.0)

        # delta_time - time
        with pytest.raises(OperandTypeError):
            dt - self.t

    def test_add(self):
        # time + time
        with pytest.raises(OperandTypeError):
            self.t2 + self.t

        # time + delta_time
        dt = self.t2 - self.t
        t2 = self.t + dt
        assert t2.iso == self.t2.iso

        # delta_time + delta_time
        dt2 = dt + self.dt
        assert allclose_sec(dt2.sec, 86501.0)

        # delta_time + time
        dt = self.t2 - self.t
        t2 = dt + self.t
        assert t2.iso == self.t2.iso

    def test_add_vector(self):
        """Check time arithmetic as well as properly keeping track of whether
        a time is a scalar or a vector"""
        t = Time(0.0, format="mjd", scale="tai")
        t2 = Time([0.0, 1.0], format="mjd", scale="tai")
        dt = TimeDelta(100.0, format="jd")
        dt2 = TimeDelta([100.0, 200.0], format="jd")

        out = t + dt
        assert allclose_jd(out.mjd, 100.0)
        assert out.isscalar

        out = t + dt2
        assert allclose_jd(out.mjd, [100.0, 200.0])
        assert not out.isscalar

        out = t2 + dt
        assert allclose_jd(out.mjd, [100.0, 101.0])
        assert not out.isscalar

        out = dt + dt
        assert allclose_jd(out.jd, 200.0)
        assert out.isscalar

        out = dt + dt2
        assert allclose_jd(out.jd, [200.0, 300.0])
        assert not out.isscalar

        # Reverse the argument order
        out = dt + t
        assert allclose_jd(out.mjd, 100.0)
        assert out.isscalar

        out = dt2 + t
        assert allclose_jd(out.mjd, [100.0, 200.0])
        assert not out.isscalar

        out = dt + t2
        assert allclose_jd(out.mjd, [100.0, 101.0])
        assert not out.isscalar

        out = dt2 + dt
        assert allclose_jd(out.jd, [200.0, 300.0])
        assert not out.isscalar

    def test_sub_vector(self):
        """Check time arithmetic as well as properly keeping track of whether
        a time is a scalar or a vector"""
        t = Time(0.0, format="mjd", scale="tai")
        t2 = Time([0.0, 1.0], format="mjd", scale="tai")
        dt = TimeDelta(100.0, format="jd")
        dt2 = TimeDelta([100.0, 200.0], format="jd")

        out = t - dt
        assert allclose_jd(out.mjd, -100.0)
        assert out.isscalar

        out = t - dt2
        assert allclose_jd(out.mjd, [-100.0, -200.0])
        assert not out.isscalar

        out = t2 - dt
        assert allclose_jd(out.mjd, [-100.0, -99.0])
        assert not out.isscalar

        out = dt - dt
        assert allclose_jd(out.jd, 0.0)
        assert out.isscalar

        out = dt - dt2
        assert allclose_jd(out.jd, [0.0, -100.0])
        assert not out.isscalar

    @pytest.mark.parametrize(
        "values", [(2455197.5, 2455198.5), ([2455197.5], [2455198.5])]
    )
    def test_copy_timedelta(self, values):
        """Test copying the values of a TimeDelta object by passing it into the
        Time initializer.
        """
        val1, val2 = values
        t = Time(val1, format="jd", scale="utc")
        t2 = Time(val2, format="jd", scale="utc")
        dt = t2 - t

        dt2 = TimeDelta(dt, copy=False)
        assert np.all(dt.jd == dt2.jd)
        assert dt._time.jd1 is dt2._time.jd1
        assert dt._time.jd2 is dt2._time.jd2

        dt2 = TimeDelta(dt, copy=True)
        assert np.all(dt.jd == dt2.jd)
        assert dt._time.jd1 is not dt2._time.jd1
        assert dt._time.jd2 is not dt2._time.jd2

        # Include initializers
        dt2 = TimeDelta(dt, format="sec")
        assert allclose_sec(dt2.value, 86400.0)

    def test_neg_abs(self):
        for dt in (self.dt, self.dt_array):
            dt2 = -dt
            assert np.all(dt2.jd == -dt.jd)
            dt3 = abs(dt)
            assert np.all(dt3.jd == dt.jd)
            dt4 = abs(dt2)
            assert np.all(dt4.jd == dt.jd)

    def test_mul_div(self):
        for dt in (self.dt, self.dt_array):
            dt2 = dt + dt + dt
            dt3 = 3.0 * dt
            assert allclose_jd(dt2.jd, dt3.jd)
            dt4 = dt3 / 3.0
            assert allclose_jd(dt4.jd, dt.jd)
        dt5 = self.dt * np.arange(3)
        assert dt5[0].jd == 0.0
        assert dt5[-1].jd == (self.dt + self.dt).jd
        dt6 = self.dt * [0, 1, 2]
        assert np.all(dt6.jd == dt5.jd)
        with pytest.raises(OperandTypeError):
            self.dt * self.t
        with pytest.raises(TypeError):
            self.dt * object()

    def test_mean(self):
        def is_consistent(time_delta: TimeDelta):
            mean_expected = (
                np.sum(time_delta.jd1) + np.sum(time_delta.jd2)
            ) / time_delta.size
            mean_test = time_delta.mean().jd1 + time_delta.mean().jd2
            return mean_test == mean_expected

        assert is_consistent(self.dt)
        assert is_consistent(self.dt_array)

    def test_keep_properties(self):
        # closes #1924 (partially)
        dt = TimeDelta(1000.0, format="sec")
        for t in (self.t, self.t3):
            ta = t + dt
            assert ta.location is t.location
            assert ta.precision == t.precision
            assert ta.in_subfmt == t.in_subfmt
            assert ta.out_subfmt == t.out_subfmt

            tr = dt + t
            assert tr.location is t.location
            assert tr.precision == t.precision
            assert tr.in_subfmt == t.in_subfmt
            assert tr.out_subfmt == t.out_subfmt

            ts = t - dt
            assert ts.location is t.location
            assert ts.precision == t.precision
            assert ts.in_subfmt == t.in_subfmt
            assert ts.out_subfmt == t.out_subfmt

        t_tdb = self.t.tdb
        assert hasattr(t_tdb, "_delta_tdb_tt")
        assert not hasattr(t_tdb, "_delta_ut1_utc")
        t_tdb_ut1 = t_tdb.ut1
        assert hasattr(t_tdb_ut1, "_delta_tdb_tt")
        assert hasattr(t_tdb_ut1, "_delta_ut1_utc")
        t_tdb_ut1_utc = t_tdb_ut1.utc
        assert hasattr(t_tdb_ut1_utc, "_delta_tdb_tt")
        assert hasattr(t_tdb_ut1_utc, "_delta_ut1_utc")
        # adding or subtracting some time should remove the delta's
        # since these are time-dependent and should be recalculated
        for op in (operator.add, operator.sub):
            t1 = op(t_tdb, dt)
            assert not hasattr(t1, "_delta_tdb_tt")
            assert not hasattr(t1, "_delta_ut1_utc")
            t2 = op(t_tdb_ut1, dt)
            assert not hasattr(t2, "_delta_tdb_tt")
            assert not hasattr(t2, "_delta_ut1_utc")
            t3 = op(t_tdb_ut1_utc, dt)
            assert not hasattr(t3, "_delta_tdb_tt")
            assert not hasattr(t3, "_delta_ut1_utc")

    def test_set_format(self):
        """
        Test basics of setting format attribute.
        """
        dt = TimeDelta(86400.0, format="sec")
        assert dt.value == 86400.0
        assert dt.format == "sec"

        dt.format = "jd"
        assert dt.value == 1.0
        assert dt.format == "jd"

        dt.format = "datetime"
        assert dt.value == timedelta(days=1)
        assert dt.format == "datetime"

    def test_from_non_float(self):
        dt = TimeDelta("1.000000000000001", format="jd")
        assert dt != TimeDelta(1.000000000000001, format="jd")  # precision loss.
        assert dt == TimeDelta(1, 0.000000000000001, format="jd")
        dt2 = TimeDelta(Decimal("1.000000000000001"), format="jd")
        assert dt2 == dt

    def test_to_value(self):
        dt = TimeDelta(86400.0, format="sec")
        assert dt.to_value("jd") == 1.0
        assert dt.to_value("jd", "str") == "1.0"
        assert dt.to_value("sec", subfmt="str") == "86400.0"
        with pytest.raises(
            ValueError,
            match="not one of the known formats.*failed to parse as a unit",
        ):
            dt.to_value("julian")

        with pytest.raises(TypeError, match="missing required format or unit"):
            dt.to_value()


class TestTimeDeltaScales:
    """Test scale conversion for Time Delta.
    Go through @taldcroft's list of expected behavior from #1932"""

    def setup_method(self):
        # pick a date that includes a leap second for better testing
        self.iso_times = [
            "2012-06-30 12:00:00",
            "2012-06-30 23:59:59",
            "2012-07-01 00:00:00",
            "2012-07-01 12:00:00",
        ]
        self.t = {
            scale: Time(self.iso_times, scale=scale, precision=9)
            for scale in TIME_SCALES
        }
        self.dt = {scale: self.t[scale] - self.t[scale][0] for scale in TIME_SCALES}

    def test_delta_scales_definition(self):
        for scale in list(TIME_DELTA_SCALES) + [None]:
            TimeDelta([0.0, 1.0, 10.0], format="sec", scale=scale)

        with pytest.raises(ScaleValueError):
            TimeDelta([0.0, 1.0, 10.0], format="sec", scale="utc")

    @pytest.mark.parametrize(
        ("scale1", "scale2"),
        list(itertools.product(STANDARD_TIME_SCALES, STANDARD_TIME_SCALES)),
    )
    def test_standard_scales_for_time_minus_time(self, scale1, scale2):
        """T(X) - T2(Y)  -- does T(X) - T2(Y).X and return dT(X)
        and T(X) +/- dT(Y)  -- does (in essence) (T(X).Y +/- dT(Y)).X

        I.e., time differences of two times should have the scale of the
        first time.  The one exception is UTC, which returns TAI.

        There are no standard timescales for which this does not work.
        """
        t1 = self.t[scale1]
        t2 = self.t[scale2]
        dt = t1 - t2
        if scale1 in TIME_DELTA_SCALES:
            assert dt.scale == scale1
        else:
            assert scale1 == "utc"
            assert dt.scale == "tai"

        # now check with delta time; also check reversibility
        t1_recover_t2_scale = t2 + dt
        assert t1_recover_t2_scale.scale == scale2
        t1_recover = getattr(t1_recover_t2_scale, scale1)
        assert allclose_jd(t1_recover.jd, t1.jd)
        t2_recover_t1_scale = t1 - dt
        assert t2_recover_t1_scale.scale == scale1
        t2_recover = getattr(t2_recover_t1_scale, scale2)
        assert allclose_jd(t2_recover.jd, t2.jd)

    def test_local_scales_for_time_minus_time(self):
        """T1(local) - T2(local) should return dT(local)
        T1(local) +/- dT(local) or T1(local) +/- Quantity(time-like) should
        also return T(local)

        I.e. Tests that time differences of two local scale times should
        return delta time with local timescale. Furthermore, checks that
        arithmetic of T(local) with dT(None) or time-like quantity does work.

        Also tests that subtracting two Time objects, one having local time
        scale and other having standard time scale should raise TypeError.
        """
        t1 = self.t["local"]
        t2 = Time("2010-01-01", scale="local")
        dt = t1 - t2
        assert dt.scale == "local"

        # now check with delta time
        t1_recover = t2 + dt
        assert t1_recover.scale == "local"
        assert allclose_jd(t1_recover.jd, t1.jd)
        # check that dT(None) can be subtracted from T(local)
        dt2 = TimeDelta([10.0], format="sec", scale=None)
        t3 = t2 - dt2
        assert t3.scale == t2.scale
        # check that time quantity can be subtracted from T(local)
        q = 10 * u.s
        assert (t2 - q).value == (t2 - dt2).value
        # Check that one cannot subtract/add times with a standard scale
        # from a local one (or vice versa)
        t1 = self.t["local"]
        for scale in STANDARD_TIME_SCALES:
            t2 = self.t[scale]
            with pytest.raises(TypeError):
                t1 - t2
            with pytest.raises(TypeError):
                t2 - t1
            with pytest.raises(TypeError):
                t2 - dt
            with pytest.raises(TypeError):
                t2 + dt
            with pytest.raises(TypeError):
                dt + t2

    def test_scales_for_delta_minus_delta(self):
        """dT(X) +/- dT2(Y) -- Add/substract JDs for dT(X) and dT(Y).X

        I.e. this will succeed if dT(Y) can be converted to scale X.
        Returns delta time in scale X
        """
        # geocentric timescales
        dt_tai = self.dt["tai"]
        dt_tt = self.dt["tt"]
        dt0 = dt_tai - dt_tt
        assert dt0.scale == "tai"
        # tai and tt have the same scale, so differences should be the same
        assert allclose_sec(dt0.sec, 0.0)

        dt_tcg = self.dt["tcg"]
        dt1 = dt_tai - dt_tcg
        assert dt1.scale == "tai"
        # tai and tcg do not have the same scale, so differences different
        assert not allclose_sec(dt1.sec, 0.0)

        t_tai_tcg = self.t["tai"].tcg
        dt_tai_tcg = t_tai_tcg - t_tai_tcg[0]
        dt2 = dt_tai - dt_tai_tcg
        assert dt2.scale == "tai"
        # but if tcg difference calculated from tai, it should roundtrip
        assert allclose_sec(dt2.sec, 0.0)
        # check that if we put TCG first, we get a TCG scale back
        dt3 = dt_tai_tcg - dt_tai
        assert dt3.scale == "tcg"
        assert allclose_sec(dt3.sec, 0.0)

        for scale in "tdb", "tcb", "ut1":
            with pytest.raises(TypeError):
                dt_tai - self.dt[scale]

        # barycentric timescales
        dt_tcb = self.dt["tcb"]
        dt_tdb = self.dt["tdb"]
        dt4 = dt_tcb - dt_tdb
        assert dt4.scale == "tcb"
        assert not allclose_sec(dt1.sec, 0.0)

        t_tcb_tdb = self.t["tcb"].tdb
        dt_tcb_tdb = t_tcb_tdb - t_tcb_tdb[0]
        dt5 = dt_tcb - dt_tcb_tdb
        assert dt5.scale == "tcb"
        assert allclose_sec(dt5.sec, 0.0)

        for scale in "utc", "tai", "tt", "tcg", "ut1":
            with pytest.raises(TypeError):
                dt_tcb - self.dt[scale]

        # rotational timescale
        dt_ut1 = self.dt["ut1"]
        dt5 = dt_ut1 - dt_ut1[-1]
        assert dt5.scale == "ut1"
        assert dt5[-1].sec == 0.0

        for scale in "utc", "tai", "tt", "tcg", "tcb", "tdb":
            with pytest.raises(TypeError):
                dt_ut1 - self.dt[scale]

        # local time scale
        dt_local = self.dt["local"]
        dt6 = dt_local - dt_local[-1]
        assert dt6.scale == "local"
        assert dt6[-1].sec == 0.0

        for scale in "utc", "tai", "tt", "tcg", "tcb", "tdb", "ut1":
            with pytest.raises(TypeError):
                dt_local - self.dt[scale]

    @pytest.mark.parametrize(
        ("scale", "op"),
        list(itertools.product(TIME_SCALES, (operator.add, operator.sub))),
    )
    def test_scales_for_delta_scale_is_none(self, scale, op):
        """T(X) +/- dT(None) or T(X) +/- Quantity(time-like)

        This is always allowed and just adds JDs, i.e., the scale of
        the TimeDelta or time-like Quantity will be taken to be X.
        The one exception is again for X=UTC, where TAI is assumed instead,
        so that a day is always defined as 86400 seconds.
        """
        dt_none = TimeDelta([0.0, 1.0, -1.0, 1000.0], format="sec")
        assert dt_none.scale is None
        q_time = dt_none.to("s")

        dt = self.dt[scale]
        dt1 = op(dt, dt_none)
        assert dt1.scale == dt.scale
        assert allclose_jd(dt1.jd, op(dt.jd, dt_none.jd))
        dt2 = op(dt_none, dt)
        assert dt2.scale == dt.scale
        assert allclose_jd(dt2.jd, op(dt_none.jd, dt.jd))
        dt3 = op(q_time, dt)
        assert dt3.scale == dt.scale
        assert allclose_jd(dt3.jd, dt2.jd)

        t = self.t[scale]
        t1 = op(t, dt_none)
        assert t1.scale == t.scale
        assert allclose_jd(t1.jd, op(t.jd, dt_none.jd))
        if op is operator.add:
            t2 = op(dt_none, t)
            assert t2.scale == t.scale
            assert allclose_jd(t2.jd, t1.jd)
        t3 = op(t, q_time)
        assert t3.scale == t.scale
        assert allclose_jd(t3.jd, t1.jd)

    @pytest.mark.parametrize("scale", TIME_SCALES)
    def test_delta_day_is_86400_seconds(self, scale):
        """TimeDelta or Quantity holding 1 day always means 24*60*60 seconds

        This holds true for all timescales but UTC, for which leap-second
        days are longer or shorter by one second.
        """
        t = self.t[scale]
        dt_day = TimeDelta(1.0, format="jd")
        q_day = dt_day.to("day")

        dt_day_leap = t[-1] - t[0]
        # ^ = exclusive or, so either equal and not UTC, or not equal and UTC
        assert allclose_jd(dt_day_leap.jd, dt_day.jd) ^ (scale == "utc")

        t1 = t[0] + dt_day
        assert allclose_jd(t1.jd, t[-1].jd) ^ (scale == "utc")
        t2 = q_day + t[0]
        assert allclose_jd(t2.jd, t[-1].jd) ^ (scale == "utc")
        t3 = t[-1] - dt_day
        assert allclose_jd(t3.jd, t[0].jd) ^ (scale == "utc")
        t4 = t[-1] - q_day
        assert allclose_jd(t4.jd, t[0].jd) ^ (scale == "utc")


def test_timedelta_setitem():
    t = TimeDelta([1, 2, 3] * u.d, format="jd")

    t[0] = 0.5
    assert allclose_jd(t.value, [0.5, 2, 3])

    t[1:] = 4.5
    assert allclose_jd(t.value, [0.5, 4.5, 4.5])

    t[:] = 86400 * u.s
    assert allclose_jd(t.value, [1, 1, 1])

    t[1] = TimeDelta(2, format="jd")
    assert allclose_jd(t.value, [1, 2, 1])

    with pytest.raises(ValueError) as err:
        t[1] = 1 * u.m
    assert "cannot convert value to a compatible TimeDelta" in str(err.value)


def test_timedelta_setitem_sec():
    t = TimeDelta([1, 2, 3], format="sec")

    t[0] = 0.5
    assert allclose_jd(t.value, [0.5, 2, 3])

    t[1:] = 4.5
    assert allclose_jd(t.value, [0.5, 4.5, 4.5])

    t[:] = 1 * u.day
    assert allclose_jd(t.value, [86400, 86400, 86400])

    t[1] = TimeDelta(2, format="jd")
    assert allclose_jd(t.value, [86400, 86400 * 2, 86400])

    with pytest.raises(ValueError) as err:
        t[1] = 1 * u.m
    assert "cannot convert value to a compatible TimeDelta" in str(err.value)


def test_timedelta_mask():
    t = TimeDelta([1, 2] * u.d, format="jd")
    t[1] = np.ma.masked
    assert np.all(t.mask == [False, True])
    assert allclose_jd(t[0].value, 1)
    assert t.value[1] is np.ma.masked


def test_python_timedelta_scalar():
    td = timedelta(days=1, seconds=1)
    td1 = TimeDelta(td, format="datetime")

    assert td1.sec == 86401.0

    td2 = TimeDelta(86401.0, format="sec")
    assert td2.datetime == td


def test_python_timedelta_vector():
    td = [
        [timedelta(days=1), timedelta(days=2)],
        [timedelta(days=3), timedelta(days=4)],
    ]

    td1 = TimeDelta(td, format="datetime")

    assert np.all(td1.jd == [[1, 2], [3, 4]])

    td2 = TimeDelta([[1, 2], [3, 4]], format="jd")
    assert np.all(td2.datetime == td)


def test_timedelta_to_datetime():
    td = TimeDelta(1, format="jd")

    assert td.to_datetime() == timedelta(days=1)

    td2 = TimeDelta([[1, 2], [3, 4]], format="jd")
    td = [
        [timedelta(days=1), timedelta(days=2)],
        [timedelta(days=3), timedelta(days=4)],
    ]

    assert np.all(td2.to_datetime() == td)


def test_insert_timedelta():
    tm = TimeDelta([1, 2], format="sec")

    # Insert a scalar using an auto-parsed string
    tm2 = tm.insert(1, TimeDelta([10, 20], format="sec"))
    assert np.all(tm2 == TimeDelta([1, 10, 20, 2], format="sec"))


def test_no_units_warning():
    with pytest.warns(TimeDeltaMissingUnitWarning):
        delta = TimeDelta(1)
        assert delta.to_value(u.day) == 1

    with pytest.warns(TimeDeltaMissingUnitWarning):
        table = Table({"t": [1, 2, 3]})
        delta = TimeDelta(table["t"])
        assert np.all(delta.to_value(u.day) == [1, 2, 3])

    with pytest.warns(TimeDeltaMissingUnitWarning):
        delta = TimeDelta(np.array([1, 2, 3]))
        assert np.all(delta.to_value(u.day) == [1, 2, 3])

    with pytest.warns(TimeDeltaMissingUnitWarning):
        t = Time("2012-01-01") + 1
        assert t.isot[:10] == "2012-01-02"

    with pytest.warns(TimeDeltaMissingUnitWarning):
        comp = TimeDelta([1, 2, 3], format="jd") >= 2
        assert np.all(comp == [False, True, True])

    with pytest.warns(TimeDeltaMissingUnitWarning):
        # 2 is also interpreted as days, not seconds
        assert (TimeDelta(5 * u.s) > 2) is False

    # with unit is ok
    assert TimeDelta(1 * u.s).to_value(u.s) == 1

    # with format is also ok
    assert TimeDelta(1, format="sec").to_value(u.s) == 1
    assert TimeDelta(1, format="jd").to_value(u.day) == 1

    # table column with units
    table = Table({"t": [1, 2, 3] * u.s})
    assert np.all(TimeDelta(table["t"]).to_value(u.s) == [1, 2, 3])

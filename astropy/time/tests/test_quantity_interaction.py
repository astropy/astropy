# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import numpy as np
import pytest

from astropy import units as u
from astropy.table import Column
from astropy.time import Time, TimeDelta

allclose_sec = functools.partial(
    np.allclose, rtol=2.0**-52, atol=2.0**-52 * 24 * 3600
)  # 20 ps atol


class TestTimeQuantity:
    """Test Interaction of Time with Quantities"""

    def test_valid_quantity_input(self):
        """Test Time formats that are allowed to take quantity input."""
        q = 2450000.125 * u.day
        t1 = Time(q, format="jd", scale="utc")
        assert t1.value == q.value
        q2 = q.to(u.second)
        t2 = Time(q2, format="jd", scale="utc")
        assert t2.value == q.value == q2.to_value(u.day)
        q3 = q - 2400000.5 * u.day
        t3 = Time(q3, format="mjd", scale="utc")
        assert t3.value == q3.value
        # test we can deal with two quantity arguments, with different units
        qs = 24.0 * 36.0 * u.second
        t4 = Time(q3, qs, format="mjd", scale="utc")
        assert t4.value == (q3 + qs).to_value(u.day)

        qy = 1990.0 * u.yr
        ty1 = Time(qy, format="jyear", scale="utc")
        assert ty1.value == qy.value
        ty2 = Time(qy.to(u.day), format="jyear", scale="utc")
        assert ty2.value == qy.value
        qy2 = 10.0 * u.yr
        tcxc = Time(qy2, format="cxcsec")
        assert tcxc.value == qy2.to_value(u.second)
        tgps = Time(qy2, format="gps")
        assert tgps.value == qy2.to_value(u.second)
        tunix = Time(qy2, format="unix")
        assert tunix.value == qy2.to_value(u.second)
        qd = 2000.0 * 365.0 * u.day
        tplt = Time(qd, format="plot_date", scale="utc")
        assert tplt.value == qd.value

    def test_invalid_quantity_input(self):
        with pytest.raises(u.UnitsError):
            Time(2450000.0 * u.m, format="jd", scale="utc")

        with pytest.raises(u.UnitsError):
            Time(2450000.0 * u.dimensionless_unscaled, format="jd", scale="utc")

    def test_column_with_and_without_units(self):
        """Ensure a Column without a unit is treated as an array [#3648]"""
        a = np.arange(50000.0, 50010.0)
        ta = Time(a, format="mjd")
        c1 = Column(np.arange(50000.0, 50010.0), name="mjd")
        tc1 = Time(c1, format="mjd")
        assert np.all(ta == tc1)
        c2 = Column(np.arange(50000.0, 50010.0), name="mjd", unit="day")
        tc2 = Time(c2, format="mjd")
        assert np.all(ta == tc2)
        c3 = Column(np.arange(50000.0, 50010.0), name="mjd", unit="m")
        with pytest.raises(u.UnitsError):
            Time(c3, format="mjd")

    def test_no_quantity_input_allowed(self):
        """Time formats that are not allowed to take Quantity input."""
        qy = 1990.0 * u.yr
        for fmt in ("iso", "yday", "datetime", "byear", "byear_str", "jyear_str"):
            with pytest.raises(ValueError):
                Time(qy, format=fmt, scale="utc")

    def test_valid_quantity_operations(self):
        """Check that adding a time-valued quantity to a Time gives a Time"""
        t0 = Time(100000.0, format="cxcsec")
        q1 = 10.0 * u.second
        t1 = t0 + q1
        assert isinstance(t1, Time)
        assert t1.value == t0.value + q1.to_value(u.second)
        q2 = 1.0 * u.day
        t2 = t0 - q2
        assert allclose_sec(t2.value, t0.value - q2.to_value(u.second))
        # check broadcasting
        q3 = np.arange(15.0).reshape(3, 5) * u.hour
        t3 = t0 - q3
        assert t3.shape == q3.shape
        assert allclose_sec(t3.value, t0.value - q3.to_value(u.second))

    def test_invalid_quantity_operations(self):
        """Check that comparisons of Time with quantities does not work
        (even for time-like, since we cannot compare Time to TimeDelta)"""
        with pytest.raises(TypeError):
            Time(100000.0, format="cxcsec") > 10.0 * u.m  # noqa: B015
        with pytest.raises(TypeError):
            Time(100000.0, format="cxcsec") > 10.0 * u.second  # noqa: B015


class TestTimeDeltaQuantity:
    """Test interaction of TimeDelta with Quantities"""

    def test_valid_quantity_input(self):
        """Test that TimeDelta can take quantity input."""
        q = 500.25 * u.day
        dt1 = TimeDelta(q, format="jd")
        assert dt1.value == q.value
        dt2 = TimeDelta(q, format="sec")
        assert dt2.value == q.to_value(u.second)
        dt3 = TimeDelta(q)
        assert dt3.value == q.value

    def test_invalid_quantity_input(self):
        with pytest.raises(u.UnitsError):
            TimeDelta(2450000.0 * u.m, format="jd")

        with pytest.raises(u.UnitsError):
            Time(2450000.0 * u.dimensionless_unscaled, format="jd", scale="utc")

        with pytest.raises(TypeError):
            TimeDelta(100, format="sec") > 10.0 * u.m  # noqa: B015

    def test_quantity_output(self):
        q = 500.25 * u.day
        dt = TimeDelta(q)
        assert dt.to(u.day) == q
        assert dt.to_value(u.day) == q.value
        assert dt.to_value("day") == q.value
        assert dt.to(u.second).value == q.to_value(u.second)
        assert dt.to_value(u.second) == q.to_value(u.second)
        assert dt.to_value("s") == q.to_value(u.second)
        # Following goes through "format", but should be the same.
        assert dt.to_value("sec") == q.to_value(u.second)

    def test_quantity_output_errors(self):
        dt = TimeDelta(250.0, format="sec")
        with pytest.raises(u.UnitsError):
            dt.to(u.m)
        with pytest.raises(u.UnitsError):
            dt.to_value(u.m)
        with pytest.raises(u.UnitsError):
            dt.to_value(unit=u.m)
        with pytest.raises(
            ValueError,
            match="not one of the known formats.*failed to parse as a unit",
        ):
            dt.to_value("parrot")
        with pytest.raises(TypeError):
            dt.to_value("sec", unit=u.s)
        with pytest.raises(
            ValueError,
            match=r"cannot specify 'subfmt' and positional arg.*not a valid format",
        ):
            dt.to_value(u.s, subfmt="str")

    def test_valid_quantity_operations1(self):
        """Check adding/subtracting/comparing a time-valued quantity works
        with a TimeDelta.  Addition/subtraction should give TimeDelta"""
        t0 = TimeDelta(106400.0, format="sec")
        q1 = 10.0 * u.second
        t1 = t0 + q1
        assert isinstance(t1, TimeDelta)
        assert t1.value == t0.value + q1.to_value(u.second)
        q2 = 1.0 * u.day
        t2 = t0 - q2
        assert isinstance(t2, TimeDelta)
        assert allclose_sec(t2.value, t0.value - q2.to_value(u.second))
        # now comparisons
        assert t0 > q1
        assert t0 < 1.0 * u.yr
        # and broadcasting
        q3 = np.arange(12.0).reshape(4, 3) * u.hour
        t3 = t0 + q3
        assert isinstance(t3, TimeDelta)
        assert t3.shape == q3.shape
        assert allclose_sec(t3.value, t0.value + q3.to_value(u.second))

    def test_valid_quantity_operations2(self):
        """Check that TimeDelta is treated as a quantity where possible."""
        t0 = TimeDelta(100000.0, format="sec")
        f = 1.0 / t0
        assert isinstance(f, u.Quantity)
        assert f.unit == 1.0 / u.day
        g = 10.0 * u.m / u.second**2
        v = t0 * g
        assert isinstance(v, u.Quantity)
        assert u.allclose(v, t0.sec * g.value * u.m / u.second)
        q = np.log10(t0 / u.second)
        assert isinstance(q, u.Quantity)
        assert q.value == np.log10(t0.sec)
        s = 1.0 * u.m
        v = s / t0
        assert isinstance(v, u.Quantity)
        assert u.allclose(v, 1.0 / t0.sec * u.m / u.s)
        t = 1.0 * u.s
        t2 = t0 * t
        assert isinstance(t2, u.Quantity)
        assert u.allclose(t2, t0.sec * u.s**2)
        t3 = [1] / t0
        assert isinstance(t3, u.Quantity)
        assert u.allclose(t3, 1 / (t0.sec * u.s))
        # broadcasting
        t1 = TimeDelta(np.arange(100000.0, 100012.0).reshape(6, 2), format="sec")
        f = np.array([1.0, 2.0]) * u.cycle * u.Hz
        phase = f * t1
        assert isinstance(phase, u.Quantity)
        assert phase.shape == t1.shape
        assert u.allclose(phase, t1.sec * f.value * u.cycle)
        q = t0 * t1
        assert isinstance(q, u.Quantity)
        assert np.all(q == t0.to(u.day) * t1.to(u.day))
        q = t1 / t0
        assert isinstance(q, u.Quantity)
        assert np.all(q == t1.to(u.day) / t0.to(u.day))

    def test_valid_quantity_operations3(self):
        """Test a TimeDelta remains one if possible."""
        t0 = TimeDelta(10.0, format="jd")
        q = 10.0 * u.one
        t1 = q * t0
        assert isinstance(t1, TimeDelta)
        assert t1 == TimeDelta(100.0, format="jd")
        t2 = t0 * q
        assert isinstance(t2, TimeDelta)
        assert t2 == TimeDelta(100.0, format="jd")
        t3 = t0 / q
        assert isinstance(t3, TimeDelta)
        assert t3 == TimeDelta(1.0, format="jd")
        q2 = 1.0 * u.percent
        t4 = t0 * q2
        assert isinstance(t4, TimeDelta)
        assert abs(t4 - TimeDelta(0.1, format="jd")) < 1.0 * u.ns
        q3 = 1.0 * u.hr / (36.0 * u.s)
        t5 = q3 * t0
        assert isinstance(t4, TimeDelta)
        assert abs(t5 - TimeDelta(1000.0, format="jd")) < 1.0 * u.ns
        # Test multiplication with a unit.
        t6 = t0 * u.one
        assert isinstance(t6, TimeDelta)
        assert t6 == TimeDelta(10.0, format="jd")
        t7 = u.one * t0
        assert isinstance(t7, TimeDelta)
        assert t7 == TimeDelta(10.0, format="jd")
        t8 = t0 * ""
        assert isinstance(t8, TimeDelta)
        assert t8 == TimeDelta(10.0, format="jd")
        t9 = "" * t0
        assert isinstance(t9, TimeDelta)
        assert t9 == TimeDelta(10.0, format="jd")
        t10 = t0 / u.one
        assert isinstance(t10, TimeDelta)
        assert t6 == TimeDelta(10.0, format="jd")
        t11 = t0 / ""
        assert isinstance(t11, TimeDelta)
        assert t11 == TimeDelta(10.0, format="jd")
        t12 = t0 / [1]
        assert isinstance(t12, TimeDelta)
        assert t12 == TimeDelta(10.0, format="jd")
        t13 = [1] * t0
        assert isinstance(t13, TimeDelta)
        assert t13 == TimeDelta(10.0, format="jd")

    def test_invalid_quantity_operations(self):
        """Check comparisons of TimeDelta with non-time quantities fails."""
        with pytest.raises(TypeError):
            TimeDelta(100000.0, format="sec") > 10.0 * u.m  # noqa: B015

    def test_invalid_quantity_operations2(self):
        """Check that operations with non-time/quantity fail."""
        td = TimeDelta(100000.0, format="sec")
        with pytest.raises(TypeError):
            td * object()
        with pytest.raises(TypeError):
            td / object()

    def test_invalid_quantity_broadcast(self):
        """Check broadcasting rules in interactions with Quantity."""
        t0 = TimeDelta(np.arange(12.0).reshape(4, 3), format="sec")
        with pytest.raises(ValueError):
            t0 + np.arange(4.0) * u.s


class TestDeltaAttributes:
    def test_delta_ut1_utc(self):
        t = Time("2010-01-01 00:00:00", format="iso", scale="utc", precision=6)
        t.delta_ut1_utc = 0.3 * u.s
        assert t.ut1.iso == "2010-01-01 00:00:00.300000"
        t.delta_ut1_utc = 0.4 / 60.0 * u.minute
        assert t.ut1.iso == "2010-01-01 00:00:00.400000"
        with pytest.raises(u.UnitsError):
            t.delta_ut1_utc = 0.4 * u.m
        # Also check that a TimeDelta works.
        t.delta_ut1_utc = TimeDelta(0.3, format="sec")
        assert t.ut1.iso == "2010-01-01 00:00:00.300000"
        t.delta_ut1_utc = TimeDelta(0.5 / 24.0 / 3600.0, format="jd")
        assert t.ut1.iso == "2010-01-01 00:00:00.500000"

    def test_delta_tdb_tt(self):
        t = Time("2010-01-01 00:00:00", format="iso", scale="tt", precision=6)
        t.delta_tdb_tt = 20.0 * u.second
        assert t.tdb.iso == "2010-01-01 00:00:20.000000"
        t.delta_tdb_tt = 30.0 / 60.0 * u.minute
        assert t.tdb.iso == "2010-01-01 00:00:30.000000"
        with pytest.raises(u.UnitsError):
            t.delta_tdb_tt = 0.4 * u.m
        # Also check that a TimeDelta works.
        t.delta_tdb_tt = TimeDelta(40.0, format="sec")
        assert t.tdb.iso == "2010-01-01 00:00:40.000000"
        t.delta_tdb_tt = TimeDelta(50.0 / 24.0 / 3600.0, format="jd")
        assert t.tdb.iso == "2010-01-01 00:00:50.000000"


@pytest.mark.parametrize(
    "q1, q2",
    (
        (5e8 * u.s, None),
        (5e17 * u.ns, None),
        (4e8 * u.s, 1e17 * u.ns),
        (4e14 * u.us, 1e17 * u.ns),
    ),
)
def test_quantity_conversion_rounding(q1, q2):
    """Check that no rounding errors are incurred by unit conversion.

    This occurred before as quantities in seconds were converted to days
    before trying to split them into two-part doubles.  See gh-7622.
    """
    t = Time("2001-01-01T00:00:00.", scale="tai")
    expected = Time("2016-11-05T00:53:20.", scale="tai")
    if q2 is None:
        t0 = t + q1
    else:
        t0 = t + q1 + q2
    assert abs(t0 - expected) < 20 * u.ps
    dt1 = TimeDelta(q1, q2)
    t1 = t + dt1
    assert abs(t1 - expected) < 20 * u.ps
    dt2 = TimeDelta(q1, q2, format="sec")
    t2 = t + dt2
    assert abs(t2 - expected) < 20 * u.ps

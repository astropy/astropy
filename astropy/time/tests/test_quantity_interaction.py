# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import numpy as np

from ...tests.helper import pytest
from .. import Time, TimeDelta, OperandTypeError
from ... import units as u

allclose_sec = functools.partial(np.allclose, rtol=2. ** -52,
                                 atol=2. ** -52 * 24 * 3600)  # 20 ps atol


class TestTimeQuantity():
    """Test Interaction of Time with Quantities"""

    def test_valid_quantity_input(self):
        """Test Time formats that are allowed to take quantity input."""
        q = 2450000.125*u.day
        t1 = Time(q, format='jd', scale='utc')
        assert t1.value == q.value
        q2 = q.to(u.second)
        t2 = Time(q2, format='jd', scale='utc')
        assert t2.value == q.value == q2.to(u.day).value
        q3 = q-2400000.5*u.day
        t3 = Time(q3, format='mjd', scale='utc')
        assert t3.value == q3.value
        # test we can deal with two quantity arguments, with different units
        qs = 24.*36.*u.second
        t4 = Time(q3, qs, format='mjd', scale='utc')
        assert t4.value == (q3+qs).to(u.day).value

        qy = 1990.*u.yr
        ty1 = Time(qy, format='jyear', scale='utc')
        assert ty1.value == qy.value
        ty2 = Time(qy.to(u.day), format='jyear', scale='utc')
        assert ty2.value == qy.value
        qy2 = 10.*u.yr
        tcxc = Time(qy2, format='cxcsec')
        assert tcxc.value == qy2.to(u.second).value
        tgps = Time(qy2, format='gps')
        assert tgps.value == qy2.to(u.second).value
        tunix = Time(qy2, format='unix')
        assert tunix.value == qy2.to(u.second).value
        qd = 2000.*365.*u.day
        tplt = Time(qd, format='plot_date', scale='utc')
        assert tplt.value == qd.value

    def test_invalid_quantity_input(self):
        with pytest.raises(u.UnitsError):
            Time(2450000.*u.m, format='jd', scale='utc')

        with pytest.raises(u.UnitsError):
            Time(2450000.*u.dimensionless_unscaled, format='jd', scale='utc')

    def test_no_quantity_input_allowed(self):
        """Time formats that are not allowed to take Quantity input."""
        qy = 1990.*u.yr
        for fmt in ('iso', 'yday', 'datetime', 'byear',
                    'byear_str', 'jyear_str'):
            with pytest.raises(ValueError):
                Time(qy, format=fmt, scale='utc')

    def test_valid_quantity_operations(self):
        """Check that adding a time-valued quantity to a Time gives a Time"""
        t0 = Time(100000., format='cxcsec')
        q1 = 10.*u.second
        t1 = t0 + q1
        assert isinstance(t1, Time)
        assert t1.value == t0.value+q1.to(u.second).value
        q2 = 1.*u.day
        t2 = t0 - q2
        assert allclose_sec(t2.value, t0.value-q2.to(u.second).value)

    def test_invalid_quantity_operations(self):
        """Check that comparisons of Time with quantities does not work
        (even for time-like, since we cannot compare Time to TimeDelta)"""
        with pytest.raises(OperandTypeError):
            Time(100000., format='cxcsec') > 10.*u.m
        with pytest.raises(OperandTypeError):
            Time(100000., format='cxcsec') > 10.*u.second


class TestTimeDeltaQuantity():
    """Test interaction of TimeDelta with Quantities"""
    def test_valid_quantity_input(self):
        """Test that TimeDelta can take quantity input."""
        q = 500.25*u.day
        dt1 = TimeDelta(q, format='jd')
        assert dt1.value == q.value
        dt2 = TimeDelta(q, format='sec')
        assert dt2.value == q.to(u.second).value
        dt3 = TimeDelta(q)
        assert dt3.value == q.value

    def test_invalid_quantity_input(self):
        with pytest.raises(u.UnitsError):
            TimeDelta(2450000.*u.m, format='jd', scale='utc')

        with pytest.raises(u.UnitsError):
            Time(2450000.*u.dimensionless_unscaled, format='jd', scale='utc')

        with pytest.raises(OperandTypeError):
            TimeDelta(100, format='sec') > 10.*u.m

    def test_quantity_output(self):
        q = 500.25*u.day
        dt = TimeDelta(q)
        assert dt.to(u.day) == q
        assert dt.to(u.second).value == q.to(u.second).value
        with pytest.raises(u.UnitsError):
            dt.to(u.m)

    def test_valid_quantity_operations1(self):
        """Check adding/substracting/comparing a time-valued quantity works
        with a TimeDelta.  Addition/subtraction should give TimeDelta"""
        t0 = TimeDelta(106400., format='sec')
        q1 = 10.*u.second
        t1 = t0 + q1
        assert isinstance(t1, TimeDelta)
        assert t1.value == t0.value+q1.to(u.second).value
        q2 = 1.*u.day
        t2 = t0 - q2
        assert allclose_sec(t2.value, t0.value-q2.to(u.second).value)
        # now comparisons
        assert t0 > q1
        assert t0 < 1.*u.yr

    def test_valid_quantity_operations2(self):
        """Check that TimeDelta is treated as a quantity where possible."""
        t0 = TimeDelta(100000., format='sec')
        f = 1./t0
        assert isinstance(f, u.Quantity)
        assert f.unit == 1./u.day
        g = 10.*u.m/u.second**2
        v = t0 * g
        assert isinstance(v, u.Quantity)
        assert v.decompose().unit == u.m / u.second
        q = np.log10(t0/u.second)
        assert isinstance(q, u.Quantity)
        assert q.value == np.log10(t0.sec)

    @pytest.mark.xfail
    def test_valid_quantity_operations3(self):
        """These should work too, but do not yet -- see #1455"""
        t0 = TimeDelta(100000., format='sec')
        s = 1.*u.m
        v = s/t0
        assert isinstance(v, u.Quantity)
        assert v.decompose().unit == u.m / u.second

    def test_invalid_quantity_operations(self):
        """Check comparisons of TimeDelta with non-time quantities fails."""
        with pytest.raises(OperandTypeError):
            TimeDelta(100000., format='sec') > 10.*u.m

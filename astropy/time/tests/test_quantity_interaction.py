# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import pytest
import numpy as np

from .. import Time, TimeDelta
from ... import units as u
from ...table import Column

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
        assert t2.value == q.value == q2.to_value(u.day)
        q3 = q-2400000.5*u.day
        t3 = Time(q3, format='mjd', scale='utc')
        assert t3.value == q3.value
        # test we can deal with two quantity arguments, with different units
        qs = 24.*36.*u.second
        t4 = Time(q3, qs, format='mjd', scale='utc')
        assert t4.value == (q3+qs).to_value(u.day)

        qy = 1990.*u.yr
        ty1 = Time(qy, format='jyear', scale='utc')
        assert ty1.value == qy.value
        ty2 = Time(qy.to(u.day), format='jyear', scale='utc')
        assert ty2.value == qy.value
        qy2 = 10.*u.yr
        tcxc = Time(qy2, format='cxcsec')
        assert tcxc.value == qy2.to_value(u.second)
        tgps = Time(qy2, format='gps')
        assert tgps.value == qy2.to_value(u.second)
        tunix = Time(qy2, format='unix')
        assert tunix.value == qy2.to_value(u.second)
        qd = 2000.*365.*u.day
        tplt = Time(qd, format='plot_date', scale='utc')
        assert tplt.value == qd.value

    def test_invalid_quantity_input(self):
        with pytest.raises(u.UnitsError):
            Time(2450000.*u.m, format='jd', scale='utc')

        with pytest.raises(u.UnitsError):
            Time(2450000.*u.dimensionless_unscaled, format='jd', scale='utc')

    def test_column_with_and_without_units(self):
        """Ensure a Column without a unit is treated as an array [#3648]"""
        a = np.arange(50000., 50010.)
        ta = Time(a, format='mjd')
        c1 = Column(np.arange(50000., 50010.), name='mjd')
        tc1 = Time(c1, format='mjd')
        assert np.all(ta == tc1)
        c2 = Column(np.arange(50000., 50010.), name='mjd', unit='day')
        tc2 = Time(c2, format='mjd')
        assert np.all(ta == tc2)
        c3 = Column(np.arange(50000., 50010.), name='mjd', unit='m')
        with pytest.raises(u.UnitsError):
            Time(c3, format='mjd')

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
        assert t1.value == t0.value+q1.to_value(u.second)
        q2 = 1.*u.day
        t2 = t0 - q2
        assert allclose_sec(t2.value, t0.value-q2.to_value(u.second))
        # check broadcasting
        q3 = np.arange(15.).reshape(3, 5) * u.hour
        t3 = t0 - q3
        assert t3.shape == q3.shape
        assert allclose_sec(t3.value, t0.value-q3.to_value(u.second))

    def test_invalid_quantity_operations(self):
        """Check that comparisons of Time with quantities does not work
        (even for time-like, since we cannot compare Time to TimeDelta)"""
        with pytest.raises(TypeError):
            Time(100000., format='cxcsec') > 10.*u.m
        with pytest.raises(TypeError):
            Time(100000., format='cxcsec') > 10.*u.second


class TestTimeDeltaQuantity():
    """Test interaction of TimeDelta with Quantities"""
    def test_valid_quantity_input(self):
        """Test that TimeDelta can take quantity input."""
        q = 500.25*u.day
        dt1 = TimeDelta(q, format='jd')
        assert dt1.value == q.value
        dt2 = TimeDelta(q, format='sec')
        assert dt2.value == q.to_value(u.second)
        dt3 = TimeDelta(q)
        assert dt3.value == q.value

    def test_invalid_quantity_input(self):
        with pytest.raises(u.UnitsError):
            TimeDelta(2450000.*u.m, format='jd')

        with pytest.raises(u.UnitsError):
            Time(2450000.*u.dimensionless_unscaled, format='jd', scale='utc')

        with pytest.raises(TypeError):
            TimeDelta(100, format='sec') > 10.*u.m

    def test_quantity_output(self):
        q = 500.25*u.day
        dt = TimeDelta(q)
        assert dt.to(u.day) == q
        assert dt.to(u.second).value == q.to_value(u.second)
        with pytest.raises(u.UnitsError):
            dt.to(u.m)

    def test_valid_quantity_operations1(self):
        """Check adding/substracting/comparing a time-valued quantity works
        with a TimeDelta.  Addition/subtraction should give TimeDelta"""
        t0 = TimeDelta(106400., format='sec')
        q1 = 10.*u.second
        t1 = t0 + q1
        assert isinstance(t1, TimeDelta)
        assert t1.value == t0.value+q1.to_value(u.second)
        q2 = 1.*u.day
        t2 = t0 - q2
        assert allclose_sec(t2.value, t0.value-q2.to_value(u.second))
        # now comparisons
        assert t0 > q1
        assert t0 < 1.*u.yr
        # and broadcasting
        q3 = np.arange(12.).reshape(4, 3) * u.hour
        t3 = t0 + q3
        assert t3.shape == q3.shape
        assert allclose_sec(t3.value, t0.value + q3.to_value(u.second))

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
        s = 1.*u.m
        v = s/t0
        assert isinstance(v, u.Quantity)
        assert v.decompose().unit == u.m / u.second
        # broadcasting
        t1 = TimeDelta(np.arange(100000., 100012.).reshape(6, 2), format='sec')
        f = np.array([1., 2.]) * u.cycle * u.Hz
        phase = f * t1
        assert isinstance(phase, u.Quantity)
        assert phase.shape == t1.shape
        assert phase.unit.is_equivalent(u.cycle)

    def test_invalid_quantity_operations(self):
        """Check comparisons of TimeDelta with non-time quantities fails."""
        with pytest.raises(TypeError):
            TimeDelta(100000., format='sec') > 10.*u.m

    def test_invalid_quantity_broadcast(self):
        """Check broadcasting rules in interactions with Quantity."""
        t0 = TimeDelta(np.arange(12.).reshape(4, 3), format='sec')
        with pytest.raises(ValueError):
            t0 + np.arange(4.) * u.s


class TestDeltaAttributes():
    def test_delta_ut1_utc(self):
        t = Time('2010-01-01 00:00:00', format='iso', scale='utc', precision=6)
        t.delta_ut1_utc = 0.3 * u.s
        assert t.ut1.iso == '2010-01-01 00:00:00.300000'
        t.delta_ut1_utc = 0.4 / 60. * u.minute
        assert t.ut1.iso == '2010-01-01 00:00:00.400000'
        with pytest.raises(u.UnitsError):
            t.delta_ut1_utc = 0.4 * u.m
        # Also check that a TimeDelta works.
        t.delta_ut1_utc = TimeDelta(0.3, format='sec')
        assert t.ut1.iso == '2010-01-01 00:00:00.300000'
        t.delta_ut1_utc = TimeDelta(0.5/24./3600., format='jd')
        assert t.ut1.iso == '2010-01-01 00:00:00.500000'

    def test_delta_tdb_tt(self):
        t = Time('2010-01-01 00:00:00', format='iso', scale='tt', precision=6)
        t.delta_tdb_tt = 20. * u.second
        assert t.tdb.iso == '2010-01-01 00:00:20.000000'
        t.delta_tdb_tt = 30. / 60. * u.minute
        assert t.tdb.iso == '2010-01-01 00:00:30.000000'
        with pytest.raises(u.UnitsError):
            t.delta_tdb_tt = 0.4 * u.m
        # Also check that a TimeDelta works.
        t.delta_tdb_tt = TimeDelta(40., format='sec')
        assert t.tdb.iso == '2010-01-01 00:00:40.000000'
        t.delta_tdb_tt = TimeDelta(50./24./3600., format='jd')
        assert t.tdb.iso == '2010-01-01 00:00:50.000000'


@pytest.mark.parametrize('q1, q2', ((5e8*u.s, None),
                                    (5e17*u.ns, None),
                                    (4e8*u.s, 1e17*u.ns),
                                    (4e14*u.us, 1e17*u.ns)))
def test_quantity_conversion_rounding(q1, q2):
    """Check that no rounding errors are incurred by unit conversion.

    This occurred before as quantities in seconds were converted to days
    before trying to split them into two-part doubles.  See gh-7622.
    """
    t = Time('2001-01-01T00:00:00.', scale='tai')
    expected = Time('2016-11-05T00:53:20.', scale='tai')
    if q2 is None:
        t0 = t + q1
    else:
        t0 = t + q1 + q2
    assert abs(t0 - expected) < 20 * u.ps
    dt1 = TimeDelta(q1, q2)
    t1 = t + dt1
    assert abs(t1 - expected) < 20 * u.ps
    dt2 = TimeDelta(q1, q2, format='sec')
    t2 = t + dt2
    assert abs(t2 - expected) < 20 * u.ps

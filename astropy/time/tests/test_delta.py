# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools
import itertools
import numpy as np
import operator

from ...tests.helper import pytest
from .. import (Time, TimeDelta, OperandTypeError, ScaleValueError,
                TIME_SCALES, TIME_DELTA_SCALES)
from ... import units as u

allclose_jd = functools.partial(np.allclose, rtol=2. ** -52, atol=0)
allclose_jd2 = functools.partial(np.allclose, rtol=2. ** -52,
                                 atol=2. ** -52)  # 20 ps atol
allclose_sec = functools.partial(np.allclose, rtol=2. ** -52,
                                 atol=2. ** -52 * 24 * 3600)  # 20 ps atol


class TestTimeDelta():
    """Test TimeDelta class"""

    def setup(self):
        self.t = Time('2010-01-01', scale='utc')
        self.t2 = Time('2010-01-02 00:00:01', scale='utc')
        self.t3 = Time('2010-01-03 01:02:03', scale='utc', precision=9,
                       in_subfmt='date_hms', out_subfmt='date_hm',
                       location=(-75.*u.degree, 30.*u.degree, 500*u.m))
        self.dt = TimeDelta(100.0, format='sec')
        self.dt_array = TimeDelta(np.arange(100, 1000, 100), format='sec')

    def test_sub(self):
        # time - time
        dt = self.t2 - self.t
        assert (repr(dt).startswith("<TimeDelta object: scale='tai' "
                                    "format='jd' value=1.00001157407"))
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
        t = Time(0.0, format='mjd', scale='utc')
        t2 = Time([0.0, 1.0], format='mjd', scale='utc')
        dt = TimeDelta(100.0, format='jd')
        dt2 = TimeDelta([100.0, 200.0], format='jd')

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
        t = Time(0.0, format='mjd', scale='utc')
        t2 = Time([0.0, 1.0], format='mjd', scale='utc')
        dt = TimeDelta(100.0, format='jd')
        dt2 = TimeDelta([100.0, 200.0], format='jd')

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

    def test_copy_timedelta(self):
        """Test copying the values of a TimeDelta object by passing it into the
        Time initializer.
        """
        t = Time(2455197.5, format='jd', scale='utc')
        t2 = Time(2455198.5, format='jd', scale='utc')
        dt = t2 - t

        dt2 = TimeDelta(dt, copy=False)
        assert dt.jd == dt2.jd
        assert dt._time.jd1 is dt2._time.jd1
        assert dt._time.jd2 is dt2._time.jd2

        dt2 = TimeDelta(dt, copy=True)
        assert dt.jd == dt2.jd
        assert dt._time.jd1 is not dt2._time.jd1
        assert dt._time.jd2 is not dt2._time.jd2

        # Include initializers
        dt2 = TimeDelta(dt, format='sec')
        assert allclose_sec(dt2.val, 86400.0)

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
            dt3 = 3. * dt
            assert allclose_jd(dt2.jd, dt3.jd)
            dt4 = dt3 / 3.
            assert allclose_jd(dt4.jd, dt.jd)
        dt5 = self.dt * np.arange(3)
        assert dt5[0].jd == 0.
        assert dt5[-1].jd == (self.dt + self.dt).jd
        with pytest.raises(OperandTypeError):
            self.dt * self.dt
        with pytest.raises(OperandTypeError):
            self.dt * self.t

    def test_keep_properties(self):
        # closes #1924 (partially)
        dt = TimeDelta(1000., format='sec')
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
        assert hasattr(t_tdb, '_delta_tdb_tt')
        assert not hasattr(t_tdb, '_delta_ut1_utc')
        t_tdb_ut1 = t_tdb.ut1
        assert hasattr(t_tdb_ut1, '_delta_tdb_tt')
        assert hasattr(t_tdb_ut1, '_delta_ut1_utc')
        t_tdb_ut1_utc = t_tdb_ut1.utc
        assert hasattr(t_tdb_ut1_utc, '_delta_tdb_tt')
        assert hasattr(t_tdb_ut1_utc, '_delta_ut1_utc')
        # adding or subtracting some time should remove the delta's
        # since these are time-dependent and should be recalculated
        for op in (operator.add, operator.sub):
            t1 = op(t_tdb, dt)
            assert not hasattr(t1, '_delta_tdb_tt')
            assert not hasattr(t1, '_delta_ut1_utc')
            t2 = op(t_tdb_ut1, dt)
            assert not hasattr(t2, '_delta_tdb_tt')
            assert not hasattr(t2, '_delta_ut1_utc')
            t3 = op(t_tdb_ut1_utc, dt)
            assert not hasattr(t3, '_delta_tdb_tt')
            assert not hasattr(t3, '_delta_ut1_utc')


class TestTimeDeltaScales():
    """Test scale conversion for Time Delta.
    Go through @taldcroft's list of expected behaviour from #1932"""

    def setup(self):
        # pick a date that includes a leap second for better testing
        self.iso_times = ['2012-06-30 12:00:00', '2012-06-30 23:59:59',
                          '2012-07-01 00:00:00', '2012-07-01 12:00:00']
        self.t = dict((scale, Time(self.iso_times, scale=scale, precision=9))
                      for scale in TIME_SCALES)
        self.dt = dict((scale, self.t[scale]-self.t[scale][0])
                       for scale in TIME_SCALES)

    def test_delta_scales_defintion(self):
        for scale in list(TIME_DELTA_SCALES) + [None]:
            TimeDelta([0., 1., 10.], format='sec', scale=scale)

        with pytest.raises(ScaleValueError):
            TimeDelta([0., 1., 10.], format='sec', scale='utc')

    @pytest.mark.parametrize(('scale1', 'scale2'),
                             list(itertools.product(TIME_SCALES, TIME_SCALES)))
    def test_scales_for_time_minus_time(self, scale1, scale2):
        """T(X) - T2(Y)  -- does T(X) - T2(Y).X and return dT(X)
        and T(X) +/- dT(Y)  -- does (in essence) (T(X).Y +/- dT(Y)).X

        I.e., time differences of two times should have the scale of the
        first time.  The one exception is UTC, which returns TAI.

        There are no timescales for which this does not work.
        """
        t1 = self.t[scale1]
        t2 = self.t[scale2]
        dt = t1 - t2
        if scale1 in TIME_DELTA_SCALES:
            assert dt.scale == scale1
        else:
            assert scale1 == 'utc'
            assert dt.scale == 'tai'

        # now check with delta time; also check reversibility
        t1_recover_t2_scale = t2 + dt
        assert t1_recover_t2_scale.scale == scale2
        t1_recover = getattr(t1_recover_t2_scale, scale1)
        assert allclose_jd(t1_recover.jd, t1.jd)
        t2_recover_t1_scale = t1 - dt
        assert t2_recover_t1_scale.scale == scale1
        t2_recover = getattr(t2_recover_t1_scale, scale2)
        assert allclose_jd(t2_recover.jd, t2.jd)

    def test_scales_for_delta_minus_delta(self):
        """dT(X) +/- dT2(Y) -- Add/substract JDs for dT(X) and dT(Y).X

        I.e. this will succeed if dT(Y) can be converted to scale X.
        Returns delta time in scale X
        """
        # geocentric timescales
        dt_tai = self.dt['tai']
        dt_tt = self.dt['tt']
        dt0 = dt_tai - dt_tt
        assert dt0.scale == 'tai'
        # tai and tt have the same scale, so differences should be the same
        assert allclose_sec(dt0.sec, 0.)

        dt_tcg = self.dt['tcg']
        dt1 = dt_tai - dt_tcg
        assert dt1.scale == 'tai'
        # tai and tcg do not have the same scale, so differences different
        assert not allclose_sec(dt1.sec, 0.)

        t_tai_tcg = self.t['tai'].tcg
        dt_tai_tcg = t_tai_tcg - t_tai_tcg[0]
        dt2 = dt_tai - dt_tai_tcg
        assert dt2.scale == 'tai'
        # but if tcg difference calculated from tai, it should roundtrip
        assert allclose_sec(dt2.sec, 0.)
        # check that if we put TCG first, we get a TCG scale back
        dt3 = dt_tai_tcg - dt_tai
        assert dt3.scale == 'tcg'
        assert allclose_sec(dt3.sec, 0.)

        for scale in 'tdb', 'tcb', 'ut1':
            with pytest.raises(TypeError):
                dt_tai - self.dt[scale]

        # barycentric timescales
        dt_tcb = self.dt['tcb']
        dt_tdb = self.dt['tdb']
        dt4 = dt_tcb - dt_tdb
        assert dt4.scale == 'tcb'
        assert not allclose_sec(dt1.sec, 0.)

        t_tcb_tdb = self.t['tcb'].tdb
        dt_tcb_tdb = t_tcb_tdb - t_tcb_tdb[0]
        dt5 = dt_tcb - dt_tcb_tdb
        assert dt5.scale == 'tcb'
        assert allclose_sec(dt5.sec, 0.)

        for scale in 'utc', 'tai', 'tt', 'tcg', 'ut1':
            with pytest.raises(TypeError):
                dt_tcb - self.dt[scale]

        # rotational timescale
        dt_ut1 = self.dt['ut1']
        dt5 = dt_ut1 - dt_ut1[-1]
        assert dt5.scale == 'ut1'
        assert dt5[-1].sec == 0.

        for scale in 'utc', 'tai', 'tt', 'tcg', 'tcb', 'tdb':
            with pytest.raises(TypeError):
                dt_ut1 - self.dt[scale]

    @pytest.mark.parametrize(
        ('scale', 'op'), list(itertools.product(TIME_SCALES,
                                                (operator.add, operator.sub))))
    def test_scales_for_delta_scale_is_none(self, scale, op):
        """T(X) +/- dT(None) or T(X) +/- Quantity(time-like)

        This is always allowed and just adds JDs, i.e., the scale of
        the TimeDelta or time-like Quantity will be taken to be X.
        The one exception is again for X=UTC, where TAI is assumed instead,
        so that a day is always defined as 86400 seconds.
        """
        dt_none = TimeDelta([0., 1., -1., 1000.], format='sec')
        assert dt_none.scale is None
        q_time = dt_none.to('s')

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

    @pytest.mark.parametrize('scale', TIME_SCALES)
    def test_delta_day_is_86400_seconds(self, scale):
        """TimeDelta or Quantity holding 1 day always means 24*60*60 seconds

        This holds true for all timescales but UTC, for which leap-second
        days are longer or shorter by one second.
        """
        t = self.t[scale]
        dt_day = TimeDelta(1., format='jd')
        q_day = dt_day.to('day')

        dt_day_leap = t[-1] - t[0]
        # ^ = exclusive or, so either equal and not UTC, or not equal and UTC
        assert allclose_jd(dt_day_leap.jd, dt_day.jd) ^ (scale == 'utc')

        t1 = t[0] + dt_day
        assert allclose_jd(t1.jd, t[-1].jd) ^ (scale == 'utc')
        t2 = q_day + t[0]
        assert allclose_jd(t2.jd, t[-1].jd) ^ (scale == 'utc')
        t3 = t[-1] - dt_day
        assert allclose_jd(t3.jd, t[0].jd) ^ (scale == 'utc')
        t4 = t[-1] - q_day
        assert allclose_jd(t4.jd, t[0].jd) ^ (scale == 'utc')

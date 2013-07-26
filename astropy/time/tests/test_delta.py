# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import numpy as np

from ...tests.helper import pytest
from .. import Time, TimeDelta, OperandTypeError

allclose_jd = functools.partial(np.allclose, rtol=1e-15, atol=0)
allclose_sec = functools.partial(np.allclose, rtol=1e-15, atol=1e-9)  # 1 nanosec atol


class TestTimeDelta():
    """Test TimeDelta class"""

    def setup(self):
        self.t = Time('2010-01-01', scale='utc')
        self.t2 = Time('2010-01-02 00:00:01', scale='utc')
        self.dt = TimeDelta(100.0, format='sec')
        self.dt_array = TimeDelta(np.arange(100,1000,100), format='sec')

    def test_sub(self):
        # time - time
        dt = self.t2 - self.t
        assert (repr(dt).startswith("<TimeDelta object: scale='tai' "
                                    "format='jd' vals=1.00001157407"))
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
        t = Time(0.0, format='mjd', scale='utc')
        t2 = Time([0.0, 1.0], format='mjd', scale='utc')
        dt = TimeDelta(100.0, format='jd')
        dt2 = TimeDelta([100.0, 200.0], format='jd')

        out = t + dt2
        assert allclose_jd(out.mjd, [100.0, 200.0])

        out = t2 + dt
        assert allclose_jd(out.mjd, [100.0, 101.0])

        out = dt + dt2
        assert allclose_jd(out.jd, [200.0, 300.0])

        # Reverse the argument order
        out = dt2 + t
        assert allclose_jd(out.mjd, [100.0, 200.0])

        out = dt + t2
        assert allclose_jd(out.mjd, [100.0, 101.0])

        out = dt2 + dt
        assert allclose_jd(out.jd, [200.0, 300.0])

    def test_sub_vector(self):
        t = Time(0.0, format='mjd', scale='utc')
        t2 = Time([0.0, 1.0], format='mjd', scale='utc')
        dt = TimeDelta(100.0, format='jd')
        dt2 = TimeDelta([100.0, 200.0], format='jd')

        out = t - dt2
        assert allclose_jd(out.mjd, [-100.0, -200.0])

        out = t2 - dt
        assert allclose_jd(out.mjd, [-100.0, -99.0])

        out = dt - dt2
        assert allclose_jd(out.jd, [0.0, -100.0])

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
            dt2 = dt+dt+dt
            dt3 = 3.*dt
            assert allclose_jd(dt2.jd, dt3.jd)
            dt4 = dt3 / 3.
            assert allclose_jd(dt4.jd, dt.jd)
        dt5 = self.dt*np.arange(3)
        assert dt5[0].jd == 0.
        assert dt5[-1].jd == (self.dt+self.dt).jd
        with pytest.raises(OperandTypeError):
            self.dt * self.dt
        with pytest.raises(OperandTypeError):
            self.dt * self.t
        with pytest.raises(TypeError):
            2. / self.dt

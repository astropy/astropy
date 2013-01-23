# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from ...tests.helper import pytest
from .. import Time, TimeDelta, OperandTypeError


class TestTimeDelta():
    """Test TimeDelta class"""

    def setup(self):
        self.t = Time('2010-01-01', scale='utc')
        self.t2 = Time('2010-01-02 00:00:01', scale='utc')
        self.dt = TimeDelta(100.0, format='sec')

    def test_sub(self):
        # time - time
        dt = self.t2 - self.t
        assert (repr(dt).startswith("<TimeDelta object: scale='tai' "
                                    "format='jd' vals=1.00001157407"))
        assert np.allclose(dt.jd, 1.00001157407)
        assert np.allclose(dt.sec, 86401.0)

        # time - delta_time
        t = self.t2 - dt
        assert t.iso == self.t.iso

        # delta_time - delta_time
        dt2 = dt - self.dt
        assert np.allclose(dt2.sec, 86301.0)

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
        assert np.allclose(dt2.sec, 86501.0)

        # delta_time + time
        dt = self.t2 - self.t
        t2 = dt + self.t
        assert t2.iso == self.t2.iso

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
        assert np.allclose(dt2.val, 86400.0)

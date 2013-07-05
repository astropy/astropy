# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from ...tests.helper import pytest
from .. import Time, TimeDelta, OperandTypeError


class TestTimeComparisons():
    """Test Comparisons of Time and TimeDelta classes"""

    def setup(self):
        self.t1 = Time(np.arange(49995,50005), format='mjd', scale='utc')
        self.t2 = Time(np.arange(49000,51000,200), format='mjd', scale='utc')

    def test_time(self):
        t1_lt_t2 = self.t1 < self.t2
        assert np.all(t1_lt_t2 == np.array([False, False, False, False, False,
                                            False,True, True, True, True]))
        t1_ge_t2 = self.t1 >= self.t2
        assert np.all(t1_ge_t2 != t1_lt_t2)

        t1_le_t2 = self.t1 <= self.t2
        assert np.all(t1_le_t2 == np.array([False, False, False, False, False,
                                            True,True, True, True, True]))
        t1_gt_t2 = self.t1 > self.t2
        assert np.all(t1_gt_t2 != t1_le_t2)

        t1_eq_t2 = self.t1 == self.t2
        assert np.all(t1_eq_t2 == np.array([False, False, False, False, False,
                                            True, False, False, False, False]))
        t1_ne_t2 = self.t1 != self.t2
        assert np.all(t1_ne_t2 != t1_eq_t2)

        t1_0_gt_t2_0 = self.t1[0] > self.t2[0]
        assert t1_0_gt_t2_0 is True
        t1_0_gt_t2 = self.t1[0] > self.t2
        assert np.all(t1_0_gt_t2 == np.array([True, True, True, True, True,
                                              False, False, False, False,
                                              False]))
        t1_gt_t2_0 = self.t1 > self.t2[0]
        assert np.all(t1_gt_t2_0 == np.array([True, True, True, True, True,
                                              True, True, True, True, True]))

    def test_timedelta(self):
        dt = self.t2 - self.t1
        with pytest.raises(OperandTypeError):
            self.t1 > dt
        dt_gt_td0 = dt > TimeDelta(0., format='sec')
        assert np.all(dt_gt_td0 == np.array([False, False, False, False, False,
                                             False, True, True, True, True]))

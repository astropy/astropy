# Test that array operations on Quantity objects work as expected

import numpy as np

from ... import units as u
from ...tests.helper import pytest


class Test1d(object):

    def setup_class(self):
        self.q1 = np.array([1, 2, 6]) * u.m
        self.q2 = np.array([4, 5, 9]) * u.s
        self.q3 = np.array([1.2, 2.2, 3.2]) * u.kg
        self.q4 = np.array([3, 4, 5]) * u.Unit(1)

    def test_add(self):
        with pytest.raises(u.UnitsException) as exc:
            self.q1 + self.q2
        assert exc.value.args[0] == ("'s' (time) and 'm' (length) are not "
                                     "convertible")

    def test_sub(self):
        with pytest.raises(u.UnitsException) as exc:
            self.q1 - self.q2
        assert exc.value.args[0] == ("'s' (time) and 'm' (length) are not "
                                     "convertible")

    def test_mult(self):
        assert np.all(self.q1 * self.q2 == np.array([4,10,54]) * u.m * u.s)

    def test_div(self):
        assert np.all(self.q2 / self.q1 == np.array([4,2.5,1.5]) * u.s / u.m)

    def test_dot(self):
        assert np.dot(self.q1, self.q2) == 68. * u.m * u.s
        assert self.q1.dot(self.q2) == 68. * u.m * u.s

    def test_min(self):
        assert self.q1.min() == 1. * u.m
        assert np.min(self.q1) == 1. * u.m
        assert self.q2.min() == 4. * u.s
        assert np.min(self.q2) == 4. * u.s

    def test_max(self):
        assert self.q1.max() == 6. * u.m
        assert np.max(self.q1) == 6. * u.m
        assert self.q2.max() == 9. * u.s
        assert np.max(self.q2) == 9. * u.s

    def test_ptp(self):
        assert self.q1.ptp() == 5. * u.m
        assert np.ptp(self.q1) == 5. * u.m
        assert self.q2.ptp() == 5. * u.s
        assert np.ptp(self.q2) == 5. * u.s

    def test_mean(self):
        assert self.q1.mean() == 3. * u.m
        assert np.mean(self.q1) == 3. * u.m
        assert self.q2.mean() == 6. * u.s
        assert np.mean(self.q2) == 6. * u.s

    def test_std(self):
        assert self.q1.std() == np.sqrt(14. / 3.) * u.m
        assert np.std(self.q1) == np.sqrt(14. / 3.) * u.m
        assert self.q2.std() == np.sqrt(14. / 3.) * u.s
        assert np.std(self.q2) == np.sqrt(14. / 3.) * u.s

    def test_var(self):
        assert self.q1.var() == (14. / 3.) * u.m
        assert np.var(self.q1) == (14. / 3.) * u.m
        assert self.q2.var() == (14. / 3.) * u.s
        assert np.var(self.q2) == (14. / 3.) * u.s

    def test_median(self):
        assert np.median(self.q1) == 2. * u.m
        assert np.median(self.q2) == 5. * u.s

    def test_round(self):
        assert np.all(np.round(self.q3) == np.array([1, 2, 3]) * u.kg)

    def test_cumsum(self):
        assert np.all(self.q1.cumsum() == np.array([1, 3, 9]) * u.m)
        assert np.all(np.cumsum(self.q1) == np.array([1, 3, 9]) * u.m)
        assert np.all(self.q2.cumsum() == np.array([4, 9, 18]) * u.s)
        assert np.all(np.cumsum(self.q2) == np.array([4, 9, 18]) * u.s)

    def test_cumprod(self):
        with pytest.raises(ValueError) as exc:
            self.q1.cumprod()
        assert exc.value.args[0] == 'cannot use cumprod on non-dimensionless Quantity arrays'
        # with pytest.raises(ValueError) as exc:
        #     np.cumprod(self.q1)
        # assert exc.value.args[0] == 'cannot use cumprod on non-dimensionless Quantity arrays'
        assert np.all(self.q4.cumprod() == np.array([3, 12, 60]) * u.Unit(1))
        # assert np.all(np.cumprod(self.q4) == np.array([3,12,60]) * u.Unit(1))


# TODO: 2d tests

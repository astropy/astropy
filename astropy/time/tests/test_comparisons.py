# Licensed under a 3-clause BSD style license - see LICENSE.rst

import operator

import pytest
import numpy as np

from astropy.time import Time, TimeDelta
import astropy.units as u


class TestTimeComparisons:
    """Test Comparisons of Time and TimeDelta classes"""

    def setup(self):
        self.t1 = Time(np.arange(49995, 50005), format='mjd', scale='utc')
        self.t2 = Time(np.arange(49000, 51000, 200), format='mjd', scale='utc')

    def test_miscompares(self):
        """
        If an incompatible object is compared to a Time object, == should
        return False and != should return True. All other comparison
        operators should raise a TypeError.
        """
        t1 = Time('J2000', scale='utc')
        for op, op_str in ((operator.ge, '>='),
                           (operator.gt, '>'),
                           (operator.le, '<='),
                           (operator.lt, '<')):
            with pytest.raises(TypeError):
                op(t1, None)
        # Keep == and != as they are specifically meant to test Time.__eq__
        # and Time.__ne__
        assert (t1 == None) is False  # noqa
        assert (t1 != None) is True  # noqa

    def test_time(self):
        t1_lt_t2 = self.t1 < self.t2
        assert np.all(t1_lt_t2 == np.array([False, False, False, False, False,
                                            False, True, True, True, True]))
        t1_ge_t2 = self.t1 >= self.t2
        assert np.all(t1_ge_t2 != t1_lt_t2)

        t1_le_t2 = self.t1 <= self.t2
        assert np.all(t1_le_t2 == np.array([False, False, False, False, False,
                                            True, True, True, True, True]))
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

    def test_time_boolean(self):
        t1_0_gt_t2_0 = self.t1[0] > self.t2[0]
        assert t1_0_gt_t2_0 is True

    def test_timedelta(self):
        dt = self.t2 - self.t1
        with pytest.raises(TypeError):
            self.t1 > dt
        dt_gt_td0 = dt > TimeDelta(0., format='sec')
        assert np.all(dt_gt_td0 == np.array([False, False, False, False, False,
                                             False, True, True, True, True]))


@pytest.mark.parametrize('swap', [True, False])
@pytest.mark.parametrize('time_delta', [True, False])
def test_isclose_time(swap, time_delta):
    """Test functionality of Time.isclose() method.

    Run every test with 2 args in original order and swapped, and using
    Quantity or TimeDelta for atol (when provided)."""

    def isclose_swap(t1, t2, **kwargs):
        if swap:
            t1, t2 = t2, t1
        if 'atol' in kwargs and time_delta:
            kwargs['atol'] = TimeDelta(kwargs['atol'])
        return t1.isclose(t2, **kwargs)

    # Start with original demonstration from #8742. In this issue both t2 == t1
    # and t3 == t1 give False, but this may change with a newer ERFA.
    t1 = Time("2018-07-24T10:41:56.807015240")
    t2 = t1 + 0.0 * u.s
    t3 = t1 + TimeDelta(0.0 * u.s)
    assert isclose_swap(t1, t2)
    assert isclose_swap(t1, t3)

    t2 = t1 + 1 * u.s
    assert isclose_swap(t1, t2, atol=1.5 / 86400 * u.day)  # Test different unit
    assert not isclose_swap(t1, t2, atol=0.5 / 86400 * u.day)

    t2 = t1 + [-1, 0, 2] * u.s
    assert np.all(isclose_swap(t1, t2, atol=1.5 * u.s) == [True, True, False])

    t2 = t1 + 3 * np.finfo(float).eps * u.day
    assert not isclose_swap(t1, t2)


def test_isclose_time_exceptions():
    t1 = Time('2020:001')
    t2 = t1 + 1 * u.s
    match = "'other' argument must support subtraction with Time"
    with pytest.raises(TypeError, match=match):
        t1.isclose(1.5)

    match = "'atol' argument must be a Quantity or TimeDelta instance, got float instead"
    with pytest.raises(TypeError, match=match):
        t1.isclose(t2, 1.5)


@pytest.mark.parametrize('swap', [True, False])
@pytest.mark.parametrize('time_delta', [True, False])
@pytest.mark.parametrize('other_quantity', [True, False])
def test_isclose_timedelta(swap, time_delta, other_quantity):
    """Test functionality of TimeDelta.isclose() method.

    Run every test with 2 args in original order and swapped, and using
    Quantity or TimeDelta for atol (when provided), and using Quantity or
    TimeDelta for the other argument."""

    def isclose_swap(t1, t2, **kwargs):
        if swap:
            t1, t2 = t2, t1
        if 'atol' in kwargs and time_delta:
            kwargs['atol'] = TimeDelta(kwargs['atol'])
        return t1.isclose(t2, **kwargs)

    def isclose_other_quantity(t1, t2, **kwargs):
        if other_quantity:
            t2 = t2.to(u.day)
        if 'atol' in kwargs and time_delta:
            kwargs['atol'] = TimeDelta(kwargs['atol'])
        return t1.isclose(t2, **kwargs)

    t1 = TimeDelta(1.0 * u.s)
    t2 = t1 + 0.0 * u.s
    t3 = t1 + TimeDelta(0.0 * u.s)
    assert isclose_swap(t1, t2)
    assert isclose_swap(t1, t3)
    assert isclose_other_quantity(t1, t2)
    assert isclose_other_quantity(t1, t3)

    t2 = t1 + 1 * u.s
    assert isclose_swap(t1, t2, atol=1.5 / 86400 * u.day)
    assert not isclose_swap(t1, t2, atol=0.5 / 86400 * u.day)
    assert isclose_other_quantity(t1, t2, atol=1.5 / 86400 * u.day)
    assert not isclose_other_quantity(t1, t2, atol=0.5 / 86400 * u.day)

    t1 = TimeDelta(0 * u.s)
    t2 = t1 + [-1, 0, 2] * u.s
    assert np.all(isclose_swap(t1, t2, atol=1.5 * u.s) == [True, True, False])
    assert np.all(isclose_other_quantity(t1, t2, atol=1.5 * u.s) == [True, True, False])

    # Check with rtol
    # 1 * 0.6 + 0.5 = 1.1 --> 1 <= 1.1 --> True
    # 0 * 0.6 + 0.5 = 0.5 --> 0 <= 0.5 --> True
    # 2 * 0.6 + 0.5 = 1.7 --> 2 <= 1.7 --> False
    assert np.all(t1.isclose(t2, atol=0.5 * u.s, rtol=0.6) == [True, True, False])

    t2 = t1 + 2 * np.finfo(float).eps * u.day
    assert not isclose_swap(t1, t2)
    assert not isclose_other_quantity(t1, t2)


def test_isclose_timedelta_exceptions():
    t1 = TimeDelta(1 * u.s)
    t2 = t1 + 1 * u.s
    match = "other' argument must support conversion to days"
    with pytest.raises(TypeError, match=match):
        t1.isclose(1.5)

    match = "'atol' argument must be a Quantity or TimeDelta instance, got float instead"
    with pytest.raises(TypeError, match=match):
        t1.isclose(t2, 1.5)

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test the propagation of info on Quantity during operations."""

import copy

import numpy as np
import pytest

from astropy import units as u


def assert_info_equal(a, b, ignore=set()):
    a_info = a.info
    b_info = b.info
    for attr in (a_info.attr_names | b_info.attr_names) - ignore:
        if attr == "unit":
            assert a_info.unit.is_equivalent(b_info.unit)
        else:
            assert getattr(a_info, attr, None) == getattr(b_info, attr, None)


def assert_no_info(a):
    __tracebackhide__ = True
    assert "info" not in a.__dict__


class TestQuantityInfo:
    @classmethod
    def setup_class(cls):
        cls.q = u.Quantity(np.arange(1.0, 5.0), "m/s")
        cls.q.info.name = "v"
        cls.q.info.description = "air speed of a african swallow"

    def test_copy(self):
        q_copy1 = self.q.copy()
        assert_info_equal(q_copy1, self.q)
        q_copy2 = copy.copy(self.q)
        assert_info_equal(q_copy2, self.q)
        q_copy3 = copy.deepcopy(self.q)
        assert_info_equal(q_copy3, self.q)

    def test_slice(self):
        q_slice = self.q[1:3]
        assert_info_equal(q_slice, self.q)
        q_take = self.q.take([0, 1])
        assert_info_equal(q_take, self.q)

    def test_item(self):
        # Scalars do not get info set (like for Column); TODO: is this OK?
        q1 = self.q[1]
        assert_no_info(q1)
        q_item = self.q.item(1)
        assert_no_info(q_item)

    def test_iter(self):
        # Scalars do not get info set.
        for q in self.q:
            assert_no_info(q)
        for q in iter(self.q):
            assert_no_info(q)

    def test_change_to_equivalent_unit(self):
        q1 = self.q.to(u.km / u.hr)
        assert_info_equal(q1, self.q)
        q2 = self.q.si
        assert_info_equal(q2, self.q)
        q3 = self.q.cgs
        assert_info_equal(q3, self.q)
        q4 = self.q.decompose()
        assert_info_equal(q4, self.q)

    def test_reshape(self):
        q = self.q.reshape(-1, 1, 2)
        assert_info_equal(q, self.q)
        q2 = q.squeeze()
        assert_info_equal(q2, self.q)

    def test_insert(self):
        q = self.q.copy()
        q.insert(1, 1 * u.cm / u.hr)
        assert_info_equal(q, self.q)

    def test_unary_op(self):
        q = -self.q
        assert_no_info(q)

    def test_binary_op(self):
        q = self.q + self.q
        assert_no_info(q)

    def test_unit_change(self):
        q = self.q * u.s
        assert_no_info(q)
        q2 = u.s / self.q
        assert_no_info(q)

    def test_inplace_unit_change(self):
        # Not sure if it is logical to keep info here!
        q = self.q.copy()
        q *= u.s
        assert_info_equal(q, self.q, ignore={"unit"})

    def test_inplace_info_name_change(self):
        # see https://github.com/astropy/astropy/issues/17449
        q = self.q.copy()

        q.info.name = "test"
        assert q.info.name == "test"

        q.info.name = None
        assert q.info.name is None

        with pytest.raises(
            TypeError, match="Expected a str value, got 2.3 with type float"
        ):
            q.info.name = 2.3


class TestStructuredQuantity:
    @classmethod
    def setup_class(cls):
        value = np.array([(1.0, 2.0), (3.0, 4.0)], dtype=[("p", "f8"), ("v", "f8")])
        cls.q = u.Quantity(value, "m, m/s")
        cls.q.info.name = "pv"
        cls.q.info.description = "Location and speed"

    def test_keying(self):
        q_p = self.q["p"]
        assert_no_info(q_p)

    def test_slicing(self):
        q = self.q[:1]
        assert_info_equal(q, self.q)

    def test_item(self):
        # Scalars do not get info set.
        q = self.q[1]
        assert_no_info(q)


class TestQuantitySubclass:
    """Regression test for gh-14514: _new_view should __array_finalize__.

    But info should be propagated only for slicing, etc.
    """

    @classmethod
    def setup_class(cls):
        class MyQuantity(u.Quantity):
            def __array_finalize__(self, obj):
                super().__array_finalize__(obj)
                if hasattr(obj, "swallow"):
                    self.swallow = obj.swallow

        cls.my_q = MyQuantity([10.0, 20.0], u.m / u.s)
        cls.my_q.swallow = "African"
        cls.my_q_w_info = cls.my_q.copy()
        cls.my_q_w_info.info.name = "swallow"

    def test_setup(self):
        assert_no_info(self.my_q)
        assert self.my_q_w_info.swallow == self.my_q.swallow
        assert self.my_q_w_info.info.name == "swallow"

    def test_slice(self):
        slc1 = self.my_q[:1]
        assert slc1.swallow == self.my_q.swallow
        assert_no_info(slc1)
        slc2 = self.my_q_w_info[1:]
        assert slc2.swallow == self.my_q.swallow
        assert_info_equal(slc2, self.my_q_w_info)

    def test_op(self):
        square1 = self.my_q**2
        assert square1.swallow == self.my_q.swallow
        assert_no_info(square1)
        square2 = self.my_q_w_info**2
        assert square2.swallow == self.my_q.swallow
        assert_no_info(square2)

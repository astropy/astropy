# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from numpy.testing import assert_array_equal

from .... import units as u
from ..core import Masked
from ....tests.helper import pytest


class TestMaskedArrayInitialization:
    def setup(self):
        self.a = np.array([1., 2.])

    def test_simple(self):
        m = np.array([True, False])
        ma = Masked(self.a, mask=m)
        assert isinstance(ma, np.ndarray)
        assert isinstance(ma, type(self.a))
        assert isinstance(ma, Masked)
        assert_array_equal(ma.data, self.a)
        assert_array_equal(ma.mask, m)


class TestMaskedQuantityInitialization(TestMaskedArrayInitialization):
    def setup(self):
        self.a = np.array([1., 2.]) * u.m


class TestFilled:
    def setup(self):
        self.a = np.arange(6.).reshape(2, 3)
        self.ma = Masked(self.a, mask=np.array([[True, False, False],
                                                [False, True, False]]))

    @pytest.mark.parametrize('fill_value', (0, 1))
    def test_filled(self, fill_value):
        expected = np.where(self.ma.mask, fill_value, self.a)
        result = self.ma.filled(fill_value)
        assert_array_equal(expected, result)


class TestMaskedArrayUfuncs:
    def setup(self):
        self.a = np.arange(6.).reshape(2, 3)
        self.b = np.array([-3., -2., -1.])
        self.c = np.array([[0.25], [0.5]])
        self.ma = Masked(self.a, mask=np.array([[True, False, False],
                                                [False, True, False]]))
        self.mb = Masked(self.b, mask=np.array([False, True, False]))
        self.mc = Masked(self.c, mask=np.array([[False], [True]]))

    def test_add(self):
        mapmb = self.ma + self.mb
        assert_array_equal(mapmb.data, self.a + self.b)
        assert_array_equal(mapmb.mask, (self.ma.mask | self.mb.mask))

    @pytest.mark.parametrize('ufunc', (np.add, np.subtract, np.divide,
                                       np.arctan2))
    def test_2op_ufunc(self, ufunc):
        ma_mb = ufunc(self.ma, self.mb)
        assert_array_equal(ma_mb.data, ufunc(self.a, self.b))
        assert_array_equal(ma_mb.mask, (self.ma.mask | self.mb.mask))

    @pytest.mark.parametrize('ufunc', (np.add, np.subtract))
    @pytest.mark.parametrize('axis', (0, 1))
    def test_reduce_filled_zero(self, ufunc, axis):
        reduction = getattr(ufunc, 'reduce')
        ma_reduce = reduction(self.ma, axis=axis)
        expected_data = reduction(self.a * (1 - self.ma.mask),
                                  axis=axis)
        expected_mask = np.logical_and.reduce(self.ma.mask, axis=axis)
        assert isinstance(ma_reduce, type(self.ma))
        assert_array_equal(ma_reduce.data, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)

    @pytest.mark.parametrize('axis', (0, 1))
    def test_sum(self, axis):
        ma_sum = self.ma.sum(axis)
        masked0 = self.ma.data * (1. - self.ma.mask)
        assert_array_equal(ma_sum.data, masked0.sum(axis))
        assert not np.any(ma_sum.mask)


class TestMaskedQuantityUfuncs(TestMaskedArrayUfuncs):
    def setup(self):
        self.a = np.arange(6.).reshape(2, 3) * u.m
        self.b = np.array([-3., -2., -1.]) * u.m
        self.c = np.array([[0.25], [0.5]]) * u.s
        self.ma = Masked(self.a, mask=np.array([[True, False, False],
                                                [False, True, False]]))
        self.mb = Masked(self.b, mask=np.array([False, True, False]))
        self.mc = Masked(self.c, mask=np.array([[False], [True]]))

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from numpy.testing import assert_array_equal

from .... import units as u
from ..core import Masked
from ....tests.helper import assert_quantity_allclose, pytest


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


class TestMaskedQuantityUfuncs(TestMaskedArrayUfuncs):
    def setup(self):
        self.a = np.arange(6.).reshape(2, 3) * u.m
        self.b = np.array([-3., -2., -1.]) * u.m
        self.c = np.array([[0.25], [0.5]]) * u.s
        self.ma = Masked(self.a, mask=np.array([[True, False, False],
                                                [False, True, False]]))
        self.mb = Masked(self.b, mask=np.array([False, True, False]))
        self.mc = Masked(self.c, mask=np.array([[False], [True]]))

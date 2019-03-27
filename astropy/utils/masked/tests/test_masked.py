# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import Longitude
from ..core import Masked
from ....tests.helper import pytest


class ArraySetup:
    def setup_arrays(self):
        self.a = np.arange(6.).reshape(2, 3)
        self.mask_a = np.array([[True, False, False],
                                [False, True, False]])
        self.b = np.array([-3., -2., -1.])
        self.mask_b = np.array([False, True, False])
        self.c = np.array([[0.25], [0.5]])
        self.mask_c = np.array([[False], [True]])

    def setup(self):
        self.setup_arrays()


class QuantitySetup(ArraySetup):
    def setup_arrays(self):
        super().setup_arrays()
        self.a = Quantity(self.a, u.m)
        self.b = Quantity(self.b, u.cm)
        self.c = Quantity(self.c, u.s)


class LongitudeSetup(ArraySetup):
    def setup_arrays(self):
        super().setup_arrays()
        self.a = Longitude(self.a, u.deg)
        self.b = Longitude(self.b, u.deg)
        self.c = Longitude(self.c, u.deg)


class TestMaskedArrayInitialization(ArraySetup):
    def setup_arrays(self):
        self.a = np.arange(6.).reshape(2, 3)
        self.mask_a = np.array([[True, False, False],
                                [False, True, False]])

    def setup(self):
        self.setup_arrays()

    def test_simple(self):
        ma = Masked(self.a, mask=self.mask_a)
        assert isinstance(ma, np.ndarray)
        assert isinstance(ma, type(self.a))
        assert isinstance(ma, Masked)
        assert_array_equal(ma.data, self.a)
        assert_array_equal(ma.mask, self.mask_a)


class TestMaskedQuantityInitialization(TestMaskedArrayInitialization):
    def setup_arrays(self):
        super().setup_arrays()
        self.a = Quantity(self.a, u.m)


# Following are tests where we trust the initializer works.


class MaskedArraySetup(ArraySetup):
    def setup(self):
        self.setup_arrays()
        self.ma = Masked(self.a, mask=self.mask_a)
        self.mb = Masked(self.b, mask=self.mask_b)
        self.mc = Masked(self.c, mask=self.mask_c)


class TestMaskedArrayFilled(MaskedArraySetup):
    @pytest.mark.parametrize('fill_value', (0, 1))
    def test_filled(self, fill_value):
        fill_value = fill_value * getattr(self.a, 'unit', 1)
        expected = self.a.copy()
        expected[self.ma.mask] = fill_value
        result = self.ma.filled(fill_value)
        assert_array_equal(expected, result)


class TestMaskedQuantityFilled(TestMaskedArrayFilled, QuantitySetup):
    pass


class TestMaskedLongitudeFilled(TestMaskedArrayFilled, LongitudeSetup):
    pass


class MaskedUfuncTests(MaskedArraySetup):
    def test_add(self):
        mapmb = self.ma + self.mb
        expected_data = self.a + self.b
        expected_mask = (self.ma.mask | self.mb.mask)
        assert_array_equal(mapmb.data, expected_data)
        assert_array_equal(mapmb.mask, expected_mask)

    @pytest.mark.parametrize('ufunc', (np.add, np.subtract, np.divide,
                                       np.arctan2))
    def test_2op_ufunc(self, ufunc):
        ma_mb = ufunc(self.ma, self.mb)
        expected_data = ufunc(self.a, self.b)
        expected_mask = (self.ma.mask | self.mb.mask)
        # Note: assert_array_equal also checks type, i.e., that, e.g.,
        # Longitude decays into an Angle.
        assert_array_equal(ma_mb.data, expected_data)
        assert_array_equal(ma_mb.mask, expected_mask)

    @pytest.mark.parametrize('ufunc', (np.add, np.subtract))
    @pytest.mark.parametrize('axis', (0, 1))
    def test_reduce_filled_zero(self, ufunc, axis):
        # axis=None for np.add tested indirectly by test_sum below.
        reduction = getattr(ufunc, 'reduce')
        ma_reduce = reduction(self.ma, axis=axis)
        expected_data = reduction(self.a * (1 - self.ma.mask),
                                  axis=axis)
        expected_mask = np.logical_and.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.data, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_sum(self, axis):
        ma_sum = self.ma.sum(axis)
        masked0 = self.ma.data * (1. - self.ma.mask)
        assert_array_equal(ma_sum.data, masked0.sum(axis))
        assert not np.any(ma_sum.mask)


class TestMaskedArrayUfuncs(MaskedUfuncTests, ArraySetup):
    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_reduce_filled_one(self, axis):
        ma_reduce = np.prod(self.ma, axis=axis)
        expected_data = np.prod(self.ma.filled(1), axis=axis)
        expected_mask = np.logical_and.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.data, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)


class TestMaskedQuantityUfuncs(MaskedUfuncTests, QuantitySetup):
    pass


class TestMaskedLongitude(MaskedUfuncTests, LongitudeSetup):
    pass

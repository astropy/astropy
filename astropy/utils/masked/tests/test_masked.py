# Licensed under a 3-clause BSD style license - see LICENSE.rst
import operator

import numpy as np
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import Longitude
from ..core import Masked
from ....tests.helper import pytest


class ArraySetup:
    _data_cls = np.ndarray

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
    _data_cls = Quantity

    def setup_arrays(self):
        super().setup_arrays()
        self.a = Quantity(self.a, u.m)
        self.b = Quantity(self.b, u.cm)
        self.c = Quantity(self.c, u.km)


class LongitudeSetup(ArraySetup):
    _data_cls = Longitude

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
        assert_array_equal(ma.unmasked, self.a)
        assert_array_equal(ma.mask, self.mask_a)
        assert ma.mask is not self.mask_a
        assert np.may_share_memory(ma.mask, self.mask_a)


class TestMaskedQuantityInitialization(TestMaskedArrayInitialization):
    def setup_arrays(self):
        super().setup_arrays()
        self.a = Quantity(self.a, u.m)


class TestMaskSetting(ArraySetup):
    def test_whole_mask_setting(self):
        ma = Masked(self.a)
        assert ma.mask.shape == ma.shape
        assert not ma.mask.any()
        ma.mask = True
        assert ma.mask.shape == ma.shape
        assert ma.mask.all()
        ma.mask = [[True], [False]]
        assert ma.mask.shape == ma.shape
        assert_array_equal(ma.mask, np.array([[True] * 3, [False] * 3]))
        ma.mask = self.mask_a
        assert ma.mask.shape == ma.shape
        assert_array_equal(ma.mask, self.mask_a)
        assert ma.mask is not self.mask_a
        assert np.may_share_memory(ma.mask, self.mask_a)

    @pytest.mark.parametrize('item', ((1, 1),
                                      slice(None, 1),
                                      (),
                                      1))
    def test_part_mask_setting(self, item):
        ma = Masked(self.a)
        ma.mask[item] = True
        expected = np.zeros(ma.shape, bool)
        expected[item] = True
        assert_array_equal(ma.mask, expected)
        ma.mask[item] = False
        assert_array_equal(ma.mask, np.zeros(ma.shape, bool))
        # Mask propagation
        mask = np.zeros(self.a.shape, bool)
        ma = Masked(self.a, mask)
        ma.mask[item] = True
        assert np.may_share_memory(ma.mask, mask)
        assert_array_equal(ma.mask, mask)


# Following are tests where we trust the initializer works.


class MaskedArraySetup(ArraySetup):
    def setup(self):
        self.setup_arrays()
        self.ma = Masked(self.a, mask=self.mask_a)
        self.mb = Masked(self.b, mask=self.mask_b)
        self.mc = Masked(self.c, mask=self.mask_c)


class TestMaskedArrayCopyFilled(MaskedArraySetup):
    def test_copy(self):
        ma_copy = self.ma.copy()
        assert type(ma_copy) is type(self.ma)
        assert_array_equal(ma_copy.unmasked, self.ma.unmasked)
        assert_array_equal(ma_copy.mask, self.ma.mask)
        assert not np.may_share_memory(ma_copy.unmasked, self.ma.unmasked)
        assert not np.may_share_memory(ma_copy.mask, self.ma.mask)

    @pytest.mark.parametrize('fill_value', (0, 1))
    def test_filled(self, fill_value):
        fill_value = fill_value * getattr(self.a, 'unit', 1)
        expected = self.a.copy()
        expected[self.ma.mask] = fill_value
        result = self.ma.unmask(fill_value)
        assert_array_equal(expected, result)

    def test_flat(self):
        ma_copy = self.ma.copy()
        ma_flat = ma_copy.flat
        # Check that single item keeps class and mask
        ma_flat1 = ma_flat[1]
        assert ma_flat1.unmasked == self.a.flat[1]
        assert ma_flat1.mask == self.mask_a.flat[1]
        # As well as getting items via iteration.
        assert all((ma.unmasked == a and ma.mask == m) for (ma, a, m)
                   in zip(self.ma.flat, self.a.flat, self.mask_a.flat))

        # check that flat works like a view of the real array
        ma_flat[1] = self.b[1]
        assert ma_flat[1] == self.b[1]
        assert ma_copy[0, 1] == self.b[1]


class TestMaskedQuantityCopyFilled(TestMaskedArrayCopyFilled, QuantitySetup):
    pass


class TestMaskedLongitudeCopyFilled(TestMaskedArrayCopyFilled, LongitudeSetup):
    pass


class TestMaskedArrayShaping(MaskedArraySetup):
    def test_reshape(self):
        ma_reshape = self.ma.reshape((6,))
        expected_data = self.a.reshape((6,))
        expected_mask = self.mask_a.reshape((6,))
        assert ma_reshape.shape == expected_data.shape
        assert_array_equal(ma_reshape.unmasked, expected_data)
        assert_array_equal(ma_reshape.mask, expected_mask)

    def test_ravel(self):
        ma_ravel = self.ma.ravel()
        expected_data = self.a.ravel()
        expected_mask = self.mask_a.ravel()
        assert ma_ravel.shape == expected_data.shape
        assert_array_equal(ma_ravel.unmasked, expected_data)
        assert_array_equal(ma_ravel.mask, expected_mask)

    def test_transpose(self):
        ma_transpose = self.ma.transpose()
        expected_data = self.a.transpose()
        expected_mask = self.mask_a.transpose()
        assert ma_transpose.shape == expected_data.shape
        assert_array_equal(ma_transpose.unmasked, expected_data)
        assert_array_equal(ma_transpose.mask, expected_mask)

    def test_iter(self):
        for ma, d, m in zip(self.ma, self.a, self.mask_a):
            assert_array_equal(ma.unmasked, d)
            assert_array_equal(ma.mask, m)


class TestMaskedArrayItems(MaskedArraySetup):
    @pytest.mark.parametrize('item', ((1, 1),
                                      slice(None, 1),
                                      (),
                                      1))
    def test_getitem(self, item):
        ma_part = self.ma[item]
        expected_data = self.a[item]
        expected_mask = self.mask_a[item]
        assert_array_equal(ma_part.unmasked, expected_data)
        assert_array_equal(ma_part.mask, expected_mask)

    @pytest.mark.parametrize('item', ((1, 1),
                                      slice(None, 1),
                                      (),
                                      1))
    def test_setitem(self, item):
        base = self.ma.copy()
        expected_data = self.a.copy()
        expected_mask = self.mask_a.copy()
        for mask in True, False:
            value = Masked(base[0, 0], mask)
            base[item] = value
            expected_data[item] = value.unmasked
            expected_mask[item] = value.mask
            assert_array_equal(base.unmasked, expected_data)
            assert_array_equal(base.mask, expected_mask)


class TestMaskedQuantityItems(TestMaskedArrayItems, QuantitySetup):
    pass


class TestMaskedLongitudeItems(TestMaskedArrayItems, LongitudeSetup):
    pass


class TestMaskedArrayOperators(MaskedArraySetup):
    @pytest.mark.parametrize('op', (operator.add, operator.sub))
    def test_add_subtract(self, op):
        mapmb = op(self.ma, self.mb)
        expected_data = op(self.a, self.b)
        expected_mask = (self.ma.mask | self.mb.mask)
        # Note: assert_array_equal also checks type, i.e., that, e.g.,
        # Longitude decays into an Angle.
        assert_array_equal(mapmb.unmasked, expected_data)
        assert_array_equal(mapmb.mask, expected_mask)


class TestMaskedQuantityOperators(TestMaskedArrayOperators, QuantitySetup):
    pass


class TestMaskedLongitudeOperators(TestMaskedArrayOperators, LongitudeSetup):
    pass


class MaskedUfuncTests(MaskedArraySetup):
    @pytest.mark.parametrize('ufunc', (np.add, np.subtract, np.divide,
                                       np.arctan2, np.minimum))
    def test_2op_ufunc(self, ufunc):
        ma_mb = ufunc(self.ma, self.mb)
        expected_data = ufunc(self.a, self.b)
        expected_mask = (self.ma.mask | self.mb.mask)
        # Note: assert_array_equal also checks type, i.e., that, e.g.,
        # Longitude decays into an Angle.
        assert_array_equal(ma_mb.unmasked, expected_data)
        assert_array_equal(ma_mb.mask, expected_mask)

    def test_3op_ufunc(self):
        ma_mb = np.clip(self.ma, self.mb, self.mc)
        expected_data = np.clip(self.a, self.b, self.c)
        expected_mask = (self.ma.mask | self.mb.mask | self.mc.mask)
        assert_array_equal(ma_mb.unmasked, expected_data)
        assert_array_equal(ma_mb.mask, expected_mask)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_add_reduce(self, axis):
        ma_reduce = np.add.reduce(self.ma, axis=axis)
        expected_data = np.add.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_minimum_reduce(self, axis):
        ma_reduce = np.minimum.reduce(self.ma, axis=axis)
        expected_data = np.minimum.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_maximum_reduce(self, axis):
        ma_reduce = np.maximum.reduce(self.ma, axis=axis)
        expected_data = np.maximum.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)


class TestMaskedArrayUfuncs(MaskedUfuncTests, ArraySetup):
    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_multiply_reduce(self, axis):
        ma_reduce = np.multiply.reduce(self.ma, axis=axis)
        expected_data = np.multiply.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)


class TestMaskedQuantityUfuncs(MaskedUfuncTests, QuantitySetup):
    pass


class TestMaskedLongitudeUfuncs(MaskedUfuncTests, LongitudeSetup):
    pass


class TestMaskedArrayMethods(MaskedArraySetup):
    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_sum(self, axis):
        ma_sum = self.ma.sum(axis)
        expected_data = self.a.sum(axis)
        expected_mask = self.ma.mask.any(axis)
        assert_array_equal(ma_sum.unmasked, expected_data)
        assert_array_equal(ma_sum.mask, expected_mask)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_mean(self, axis):
        ma_mean = self.ma.mean(axis)
        filled = self.a.copy()
        filled[self.mask_a] = 0.
        count = 1 - self.ma.mask.astype(int)
        expected_data = filled.sum(axis) / count.sum(axis)
        expected_mask = self.ma.mask.all(axis)
        assert_array_equal(ma_mean.unmasked, expected_data)
        assert_array_equal(ma_mean.mask, expected_mask)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_var(self, axis):
        ma_var = self.ma.var(axis)
        filled = (self.a - self.ma.mean(axis, keepdims=True))**2
        filled[self.mask_a] = 0.
        count = (1 - self.ma.mask.astype(int)).sum(axis)
        expected_data = filled.sum(axis) / count
        expected_mask = self.ma.mask.all(axis)
        assert_array_equal(ma_var.unmasked, expected_data)
        assert_array_equal(ma_var.mask, expected_mask)
        ma_var1 = self.ma.var(axis, ddof=1)
        expected_data1 = filled.sum(axis) / (count - 1)
        expected_mask1 = self.ma.mask.all(axis) | (count <= 1)
        assert_array_equal(ma_var1.unmasked, expected_data1)
        assert_array_equal(ma_var1.mask, expected_mask1)
        ma_var5 = self.ma.var(axis, ddof=5)
        assert np.all(~np.isfinite(ma_var5.unmasked))
        assert ma_var5.mask.all()

    def test_std(self):
        ma_std = self.ma.std(1, ddof=1)
        ma_var1 = self.ma.var(1, ddof=1)
        expected = np.sqrt(ma_var1)
        assert_array_equal(ma_std.unmasked, expected.unmasked)
        assert_array_equal(ma_std.mask, expected.mask)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_min(self, axis):
        ma_min = self.ma.min(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.max()
        expected_data = filled.min(axis)
        assert_array_equal(ma_min.unmasked, expected_data)
        assert not np.any(ma_min.mask)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_max(self, axis):
        ma_max = self.ma.max(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.min()
        expected_data = filled.max(axis)
        assert_array_equal(ma_max.unmasked, expected_data)
        assert not np.any(ma_max.mask)

    def test_all_explicit(self):
        a1 = np.array([[1., 2.],
                       [3., 4.]])
        a2 = np.array([[1., 0.],
                       [3., 4.]])
        if self._data_cls is not np.ndarray:
            a1 = self._data_cls(a1, self.a.unit)
            a2 = self._data_cls(a2, self.a.unit)
        ma1 = Masked(a1, mask=[[False, False],
                               [True, True]])
        ma2 = Masked(a2, mask=[[False, True],
                               [False, True]])
        ma1_eq_ma2 = ma1 == ma2
        assert_array_equal(ma1_eq_ma2.unmasked, np.array([[True, False],
                                                          [True, True]]))
        assert_array_equal(ma1_eq_ma2.mask, np.array([[False, True],
                                                      [True, True]]))
        assert ma1_eq_ma2.all()
        assert not (ma1 != ma2).all()
        ma_eq1 = ma1_eq_ma2.all(1)
        assert_array_equal(ma_eq1.mask, np.array([False, True]))
        assert bool(ma_eq1[0]) is True
        assert bool(ma_eq1[1]) is False
        ma_eq0 = ma1_eq_ma2.all(0)
        assert_array_equal(ma_eq0.mask, np.array([False, True]))
        assert bool(ma_eq1[0]) is True
        assert bool(ma_eq1[1]) is False

    @pytest.mark.parametrize('method', ['any', 'all'])
    @pytest.mark.parametrize('array,axis', [
        ('a', 0), ('a', 1), ('a', None),
        ('b', None),
        ('c', 0), ('c', 1), ('c', None)])
    def test_all_and_any(self, array, axis, method):
        ma = getattr(self, 'm'+array)
        ma_eq = ma == ma
        a_eq = ma.unmasked == ma.unmasked
        ma_all_or_any = getattr(ma_eq, method)(axis=axis)
        filled = ma_eq.unmasked.copy()
        filled[ma_eq.mask] = method == 'all'
        a_all_or_any = getattr(filled, method)(axis=axis)
        all_masked = ma.mask.all(axis)
        assert_array_equal(ma_all_or_any.mask, all_masked)
        assert_array_equal(ma_all_or_any.unmasked, a_all_or_any)
        # interpretation as bool
        as_bool = [bool(a) for a in ma_all_or_any.ravel()]
        expected = [bool(a) for a in (a_all_or_any & ~all_masked).ravel()]
        assert as_bool == expected


class TestMaskedQuantityMethods(TestMaskedArrayMethods, QuantitySetup):
    pass


class TestMaskedLongitudeMethods(TestMaskedArrayMethods, LongitudeSetup):
    pass


class TestMaskedArrayRepr(MaskedArraySetup):
    def test_array_str(self):
        # very blunt check they work at all.
        str(self.ma)
        str(self.mb)
        str(self.mb)

    def test_scalar_str(self):
        assert self.mb[0].shape == ()
        str(self.mb[0])

    def test_array_repr(self):
        repr(self.ma)
        repr(self.mb)
        repr(self.mc)

    def test_scalar_repr(self):
        repr(self.mb[0])


class TestMaskedQuantityRepr(TestMaskedArrayRepr, QuantitySetup):
    pass


class TestMaskedArrayConcatenation(MaskedArraySetup):
    def test_concatenate(self):
        mb = self.mb[np.newaxis]
        concat_a_b = np.concatenate((self.ma, mb), axis=0)
        expected_data = np.concatenate((self.a, self.b[np.newaxis]), axis=0)
        expected_mask = np.concatenate((self.mask_a, self.mask_b[np.newaxis]),
                                       axis=0)
        assert_array_equal(concat_a_b.unmasked, expected_data)
        assert_array_equal(concat_a_b.mask, expected_mask)


class TestMaskedQuantityConcatenation(TestMaskedArrayConcatenation,
                                      QuantitySetup):
    pass


class TestMaskedLongitudeConcatenation(TestMaskedArrayConcatenation,
                                       LongitudeSetup):
    pass

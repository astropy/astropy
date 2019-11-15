# Licensed under a 3-clause BSD style license - see LICENSE.rst
import operator

import numpy as np
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import Longitude
from ..core import Masked
from ....tests.helper import pytest


def assert_masked_equal(a, b):
    assert_array_equal(a.unmasked, b.unmasked)
    assert_array_equal(a.mask, b.mask)


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


class MaskedItemTests(MaskedArraySetup):
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

    @pytest.mark.parametrize('indices,axis', [
        ([0, 1], 1), ([0, 1], 0), ([0, 1], None), ([[0, 1], [2, 3]], None)])
    def test_take(self, indices, axis):
        ma_take = self.ma.take(indices, axis=axis)
        expected_data = self.a.take(indices, axis=axis)
        expected_mask = self.mask_a.take(indices, axis=axis)
        assert_array_equal(ma_take.unmasked, expected_data)
        assert_array_equal(ma_take.mask, expected_mask)
        ma_take2 = np.take(self.ma, indices, axis=axis)
        assert_masked_equal(ma_take2, ma_take)

    @pytest.mark.parametrize('item', ((1, 1),
                                      slice(None, 1),
                                      (),
                                      1))
    def test_setitem(self, item):
        base = self.ma.copy()
        expected_data = self.a.copy()
        expected_mask = self.mask_a.copy()
        for mask in True, False:
            value = Masked(self.a[0, 0], mask)
            base[item] = value
            expected_data[item] = value.unmasked
            expected_mask[item] = value.mask
            assert_array_equal(base.unmasked, expected_data)
            assert_array_equal(base.mask, expected_mask)

    @pytest.mark.parametrize('item', ((1, 1),
                                      slice(None, 1),
                                      (),
                                      1))
    def test_setitem_without_mask(self, item):
        base = self.ma.copy()
        expected_data = self.a.copy()
        expected_mask = self.mask_a.copy()
        value = self.a[0, 0]
        base[item] = value
        expected_data[item] = value
        expected_mask[item] = False
        assert_array_equal(base.unmasked, expected_data)
        assert_array_equal(base.mask, expected_mask)


class TestMaskedArrayItems(MaskedItemTests):
    # TODO: make this work for Quantity & Longitude as well.
    # Or decide that it just shouldn't work...
    @pytest.mark.parametrize('item', ((1, 1),
                                      slice(None, 1),
                                      (),
                                      1))
    def test_setitem_np_ma_masked(self, item):
        base = self.ma.copy()
        expected_mask = self.mask_a.copy()
        base[item] = np.ma.masked
        expected_mask[item] = True
        assert_array_equal(base.unmasked, self.a)
        assert_array_equal(base.mask, expected_mask)


class TestMaskedQuantityItems(MaskedItemTests, QuantitySetup):
    pass


class TestMaskedLongitudeItems(MaskedItemTests, LongitudeSetup):
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

    @pytest.mark.parametrize('op', (operator.eq, operator.ne))
    def test_equality(self, op):
        mapmb = op(self.ma, self.mb)
        expected_data = op(self.a, self.b)
        expected_mask = (self.ma.mask | self.mb.mask)
        # Note: assert_array_equal also checks type, i.e., that boolean
        # output is represented as plain Masked ndarray.
        assert_array_equal(mapmb.unmasked, expected_data)
        assert_array_equal(mapmb.mask, expected_mask)

    def test_matmul(self):
        result = self.ma.T @ self.ma
        assert_array_equal(result.unmasked, self.a.T @ self.a)
        mask1 = np.any(self.mask_a, axis=0)
        expected_mask = np.logical_or.outer(mask1, mask1)
        assert_array_equal(result.mask, expected_mask)
        result2 = self.ma.T @ self.a
        assert_array_equal(result2.unmasked, self.a.T @ self.a)
        expected_mask2 = np.logical_or.outer(mask1, np.zeros(3, bool))
        assert_array_equal(result2.mask, expected_mask2)
        result3 = self.a.T @ self.ma
        assert_array_equal(result3.unmasked, self.a.T @ self.a)
        expected_mask3 = np.logical_or.outer(np.zeros(3, bool), mask1)
        assert_array_equal(result3.mask, expected_mask3)

    def test_matvec(self):
        result = self.ma @ self.mb
        assert np.all(result.mask)
        assert_array_equal(result.unmasked, self.a @ self.b)
        # Just using the masked vector still has all elements masked.
        result2 = self.a @ self.mb
        assert np.all(result2.mask)
        assert_array_equal(result2.unmasked, self.a @ self.b)
        new_ma = self.ma.copy()
        new_ma.mask[0, 0] = False
        result3 = new_ma @ self.b
        assert_array_equal(result3.unmasked, self.a @ self.b)
        assert_array_equal(result3.mask, new_ma.mask.any(-1))

    def test_vecmat(self):
        result = self.mb @ self.ma.T
        assert np.all(result.mask)
        assert_array_equal(result.unmasked, self.b @ self.a.T)
        result2 = self.b @ self.ma.T
        assert np.all(result2.mask)
        assert_array_equal(result2.unmasked, self.b @ self.a.T)
        new_ma = self.ma.T.copy()
        new_ma.mask[0, 0] = False
        result3 = self.b @ new_ma
        assert_array_equal(result3.unmasked, self.b @ self.a.T)
        assert_array_equal(result3.mask, new_ma.mask.any(0))

    def test_vecvec(self):
        result = self.mb @ self.mb
        assert result.shape == ()
        assert result.mask
        assert result.unmasked == self.b @ self.b
        mb_no_mask = Masked(self.b, False)
        result2 = mb_no_mask @ mb_no_mask
        assert not result2.mask


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

    @pytest.mark.parametrize('ufunc', (np.add, np.subtract, np.divide,
                                       np.arctan2, np.minimum))
    def test_ufunc_inplace(self, ufunc):
        ma_mb = ufunc(self.ma, self.mb)
        out = Masked(np.zeros_like(ma_mb.unmasked))
        result = ufunc(self.ma, self.mb, out=out)
        assert result is out
        assert_masked_equal(result, ma_mb)

    def test_ufunc_inplace_error(self):
        out = np.zeros(self.ma.shape)
        with pytest.raises(TypeError):
            np.add(self.ma, self.mb, out=out)

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
    # multiply.reduce does not work with units, so test only for plain array.
    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_multiply_reduce(self, axis):
        ma_reduce = np.multiply.reduce(self.ma, axis=axis)
        expected_data = np.multiply.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)


class TestMaskedQuantityUfuncs(MaskedUfuncTests, QuantitySetup):
    def test_ufunc_inplace_error2(self):
        out = Masked(np.zeros(self.ma.shape))
        with pytest.raises(TypeError):
            np.add(self.ma, self.mb, out=out)


class TestMaskedLongitudeUfuncs(MaskedUfuncTests, LongitudeSetup):
    def test_ufunc_inplace_quantity_initial(self):
        out = Masked(np.zeros(self.ma.shape) << u.m)
        result = np.add(self.ma, self.mb, out=out)
        assert result is out
        expected = np.add(self.ma, self.mb).view(Quantity)
        assert_masked_equal(result, expected)


class TestMaskedArrayMethods(MaskedArraySetup):
    def test_round(self):
        # Goes via ufunc, hence easy.
        mrc = self.mc.round()
        expected = Masked(self.c.round(), self.mask_c)
        assert_masked_equal(mrc, expected)

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

    def test_mean_float16(self):
        ma = self.ma.astype('f2')
        ma_mean = ma.mean()
        expected = self.ma.mean().astype('f2')
        assert_masked_equal(ma_mean, expected)

    def test_mean_inplace(self):
        expected = self.ma.mean(1)
        out = Masked(np.zeros_like(expected.unmasked))
        result = self.ma.mean(1, out=out)
        assert result is out
        assert_masked_equal(out, expected)

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
        assert_masked_equal(ma_std, expected)

    def test_std_inplace(self):
        expected = self.ma.std(1, ddof=1)
        out = Masked(np.zeros_like(expected.unmasked))
        result = self.ma.std(1, ddof=1, out=out)
        assert result is out
        assert_masked_equal(result, expected)

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

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_argmin(self, axis):
        ma_argmin = self.ma.argmin(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.max()
        expected_data = filled.argmin(axis)
        assert_array_equal(ma_argmin, expected_data)

    @pytest.mark.parametrize('axis', (0, 1, None))
    def test_argmax(self, axis):
        ma_argmax = self.ma.argmax(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.min()
        expected_data = filled.argmax(axis)
        assert_array_equal(ma_argmax, expected_data)

    @pytest.mark.parametrize('axis', (0, 1))
    def test_argsort(self, axis):
        ma_argsort = self.ma.argsort(axis)
        filled = self.a.copy()
        filled[self.mask_a] = self.a.max() * 1.1
        expected_data = filled.argsort(axis)
        assert_array_equal(ma_argsort, expected_data)

    @pytest.mark.parametrize('axis', (0, 1))
    def test_sort(self, axis):
        ma_sort = self.ma.copy()
        ma_sort.sort(axis)
        indices = self.ma.argsort(axis)
        expected_data = np.take_along_axis(self.ma.unmasked, indices, axis)
        expected_mask = np.take_along_axis(self.ma.mask, indices, axis)
        assert_array_equal(ma_sort.unmasked, expected_data)
        assert_array_equal(ma_sort.mask, expected_mask)

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

    def test_any_inplace(self):
        ma_eq = self.ma == self.ma
        expected = ma_eq.any(1)
        out = Masked(np.zeros_like(expected.unmasked))
        result = ma_eq.any(1, out=out)
        assert result is out
        assert_masked_equal(result, expected)

    @pytest.mark.parametrize('offset', (0, 1))
    def test_diagonal(self, offset):
        mda = self.ma.diagonal(offset=offset)
        expected = Masked(self.a.diagonal(offset=offset),
                          self.mask_a.diagonal(offset=offset))
        assert_masked_equal(mda, expected)

    @pytest.mark.parametrize('offset', (0, 1))
    def test_trace(self, offset):
        mta = self.ma.trace(offset=offset)
        expected = Masked(self.a.trace(offset=offset),
                          self.mask_a.trace(offset=offset, dtype=bool))
        assert_masked_equal(mta, expected)


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

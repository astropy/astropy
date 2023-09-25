# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test all functions covered by __array_function__.

Here, run through all functions, with simple tests just to check the helpers.
More complicated tests of functionality, including with subclasses, are done
in test_functions.

TODO: finish full coverage (see also `~astropy.utils.masked.function_helpers`)
- np.linalg
- np.fft (is there any point?)
- np.lib.nanfunctions

"""
import inspect
import itertools

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy.units.tests.test_quantity_non_ufuncs import get_wrapped_functions
from astropy.utils.compat import NUMPY_LT_1_23, NUMPY_LT_1_24, NUMPY_LT_1_25
from astropy.utils.masked import Masked, MaskedNDArray
from astropy.utils.masked.function_helpers import (
    APPLY_TO_BOTH_FUNCTIONS,
    DISPATCHED_FUNCTIONS,
    IGNORED_FUNCTIONS,
    MASKED_SAFE_FUNCTIONS,
    UNSUPPORTED_FUNCTIONS,
)

from .test_masked import MaskedArraySetup, assert_masked_equal

all_wrapped_functions = get_wrapped_functions(np)
all_wrapped = set(all_wrapped_functions.values())


class BasicTestSetup(MaskedArraySetup):
    def check(self, func, *args, **kwargs):
        out = func(self.ma, *args, **kwargs)
        expected = Masked(
            func(self.a, *args, **kwargs), mask=func(self.mask_a, *args, **kwargs)
        )
        assert_masked_equal(out, expected)

    def check2(self, func, *args, **kwargs):
        out = func(self.ma, self.mb, *args, **kwargs)
        expected = Masked(
            func(self.a, self.b, *args, **kwargs),
            mask=func(self.mask_a, self.mask_b, *args, **kwargs),
        )
        if isinstance(out, (tuple, list)):
            for o, x in zip(out, expected):
                assert_masked_equal(o, x)
        else:
            assert_masked_equal(out, expected)


class NoMaskTestSetup(MaskedArraySetup):
    def check(self, func, *args, **kwargs):
        o = func(self.ma, *args, **kwargs)
        expected = func(self.a, *args, **kwargs)
        assert_array_equal(o, expected)


class InvariantMaskTestSetup(MaskedArraySetup):
    def check(self, func, *args, **kwargs):
        o = func(self.ma, *args, **kwargs)
        expected = func(self.a, *args, **kwargs)
        assert_array_equal(o.unmasked, expected)
        assert_array_equal(o.mask, self.mask_a)


class TestShapeInformation(BasicTestSetup):
    def test_shape(self):
        assert np.shape(self.ma) == (2, 3)

    def test_size(self):
        assert np.size(self.ma) == 6

    def test_ndim(self):
        assert np.ndim(self.ma) == 2


class TestShapeManipulation(BasicTestSetup):
    # Note: do not parametrize the below, since test names are used
    # to check coverage.
    def test_reshape(self):
        self.check(np.reshape, (6, 1))

    def test_ravel(self):
        self.check(np.ravel)

    def test_moveaxis(self):
        self.check(np.moveaxis, 0, 1)

    def test_rollaxis(self):
        self.check(np.rollaxis, 0, 2)

    def test_swapaxes(self):
        self.check(np.swapaxes, 0, 1)

    def test_transpose(self):
        self.check(np.transpose)

    def test_atleast_1d(self):
        self.check(np.atleast_1d)
        o, so = np.atleast_1d(self.mb[0], self.mc[0])
        assert o.shape == o.mask.shape == so.shape == so.mask.shape == (1,)

    def test_atleast_2d(self):
        self.check(np.atleast_2d)
        o, so = np.atleast_2d(self.mb[0], self.mc[0])
        assert o.shape == o.mask.shape == so.shape == so.mask.shape == (1, 1)

    def test_atleast_3d(self):
        self.check(np.atleast_3d)
        o, so = np.atleast_3d(self.mb[0], self.mc[0])
        assert o.shape == o.mask.shape == so.shape == so.mask.shape == (1, 1, 1)

    def test_expand_dims(self):
        self.check(np.expand_dims, 1)

    def test_squeeze(self):
        o = np.squeeze(self.mc)
        assert o.shape == o.mask.shape == (2,)
        assert_array_equal(o.unmasked, self.c.squeeze())
        assert_array_equal(o.mask, self.mask_c.squeeze())

    def test_flip(self):
        self.check(np.flip)

    def test_fliplr(self):
        self.check(np.fliplr)

    def test_flipud(self):
        self.check(np.flipud)

    def test_rot90(self):
        self.check(np.rot90)

    def test_broadcast_to(self):
        self.check(np.broadcast_to, (3, 2, 3))
        self.check(np.broadcast_to, (3, 2, 3), subok=False)

    def test_broadcast_arrays(self):
        self.check2(np.broadcast_arrays)
        self.check2(np.broadcast_arrays, subok=False)


class TestArgFunctions(MaskedArraySetup):
    def check(self, function, *args, fill_value=np.nan, **kwargs):
        o = function(self.ma, *args, **kwargs)
        a_filled = self.ma.filled(fill_value=fill_value)
        expected = function(a_filled, *args, **kwargs)
        assert_array_equal(o, expected)

    def test_argmin(self):
        self.check(np.argmin, fill_value=np.inf)

    def test_argmax(self):
        self.check(np.argmax, fill_value=-np.inf)

    def test_argsort(self):
        self.check(np.argsort, fill_value=np.nan)

    def test_lexsort(self):
        self.check(np.lexsort, fill_value=np.nan)

    def test_nonzero(self):
        self.check(np.nonzero, fill_value=0.0)

    @pytest.mark.filterwarnings("ignore:Calling nonzero on 0d arrays is deprecated")
    def test_nonzero_0d(self):
        res1 = Masked(1, mask=False).nonzero()
        assert len(res1) == 1
        assert_array_equal(res1[0], np.ones(()).nonzero()[0])
        res2 = Masked(1, mask=True).nonzero()
        assert len(res2) == 1
        assert_array_equal(res2[0], np.zeros(()).nonzero()[0])

    def test_argwhere(self):
        self.check(np.argwhere, fill_value=0.0)

    def test_argpartition(self):
        self.check(np.argpartition, 2, fill_value=np.inf)

    def test_flatnonzero(self):
        self.check(np.flatnonzero, fill_value=0.0)


class TestAlongAxis(MaskedArraySetup):
    def test_take_along_axis(self):
        indices = np.expand_dims(np.argmax(self.ma, axis=0), axis=0)
        out = np.take_along_axis(self.ma, indices, axis=0)
        expected = np.take_along_axis(self.a, indices, axis=0)
        expected_mask = np.take_along_axis(self.mask_a, indices, axis=0)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_put_along_axis(self):
        ma = self.ma.copy()
        indices = np.expand_dims(np.argmax(self.ma, axis=0), axis=0)
        np.put_along_axis(ma, indices, axis=0, values=-1)
        expected = self.a.copy()
        np.put_along_axis(expected, indices, axis=0, values=-1)
        assert_array_equal(ma.unmasked, expected)
        assert_array_equal(ma.mask, self.mask_a)
        np.put_along_axis(ma, indices, axis=0, values=np.ma.masked)
        assert_array_equal(ma.unmasked, expected)
        expected_mask = self.mask_a.copy()
        np.put_along_axis(expected_mask, indices, axis=0, values=True)
        assert_array_equal(ma.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1))
    def test_apply_along_axis(self, axis):
        out = np.apply_along_axis(np.square, axis, self.ma)
        expected = np.apply_along_axis(np.square, axis, self.a)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, self.mask_a)

    @pytest.mark.parametrize("axes", [(1,), 0, (0, -1)])
    def test_apply_over_axes(self, axes):
        def function(x, axis):
            return np.mean(np.square(x), axis)

        out = np.apply_over_axes(function, self.ma, axes)
        expected = self.ma
        for axis in axes if isinstance(axes, tuple) else (axes,):
            expected = (expected**2).mean(axis, keepdims=True)
        assert_array_equal(out.unmasked, expected.unmasked)
        assert_array_equal(out.mask, expected.mask)

    def test_apply_over_axes_no_reduction(self):
        out = np.apply_over_axes(np.cumsum, self.ma, 0)
        expected = self.ma.cumsum(axis=0)
        assert_masked_equal(out, expected)

    def test_apply_over_axes_wrong_size(self):
        with pytest.raises(ValueError, match="not.*correct shape"):
            np.apply_over_axes(lambda x, axis: x[..., np.newaxis], self.ma, 0)


class TestIndicesFrom(NoMaskTestSetup):
    @classmethod
    def setup_class(self):
        self.a = np.arange(9).reshape(3, 3)
        self.mask_a = np.eye(3, dtype=bool)
        self.ma = Masked(self.a, self.mask_a)

    def test_diag_indices_from(self):
        self.check(np.diag_indices_from)

    def test_triu_indices_from(self):
        self.check(np.triu_indices_from)

    def test_tril_indices_from(self):
        self.check(np.tril_indices_from)


class TestRealImag(InvariantMaskTestSetup):
    @classmethod
    def setup_class(self):
        self.a = np.array([1 + 2j, 3 + 4j])
        self.mask_a = np.array([True, False])
        self.ma = Masked(self.a, mask=self.mask_a)

    def test_real(self):
        self.check(np.real)

    def test_imag(self):
        self.check(np.imag)


class TestCopyAndCreation(InvariantMaskTestSetup):
    def test_copy(self):
        self.check(np.copy)
        # Also as kwarg
        copy = np.copy(a=self.ma)
        assert_array_equal(copy, self.ma)

    def test_asfarray(self):
        self.check(np.asfarray)
        farray = np.asfarray(a=self.ma)
        assert_array_equal(farray, self.ma)


class TestArrayCreation(MaskedArraySetup):
    def test_empty_like(self):
        o = np.empty_like(self.ma)
        assert o.shape == (2, 3)
        assert isinstance(o, Masked)
        assert isinstance(o, np.ndarray)
        o2 = np.empty_like(prototype=self.ma)
        assert o2.shape == (2, 3)
        assert isinstance(o2, Masked)
        assert isinstance(o2, np.ndarray)
        o3 = np.empty_like(self.ma, subok=False)
        assert type(o3) is MaskedNDArray

    def test_zeros_like(self):
        o = np.zeros_like(self.ma)
        assert_array_equal(o.unmasked, np.zeros_like(self.a))
        assert_array_equal(o.mask, np.zeros_like(self.mask_a))
        o2 = np.zeros_like(a=self.ma)
        assert_array_equal(o2.unmasked, np.zeros_like(self.a))
        assert_array_equal(o2.mask, np.zeros_like(self.mask_a))

    def test_ones_like(self):
        o = np.ones_like(self.ma)
        assert_array_equal(o.unmasked, np.ones_like(self.a))
        assert_array_equal(o.mask, np.zeros_like(self.mask_a))
        o2 = np.ones_like(a=self.ma)
        assert_array_equal(o2.unmasked, np.ones_like(self.a))
        assert_array_equal(o2.mask, np.zeros_like(self.mask_a))

    @pytest.mark.parametrize("value", [0.5, Masked(0.5, mask=True), np.ma.masked])
    def test_full_like(self, value):
        o = np.full_like(self.ma, value)
        if value is np.ma.masked:
            expected = Masked(o.unmasked, True)
        else:
            expected = Masked(np.empty_like(self.a))
            expected[...] = value
        assert_array_equal(o.unmasked, expected.unmasked)
        assert_array_equal(o.mask, expected.mask)


class TestAccessingParts(BasicTestSetup):
    def test_diag(self):
        self.check(np.diag)

    def test_diag_1d_input(self):
        ma = self.ma.ravel()
        o = np.diag(ma)
        assert_array_equal(o.unmasked, np.diag(self.a.ravel()))
        assert_array_equal(o.mask, np.diag(self.mask_a.ravel()))

    def test_diagonal(self):
        self.check(np.diagonal)

    def test_diagflat(self):
        self.check(np.diagflat)

    def test_compress(self):
        o = np.compress([True, False], self.ma, axis=0)
        expected = np.compress([True, False], self.a, axis=0)
        expected_mask = np.compress([True, False], self.mask_a, axis=0)
        assert_array_equal(o.unmasked, expected)
        assert_array_equal(o.mask, expected_mask)

    def test_extract(self):
        o = np.extract([True, False, True], self.ma)
        expected = np.extract([True, False, True], self.a)
        expected_mask = np.extract([True, False, True], self.mask_a)
        assert_array_equal(o.unmasked, expected)
        assert_array_equal(o.mask, expected_mask)

    def test_delete(self):
        self.check(np.delete, slice(1, 2), 0)
        self.check(np.delete, [0, 2], 1)

    def test_roll(self):
        self.check(np.roll, 1)
        self.check(np.roll, 1, axis=0)

    def test_take(self):
        self.check(np.take, [0, 1], axis=1)
        self.check(np.take, 1)


class TestSettingParts(MaskedArraySetup):
    def test_put(self):
        ma = self.ma.copy()
        v = Masked([50, 150], [False, True])
        np.put(ma, [0, 2], v)
        expected = self.a.copy()
        np.put(expected, [0, 2], [50, 150])
        expected_mask = self.mask_a.copy()
        np.put(expected_mask, [0, 2], [False, True])
        assert_array_equal(ma.unmasked, expected)
        assert_array_equal(ma.mask, expected_mask)

        with pytest.raises(TypeError):
            # Indices cannot be masked.
            np.put(ma, Masked([0, 2]), v)

        with pytest.raises(TypeError):
            # Array to put masked values in must be masked.
            np.put(self.a.copy(), [0, 2], v)

    def test_putmask(self):
        ma = self.ma.flatten()
        mask = [True, False, False, False, True, False]
        values = Masked(
            np.arange(100, 650, 100), mask=[False, True, True, True, False, False]
        )
        np.putmask(ma, mask, values)
        expected = self.a.flatten()
        np.putmask(expected, mask, values.unmasked)
        expected_mask = self.mask_a.flatten()
        np.putmask(expected_mask, mask, values.mask)
        assert_array_equal(ma.unmasked, expected)
        assert_array_equal(ma.mask, expected_mask)

        with pytest.raises(TypeError):
            np.putmask(self.a.flatten(), mask, values)

    def test_place(self):
        ma = self.ma.flatten()
        mask = [True, False, False, False, True, False]
        values = Masked([100, 200], mask=[False, True])
        np.place(ma, mask, values)
        expected = self.a.flatten()
        np.place(expected, mask, values.unmasked)
        expected_mask = self.mask_a.flatten()
        np.place(expected_mask, mask, values.mask)
        assert_array_equal(ma.unmasked, expected)
        assert_array_equal(ma.mask, expected_mask)

        with pytest.raises(TypeError):
            np.place(self.a.flatten(), mask, values)

    def test_copyto(self):
        ma = self.ma.flatten()
        mask = [True, False, False, False, True, False]
        values = Masked(
            np.arange(100, 650, 100), mask=[False, True, True, True, False, False]
        )
        np.copyto(ma, values, where=mask)
        expected = self.a.flatten()
        np.copyto(expected, values.unmasked, where=mask)
        expected_mask = self.mask_a.flatten()
        np.copyto(expected_mask, values.mask, where=mask)
        assert_array_equal(ma.unmasked, expected)
        assert_array_equal(ma.mask, expected_mask)

        with pytest.raises(TypeError):
            np.copyto(self.a.flatten(), values, where=mask)

    @pytest.mark.parametrize("value", [0.25, np.ma.masked])
    def test_fill_diagonal(self, value):
        ma = self.ma[:2, :2].copy()
        np.fill_diagonal(ma, value)
        expected = ma.copy()
        expected[np.diag_indices_from(expected)] = value
        assert_array_equal(ma.unmasked, expected.unmasked)
        assert_array_equal(ma.mask, expected.mask)


class TestRepeat(BasicTestSetup):
    def test_tile(self):
        self.check(np.tile, 2)

    def test_repeat(self):
        self.check(np.repeat, 2)

    def test_resize(self):
        self.check(np.resize, (4, 4))


class TestConcatenate(MaskedArraySetup):
    # More tests at TestMaskedArrayConcatenation in test_functions.
    def check(self, func, *args, **kwargs):
        ma_list = kwargs.pop("ma_list", [self.ma, self.ma])
        a_list = [Masked(ma).unmasked for ma in ma_list]
        m_list = [Masked(ma).mask for ma in ma_list]
        o = func(ma_list, *args, **kwargs)
        expected = func(a_list, *args, **kwargs)
        expected_mask = func(m_list, *args, **kwargs)
        assert_array_equal(o.unmasked, expected)
        assert_array_equal(o.mask, expected_mask)

    def test_concatenate(self):
        self.check(np.concatenate)
        self.check(np.concatenate, axis=1)
        self.check(np.concatenate, ma_list=[self.a, self.ma])
        self.check(np.concatenate, dtype="f4")

        out = Masked(np.empty((4, 3)))
        result = np.concatenate([self.ma, self.ma], out=out)
        assert out is result
        expected = np.concatenate([self.a, self.a])
        expected_mask = np.concatenate([self.mask_a, self.mask_a])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

        with pytest.raises(TypeError):
            np.concatenate([self.ma, self.ma], out=np.empty((4, 3)))

    def test_stack(self):
        self.check(np.stack)

    def test_column_stack(self):
        self.check(np.column_stack)

    def test_hstack(self):
        self.check(np.hstack)

    def test_vstack(self):
        self.check(np.vstack)

    def test_dstack(self):
        self.check(np.dstack)

    def test_block(self):
        self.check(np.block)

        out = np.block([[0.0, Masked(1.0, True)], [Masked(1, False), Masked(2, False)]])
        expected = np.array([[0, 1.0], [1, 2]])
        expected_mask = np.array([[False, True], [False, False]])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_append(self):
        out = np.append(self.ma, self.mc, axis=1)
        expected = np.append(self.a, self.c, axis=1)
        expected_mask = np.append(self.mask_a, self.mask_c, axis=1)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_insert(self):
        obj = (1, 1)
        values = Masked([50.0, 25.0], mask=[True, False])
        out = np.insert(self.ma.flatten(), obj, values)
        expected = np.insert(self.a.flatten(), obj, [50.0, 25.0])
        expected_mask = np.insert(self.mask_a.flatten(), obj, [True, False])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

        with pytest.raises(TypeError):
            np.insert(self.a.flatten(), obj, values)

        with pytest.raises(TypeError):
            np.insert(self.ma.flatten(), Masked(obj), values)


class TestSplit:
    @classmethod
    def setup_class(self):
        self.a = np.arange(54.0).reshape(3, 3, 6)
        self.mask_a = np.zeros(self.a.shape, dtype=bool)
        self.mask_a[1, 1, 1] = True
        self.mask_a[0, 1, 4] = True
        self.mask_a[1, 2, 5] = True
        self.ma = Masked(self.a, mask=self.mask_a)

    def check(self, func, *args, **kwargs):
        out = func(self.ma, *args, **kwargs)
        expected = func(self.a, *args, **kwargs)
        expected_mask = func(self.mask_a, *args, **kwargs)
        assert len(out) == len(expected)
        for o, x, xm in zip(out, expected, expected_mask):
            assert_array_equal(o.unmasked, x)
            assert_array_equal(o.mask, xm)

    def test_split(self):
        self.check(np.split, [1])

    def test_array_split(self):
        self.check(np.array_split, 2)

    def test_hsplit(self):
        self.check(np.hsplit, [1, 4])

    def test_vsplit(self):
        self.check(np.vsplit, [1])

    def test_dsplit(self):
        self.check(np.dsplit, [1])


class TestMethodLikes(MaskedArraySetup):
    def check(self, function, *args, method=None, **kwargs):
        if method is None:
            method = function.__name__

        o = function(self.ma, *args, **kwargs)
        x = getattr(self.ma, method)(*args, **kwargs)
        assert_masked_equal(o, x)

    def test_max(self):
        self.check(np.max, method="max")

    def test_min(self):
        self.check(np.min, method="min")

    def test_amax(self):
        self.check(np.amax, method="max")

    def test_amin(self):
        self.check(np.amin, method="min")

    def test_sum(self):
        self.check(np.sum)

    def test_cumsum(self):
        self.check(np.cumsum)

    def test_any(self):
        self.check(np.any)

    def test_all(self):
        self.check(np.all)

    # NUMPY_LT_1_25
    @pytest.mark.filterwarnings("ignore:`sometrue` is deprecated as of NumPy 1.25.0")
    def test_sometrue(self):
        self.check(np.sometrue, method="any")

    # NUMPY_LT_1_25
    @pytest.mark.filterwarnings("ignore:`alltrue` is deprecated as of NumPy 1.25.0")
    def test_alltrue(self):
        self.check(np.alltrue, method="all")

    def test_prod(self):
        self.check(np.prod)

    # NUMPY_LT_1_25
    @pytest.mark.filterwarnings("ignore:`product` is deprecated as of NumPy 1.25.0")
    def test_product(self):
        self.check(np.product, method="prod")

    def test_cumprod(self):
        self.check(np.cumprod)

    # NUMPY_LT_1_25
    @pytest.mark.filterwarnings("ignore:`cumproduct` is deprecated as of NumPy 1.25.0")
    def test_cumproduct(self):
        self.check(np.cumproduct, method="cumprod")

    def test_round(self):
        self.check(np.round, method="round")

    # NUMPY_LT_1_25
    @pytest.mark.filterwarnings("ignore:`round_` is deprecated as of NumPy 1.25.0")
    def test_round_(self):
        self.check(np.round_, method="round")

    def test_around(self):
        self.check(np.around, method="round")

    def test_clip(self):
        self.check(np.clip, 2.0, 4.0)
        self.check(np.clip, self.mb, self.mc)

    def test_mean(self):
        self.check(np.mean)

    def test_std(self):
        self.check(np.std)

    def test_var(self):
        self.check(np.var)


class TestUfuncLike(InvariantMaskTestSetup):
    def test_fix(self):
        self.check(np.fix)
        # Check np.fix with out argument for completeness
        # (Note: could be done in self.check, but np.fix is the only
        # invariant mask function that has `out`, so no point.)
        out = np.zeros_like(self.ma)
        result = np.fix(self.ma, out=out)
        assert result is out
        expected = np.fix(self.a)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, self.mask_a)

    def test_angle(self):
        a = np.array([1 + 0j, 0 + 1j, 1 + 1j, 0 + 0j])
        mask_a = np.array([True, False, True, False])
        ma = Masked(a, mask=mask_a)
        out = np.angle(ma)
        expected = np.angle(ma.unmasked)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, mask_a)

    def test_i0(self):
        self.check(np.i0)

    def test_sinc(self):
        self.check(np.sinc)

    def test_where(self):
        mask = [True, False, True]
        out = np.where(mask, self.ma, 1000.0)
        expected = np.where(mask, self.a, 1000.0)
        expected_mask = np.where(mask, self.mask_a, False)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

        mask2 = Masked(mask, [True, False, False])
        out2 = np.where(mask2, self.ma, 1000.0)
        expected2 = np.where(mask, self.a, 1000.0)
        expected_mask2 = np.where(mask, self.mask_a, False) | mask2.mask
        assert_array_equal(out2.unmasked, expected2)
        assert_array_equal(out2.mask, expected_mask2)

    def test_where_single_arg(self):
        m = Masked(np.arange(3), mask=[True, False, False])
        out = np.where(m)
        expected = m.nonzero()
        assert isinstance(out, tuple) and len(out) == 1
        assert_array_equal(out[0], expected[0])

    def test_where_wrong_number_of_arg(self):
        with pytest.raises(ValueError, match="either both or neither"):
            np.where([True, False, False], self.a)

    def test_choose(self):
        a = np.array([0, 1]).reshape((2, 1))
        result = np.choose(a, (self.ma, self.mb))
        expected = np.choose(a, (self.a, self.b))
        expected_mask = np.choose(a, (self.mask_a, self.mask_b))
        assert_array_equal(result.unmasked, expected)
        assert_array_equal(result.mask, expected_mask)

        out = np.zeros_like(result)
        result2 = np.choose(a, (self.ma, self.mb), out=out)
        assert result2 is out
        assert_array_equal(result2, result)

        with pytest.raises(TypeError):
            np.choose(a, (self.ma, self.mb), out=np.zeros_like(expected))

    def test_choose_masked(self):
        ma = Masked(np.array([-1, 1]), mask=[True, False]).reshape((2, 1))
        out = ma.choose((self.ma, self.mb))
        expected = np.choose(ma.filled(0), (self.a, self.b))
        expected_mask = np.choose(ma.filled(0), (self.mask_a, self.mask_b)) | ma.mask
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

        with pytest.raises(ValueError):
            ma.unmasked.choose((self.ma, self.mb))

    @pytest.mark.parametrize("default", [-1.0, np.ma.masked, Masked(-1, mask=True)])
    def test_select(self, default):
        a, mask_a, ma = self.a, self.mask_a, self.ma
        out = np.select([a < 1.5, a > 3.5], [ma, ma + 1], default=default)
        expected = np.select(
            [a < 1.5, a > 3.5],
            [a, a + 1],
            default=-1 if default is not np.ma.masked else 0,
        )
        expected_mask = np.select(
            [a < 1.5, a > 3.5],
            [mask_a, mask_a],
            default=getattr(default, "mask", False),
        )
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_real_if_close(self):
        a = np.array([1 + 0j, 0 + 1j, 1 + 1j, 0 + 0j])
        mask_a = np.array([True, False, True, False])
        ma = Masked(a, mask=mask_a)
        out = np.real_if_close(ma)
        expected = np.real_if_close(a)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, mask_a)

    def test_tril(self):
        self.check(np.tril)

    def test_triu(self):
        self.check(np.triu)

    def test_unwrap(self):
        self.check(np.unwrap)

    def test_nan_to_num(self):
        self.check(np.nan_to_num)
        ma = Masked([np.nan, 1.0], mask=[True, False])
        o = np.nan_to_num(ma, copy=False)
        assert_masked_equal(o, Masked([0.0, 1.0], mask=[True, False]))
        assert ma is o


class TestUfuncLikeTests:
    @classmethod
    def setup_class(self):
        self.a = np.array([[-np.inf, +np.inf, np.nan, 3.0, 4.0]] * 2)
        self.mask_a = np.array([[False] * 5, [True] * 4 + [False]])
        self.ma = Masked(self.a, mask=self.mask_a)
        self.b = np.array([[3.0001], [3.9999]])
        self.mask_b = np.array([[True], [False]])
        self.mb = Masked(self.b, mask=self.mask_b)

    def check(self, func):
        out = func(self.ma)
        expected = func(self.a)
        assert type(out) is MaskedNDArray
        assert out.dtype.kind == "b"
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, self.mask_a)
        assert not np.may_share_memory(out.mask, self.mask_a)

    def test_isposinf(self):
        self.check(np.isposinf)

    def test_isneginf(self):
        self.check(np.isneginf)

    def test_isreal(self):
        self.check(np.isreal)
        o = np.isreal(Masked([1.0 + 1j], mask=False))
        assert not o.unmasked and not o.mask
        o = np.isreal(Masked([1.0 + 1j], mask=True))
        assert not o.unmasked and o.mask

    def test_iscomplex(self):
        self.check(np.iscomplex)
        o = np.iscomplex(Masked([1.0 + 1j], mask=False))
        assert o.unmasked and not o.mask
        o = np.iscomplex(Masked([1.0 + 1j], mask=True))
        assert o.unmasked and o.mask

    def test_isclose(self):
        out = np.isclose(self.ma, self.mb, atol=0.01)
        expected = np.isclose(self.ma, self.mb, atol=0.01)
        expected_mask = self.mask_a | self.mask_b
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_allclose(self):
        out = np.allclose(self.ma, self.mb, atol=0.01)
        expected = np.isclose(self.ma, self.mb, atol=0.01)[
            self.mask_a | self.mask_b
        ].all()
        assert_array_equal(out, expected)

    def test_array_equal(self):
        assert not np.array_equal(self.ma, self.ma)
        assert not np.array_equal(self.ma, self.a)
        assert np.array_equal(self.ma, self.ma, equal_nan=True)
        assert np.array_equal(self.ma, self.a, equal_nan=True)
        assert not np.array_equal(self.ma, self.mb)
        ma2 = self.ma.copy()
        ma2.mask |= np.isnan(self.a)
        assert np.array_equal(ma2, self.ma)

    def test_array_equiv(self):
        assert np.array_equiv(self.mb, self.mb)
        assert np.array_equiv(self.mb, self.b)
        assert not np.array_equiv(self.ma, self.mb)
        assert np.array_equiv(self.mb, np.stack([self.mb, self.mb]))


class TestOuterLikeFunctions(MaskedArraySetup):
    def test_outer(self):
        result = np.outer(self.ma, self.mb)
        expected_data = np.outer(self.a.ravel(), self.b.ravel())
        expected_mask = np.logical_or.outer(self.mask_a.ravel(), self.mask_b.ravel())
        assert_array_equal(result.unmasked, expected_data)
        assert_array_equal(result.mask, expected_mask)

        out = np.zeros_like(result)
        result2 = np.outer(self.ma, self.mb, out=out)
        assert result2 is out
        assert result2 is not result
        assert_masked_equal(result2, result)

        out2 = np.zeros_like(result.unmasked)
        with pytest.raises(TypeError):
            np.outer(self.ma, self.mb, out=out2)

    def test_kron(self):
        result = np.kron(self.ma, self.mb)
        expected_data = np.kron(self.a, self.b)
        expected_mask = np.logical_or.outer(self.mask_a, self.mask_b).reshape(
            result.shape
        )
        assert_array_equal(result.unmasked, expected_data)
        assert_array_equal(result.mask, expected_mask)


class TestReductionLikeFunctions(MaskedArraySetup):
    def test_average(self):
        o = np.average(self.ma)
        assert_masked_equal(o, self.ma.mean())

        o = np.average(self.ma, weights=self.mb, axis=-1)
        expected = np.average(self.a, weights=self.b, axis=-1)
        expected_mask = (self.mask_a | self.mask_b).any(-1)
        assert_array_equal(o.unmasked, expected)
        assert_array_equal(o.mask, expected_mask)

    @pytest.mark.parametrize("kwargs", [{}, {"axis": 0}])
    def test_ptp(self, kwargs):
        o = np.ptp(self.ma, **kwargs)
        expected = self.ma.max(**kwargs) - self.ma.min(**kwargs)
        assert_array_equal(o.unmasked, expected.unmasked)
        assert_array_equal(o.mask, expected.mask)
        out = np.zeros_like(expected)
        o2 = np.ptp(self.ma, out=out, **kwargs)
        assert o2 is out
        assert_array_equal(o2.unmasked, expected.unmasked)
        assert_array_equal(o2.mask, expected.mask)
        # Check method does the same thing.
        o3 = self.ma.ptp(**kwargs)
        assert_array_equal(o3.unmasked, expected.unmasked)
        assert_array_equal(o3.mask, expected.mask)

    def test_trace(self):
        o = np.trace(self.ma)
        expected = np.trace(self.a)
        expected_mask = np.trace(self.mask_a).astype(bool)
        assert_array_equal(o.unmasked, expected)
        assert_array_equal(o.mask, expected_mask)

    @pytest.mark.parametrize("axis", [0, 1, None])
    def test_count_nonzero(self, axis):
        o = np.count_nonzero(self.ma, axis=axis)
        expected = np.count_nonzero(self.ma.filled(0), axis=axis)
        assert_array_equal(o, expected)


@pytest.mark.filterwarnings("ignore:all-nan")
class TestPartitionLikeFunctions:
    @classmethod
    def setup_class(self):
        self.a = np.arange(36.0).reshape(6, 6)
        self.mask_a = np.zeros_like(self.a, bool)
        # On purpose fill diagonal, so we get all masked elements.
        self.mask_a[np.tril_indices_from(self.a)] = True
        self.ma = Masked(self.a, mask=self.mask_a)

    def check(self, function, *args, **kwargs):
        # Check function by comparing to nan-equivalent, with masked
        # values set to NaN.
        o = function(self.ma, *args, **kwargs)
        nanfunc = getattr(np, "nan" + function.__name__)
        nanfilled = self.ma.filled(np.nan)
        expected = nanfunc(nanfilled, *args, **kwargs)
        assert_array_equal(o.filled(np.nan), expected)
        assert_array_equal(o.mask, np.isnan(expected))
        # Also check that we can give an output MaskedArray.
        if NUMPY_LT_1_25 and kwargs.get("keepdims", False):
            # numpy bug gh-22714 prevents using out with keepdims=True.
            # This is fixed in numpy 1.25.
            return

        out = np.zeros_like(o)
        o2 = function(self.ma, *args, out=out, **kwargs)
        assert o2 is out
        assert_masked_equal(o2, o)
        # But that a regular array cannot be used since it has no mask.
        with pytest.raises(TypeError):
            function(self.ma, *args, out=np.zeros_like(expected), **kwargs)

    @pytest.mark.parametrize("keepdims", [False, True])
    @pytest.mark.parametrize("axis", [None, 0, 1])
    def test_median(self, axis, keepdims):
        self.check(np.median, axis=axis, keepdims=keepdims)

    @pytest.mark.parametrize("keepdims", [False, True])
    @pytest.mark.parametrize("axis", [None, 0, 1])
    def test_quantile(self, axis, keepdims):
        self.check(np.quantile, q=[0.25, 0.5], axis=axis, keepdims=keepdims)

    def test_quantile_out_of_range(self):
        with pytest.raises(ValueError, match="must be in the range"):
            np.quantile(self.ma, q=1.5)

    @pytest.mark.parametrize("axis", [None, 0, 1])
    def test_percentile(self, axis):
        self.check(np.percentile, q=50, axis=axis)


class TestIntDiffFunctions(MaskedArraySetup):
    def test_diff(self):
        out = np.diff(self.ma)
        expected = np.diff(self.a)
        expected_mask = self.mask_a[:, 1:] | self.mask_a[:, :-1]
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_diff_prepend_append(self):
        out = np.diff(self.ma, prepend=Masked(-1, mask=True), append=1)
        expected = np.diff(self.a, prepend=-1, append=1.0)
        mask = np.concatenate(
            [np.ones((2, 1), bool), self.mask_a, np.zeros((2, 1), bool)], axis=-1
        )
        expected_mask = mask[:, 1:] | mask[:, :-1]
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_trapz(self):
        ma = self.ma.copy()
        ma.mask[1] = False
        out = np.trapz(ma)
        assert_array_equal(out.unmasked, np.trapz(self.a))
        assert_array_equal(out.mask, np.array([True, False]))

    def test_gradient(self):
        out = np.gradient(self.ma)
        expected = np.gradient(self.a)
        expected_mask = [
            (self.mask_a[1:] | self.mask_a[:-1]).repeat(2, axis=0),
            np.stack(
                [
                    self.mask_a[:, 0] | self.mask_a[:, 1],
                    self.mask_a[:, 0] | self.mask_a[:, 2],
                    self.mask_a[:, 1] | self.mask_a[:, 2],
                ],
                axis=-1,
            ),
        ]

        for o, x, m in zip(out, expected, expected_mask):
            assert_array_equal(o.unmasked, x)
            assert_array_equal(o.mask, m)


class TestSpaceFunctions:
    @classmethod
    def setup_class(self):
        self.a = np.arange(1.0, 7.0).reshape(2, 3)
        self.mask_a = np.array(
            [
                [True, False, False],
                [False, True, False],
            ]
        )
        self.ma = Masked(self.a, mask=self.mask_a)
        self.b = np.array([2.5, 10.0, 3.0])
        self.mask_b = np.array([False, True, False])
        self.mb = Masked(self.b, mask=self.mask_b)

    def check(self, function, *args, **kwargs):
        out = function(self.ma, self.mb, 5)
        expected = function(self.a, self.b, 5)
        expected_mask = np.broadcast_to(
            self.mask_a | self.mask_b, expected.shape
        ).copy()
        # TODO: make implementation that also ensures start point mask is
        # determined just by start point? (as for geomspace in numpy 1.20)?
        expected_mask[-1] = self.mask_b
        if function is np.geomspace:
            expected_mask[0] = self.mask_a

        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_linspace(self):
        self.check(np.linspace, 5)

    def test_logspace(self):
        self.check(np.logspace, 10)

    def test_geomspace(self):
        self.check(np.geomspace, 5)


class TestInterpolationFunctions(MaskedArraySetup):
    def test_interp(self):
        xp = np.arange(5.0)
        fp = np.array([1.0, 5.0, 6.0, 19.0, 20.0])
        mask_fp = np.array([False, False, False, True, False])
        mfp = Masked(fp, mask=mask_fp)
        x = np.array([1.5, 17.0])
        mask_x = np.array([False, True])
        mx = Masked(x, mask=mask_x)
        out = np.interp(mx, xp, mfp)
        expected = np.interp(x, xp[~mask_fp], fp[~mask_fp])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, mask_x)

    def test_piecewise(self):
        condlist = [self.a < 1, self.a >= 1]
        out = np.piecewise(self.ma, condlist, [Masked(-1, mask=True), 1.0])
        expected = np.piecewise(self.a, condlist, [-1, 1.0])
        expected_mask = np.piecewise(self.mask_a, condlist, [True, False])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)
        condlist2 = [self.a < 1, self.a >= 3]
        out2 = np.piecewise(
            self.ma,
            condlist2,
            [Masked(-1, True), 1, lambda x: Masked(np.full_like(x, 2.0), mask=~x.mask)],
        )
        expected = np.piecewise(self.a, condlist2, [-1, 1, 2])
        expected_mask = np.piecewise(
            self.mask_a, condlist2, [True, False, lambda x: ~x]
        )
        assert_array_equal(out2.unmasked, expected)
        assert_array_equal(out2.mask, expected_mask)

        with pytest.raises(ValueError, match="with 2 condition"):
            np.piecewise(self.ma, condlist2, [])

    def test_regression_12978(self):
        """Regression tests for https://github.com/astropy/astropy/pull/12978"""
        # This case produced incorrect results
        mask = [False, True, False]
        x = np.array([1, 2, 3])
        xp = Masked(np.array([1, 2, 3]), mask=mask)
        fp = Masked(np.array([1, 2, 3]), mask=mask)
        result = np.interp(x, xp, fp)
        assert_array_equal(result, x)

        # This case raised a ValueError
        xp = np.array([1, 3])
        fp = Masked(np.array([1, 3]))
        result = np.interp(x, xp, fp)
        assert_array_equal(result, x)


class TestBincount(MaskedArraySetup):
    def test_bincount(self):
        i = np.array([1, 1, 2, 3, 2, 4])
        mask_i = np.array([True, False, False, True, False, False])
        mi = Masked(i, mask=mask_i)
        out = np.bincount(mi)
        expected = np.bincount(i[~mask_i])
        assert_array_equal(out, expected)
        w = np.arange(len(i))
        mask_w = np.array([True] + [False] * 5)
        mw = Masked(w, mask=mask_w)
        out2 = np.bincount(i, mw)
        expected = np.bincount(i, w)
        expected_mask = np.array([False, True, False, False, False])
        assert_array_equal(out2.unmasked, expected)
        assert_array_equal(out2.mask, expected_mask)

        out3 = np.bincount(mi, mw)
        expected = np.bincount(i[~mask_i], w[~mask_i])
        expected_mask = np.array([False, False, False, False, False])
        assert_array_equal(out3.unmasked, expected)
        assert_array_equal(out3.mask, expected_mask)


class TestSortFunctions(MaskedArraySetup):
    def test_sort(self):
        o = np.sort(self.ma)
        expected = self.ma.copy()
        expected.sort()
        assert_masked_equal(o, expected)

    def test_sort_complex(self):
        ma = Masked(
            np.array([1 + 2j, 0 + 4j, 3 + 0j, -1 - 1j]),
            mask=[True, False, False, False],
        )
        o = np.sort_complex(ma)
        indx = np.lexsort((ma.unmasked.imag, ma.unmasked.real, ma.mask))
        expected = ma[indx]
        assert_masked_equal(o, expected)

    @pytest.mark.skipif(not NUMPY_LT_1_24, reason="np.msort is deprecated")
    def test_msort(self):
        o = np.msort(self.ma)
        expected = np.sort(self.ma, axis=0)
        assert_masked_equal(o, expected)

    def test_partition(self):
        o = np.partition(self.ma, 1)
        expected = self.ma.copy()
        expected.partition(1)
        assert_masked_equal(o, expected)


class TestStringFunctions:
    # More elaborate tests done in test_masked.py
    @classmethod
    def setup_class(self):
        self.ma = Masked(np.arange(3), mask=[True, False, False])

    def test_array2string(self):
        out0 = np.array2string(self.ma)
        assert out0 == "[— 1 2]"
        # Arguments are interpreted as usual.
        out1 = np.array2string(self.ma, separator=", ")
        assert out1 == "[—, 1, 2]"
        # If we do pass in a formatter, though, it should be used.
        out2 = np.array2string(self.ma, separator=", ", formatter={"all": hex})
        assert out2 == "[———, 0x1, 0x2]"
        # Also as positional argument (no, nobody will do this!)
        out3 = np.array2string(
            self.ma, None, None, None, ", ", "", np._NoValue, {"int": hex}
        )
        assert out3 == out2
        # But not if the formatter is not relevant for us.
        out4 = np.array2string(self.ma, separator=", ", formatter={"float": hex})
        assert out4 == out1

    def test_array_repr(self):
        out = np.array_repr(self.ma)
        assert out == "MaskedNDArray([—, 1, 2])"
        ma2 = self.ma.astype("f4")
        out2 = np.array_repr(ma2)
        assert out2 == "MaskedNDArray([——, 1., 2.], dtype=float32)"

    def test_array_str(self):
        out = np.array_str(self.ma)
        assert out == "[— 1 2]"


class TestBitFunctions:
    @classmethod
    def setup_class(self):
        self.a = np.array([15, 255, 0], dtype="u1")
        self.mask_a = np.array([False, True, False])
        self.ma = Masked(self.a, mask=self.mask_a)
        self.b = np.unpackbits(self.a).reshape(6, 4)
        self.mask_b = np.array([False] * 15 + [True, True] + [False] * 7).reshape(6, 4)
        self.mb = Masked(self.b, mask=self.mask_b)

    @pytest.mark.parametrize("axis", [None, 1, 0])
    def test_packbits(self, axis):
        out = np.packbits(self.mb, axis=axis)
        if axis is None:
            expected = self.a
        else:
            expected = np.packbits(self.b, axis=axis)
        expected_mask = np.packbits(self.mask_b, axis=axis) > 0
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_unpackbits(self):
        out = np.unpackbits(self.ma)
        mask = np.where(self.mask_a, np.uint8(255), np.uint8(0))
        expected_mask = np.unpackbits(mask) > 0
        assert_array_equal(out.unmasked, self.b.ravel())
        assert_array_equal(out.mask, expected_mask)


class TestIndexFunctions(MaskedArraySetup):
    """Does not seem much sense to support these..."""

    def test_unravel_index(self):
        with pytest.raises(TypeError):
            np.unravel_index(self.ma, 3)

    def test_ravel_multi_index(self):
        with pytest.raises(TypeError):
            np.ravel_multi_index((self.ma,), 3)

    def test_ix_(self):
        with pytest.raises(TypeError):
            np.ix_(self.ma)


class TestDtypeFunctions(MaskedArraySetup):
    def check(self, function, *args, **kwargs):
        out = function(self.ma, *args, **kwargs)
        expected = function(self.a, *args, **kwargs)
        assert out == expected

    def test_common_type(self):
        self.check(np.common_type)

    def test_result_type(self):
        self.check(np.result_type)

    def test_can_cast(self):
        self.check(np.can_cast, self.a.dtype)
        self.check(np.can_cast, "f4")

    def test_min_scalar_type(self):
        out = np.min_scalar_type(self.ma[0, 0])
        expected = np.min_scalar_type(self.a[0, 0])
        assert out == expected

    def test_iscomplexobj(self):
        self.check(np.iscomplexobj)

    def test_isrealobj(self):
        self.check(np.isrealobj)


class TestMeshGrid(MaskedArraySetup):
    def test_meshgrid(self):
        a = np.arange(1.0, 4.0)
        mask_a = np.array([True, False, False])
        ma = Masked(a, mask=mask_a)
        b = np.array([2.5, 10.0, 3.0, 4.0])
        mask_b = np.array([False, True, False, True])
        mb = Masked(b, mask=mask_b)
        oa, ob = np.meshgrid(ma, mb)
        xa, xb = np.broadcast_arrays(a, b[:, np.newaxis])
        ma, mb = np.broadcast_arrays(mask_a, mask_b[:, np.newaxis])
        for o, x, m in ((oa, xa, ma), (ob, xb, mb)):
            assert_array_equal(o.unmasked, x)
            assert_array_equal(o.mask, m)


class TestMemoryFunctions(MaskedArraySetup):
    def test_shares_memory(self):
        assert np.shares_memory(self.ma, self.ma.unmasked)
        assert not np.shares_memory(self.ma, self.ma.mask)

    def test_may_share_memory(self):
        assert np.may_share_memory(self.ma, self.ma.unmasked)
        assert not np.may_share_memory(self.ma, self.ma.mask)


class TestDatetimeFunctions:
    # Could in principle support np.is_busday, np.busday_count, np.busday_offset.
    @classmethod
    def setup_class(self):
        self.a = np.array(["2020-12-31", "2021-01-01", "2021-01-02"], dtype="M")
        self.mask_a = np.array([False, True, False])
        self.ma = Masked(self.a, mask=self.mask_a)
        self.b = np.array([["2021-01-07"], ["2021-01-31"]], dtype="M")
        self.mask_b = np.array([[False], [True]])
        self.mb = Masked(self.b, mask=self.mask_b)

    def test_datetime_as_string(self):
        out = np.datetime_as_string(self.ma)
        expected = np.datetime_as_string(self.a)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, self.mask_a)


@pytest.mark.filterwarnings("ignore:all-nan")
class TestNaNFunctions:
    def setup_class(self):
        self.a = np.array(
            [
                [np.nan, np.nan, 3.0],
                [4.0, 5.0, 6.0],
            ]
        )
        self.mask_a = np.array(
            [
                [True, False, False],
                [False, True, False],
            ]
        )
        self.b = np.arange(1, 7).reshape(2, 3)
        self.mask_b = self.mask_a
        self.ma = Masked(self.a, mask=self.mask_a)
        self.mb = Masked(self.b, mask=self.mask_b)

    def check(self, function, exact_fill_value=None, masked_result=True, **kwargs):
        result = function(self.ma, **kwargs)
        expected_data = function(self.ma.filled(np.nan), **kwargs)
        expected_mask = np.isnan(expected_data)
        if masked_result:
            assert isinstance(result, Masked)
            assert_array_equal(result.mask, expected_mask)
            assert np.all(result == expected_data)
        else:
            assert not isinstance(result, Masked)
            assert_array_equal(result, expected_data)
            assert not np.any(expected_mask)
        out = np.zeros_like(result)
        result2 = function(self.ma, out=out, **kwargs)
        assert result2 is out
        assert_array_equal(result2, result)

    def check_arg(self, function, **kwargs):
        # arg functions do not have an 'out' argument, so just test directly.
        result = function(self.ma, **kwargs)
        assert not isinstance(result, Masked)
        expected = function(self.ma.filled(np.nan), **kwargs)
        assert_array_equal(result, expected)

    def test_nanmin(self):
        self.check(np.nanmin)
        self.check(np.nanmin, axis=0)
        self.check(np.nanmin, axis=1)
        resi = np.nanmin(self.mb, axis=1)
        assert_array_equal(resi.unmasked, np.array([2, 4]))
        assert_array_equal(resi.mask, np.array([False, False]))

    def test_nanmax(self):
        self.check(np.nanmax)

    def test_nanargmin(self):
        self.check_arg(np.nanargmin)
        self.check_arg(np.nanargmin, axis=1)

    def test_nanargmax(self):
        self.check_arg(np.nanargmax)

    def test_nansum(self):
        self.check(np.nansum, masked_result=False)
        resi = np.nansum(self.mb, axis=1)
        assert not isinstance(resi, Masked)
        assert_array_equal(resi, np.array([5, 10]))

    def test_nanprod(self):
        self.check(np.nanprod, masked_result=False)
        resi = np.nanprod(self.mb, axis=1)
        assert not isinstance(resi, Masked)
        assert_array_equal(resi, np.array([6, 24]))

    def test_nancumsum(self):
        self.check(np.nancumsum, masked_result=False)
        resi = np.nancumsum(self.mb, axis=1)
        assert not isinstance(resi, Masked)
        assert_array_equal(resi, np.array([[0, 2, 5], [4, 4, 10]]))

    def test_nancumprod(self):
        self.check(np.nancumprod, masked_result=False)
        resi = np.nancumprod(self.mb, axis=1)
        assert not isinstance(resi, Masked)
        assert_array_equal(resi, np.array([[1, 2, 6], [4, 4, 24]]))

    def test_nanmean(self):
        self.check(np.nanmean)
        resi = np.nanmean(self.mb, axis=1)
        assert_array_equal(resi.unmasked, np.mean(self.mb, axis=1).unmasked)
        assert_array_equal(resi.mask, np.array([False, False]))

    def test_nanvar(self):
        self.check(np.nanvar)
        self.check(np.nanvar, ddof=1)

    def test_nanstd(self):
        self.check(np.nanstd)

    def test_nanmedian(self):
        self.check(np.nanmedian)

    def test_nanquantile(self):
        self.check(np.nanquantile, q=0.5)

    def test_nanpercentile(self):
        self.check(np.nanpercentile, q=50)


untested_functions = set()

if NUMPY_LT_1_23:
    deprecated_functions = {
        # Deprecated, removed in numpy 1.23
        np.asscalar,
        np.alen,
    }
else:
    deprecated_functions = set()

untested_functions |= deprecated_functions
io_functions = {np.save, np.savez, np.savetxt, np.savez_compressed}
untested_functions |= io_functions
poly_functions = {
    np.poly, np.polyadd, np.polyder, np.polydiv, np.polyfit, np.polyint,
    np.polymul, np.polysub, np.polyval, np.roots, np.vander,
}  # fmt: skip
untested_functions |= poly_functions


# Get covered functions
tested_functions = set()
for cov_cls in list(filter(inspect.isclass, locals().values())):
    for k, v in cov_cls.__dict__.items():
        if inspect.isfunction(v) and k.startswith("test"):
            f = k.replace("test_", "")
            if f in all_wrapped_functions:
                tested_functions.add(all_wrapped_functions[f])


def test_basic_testing_completeness():
    assert all_wrapped == (tested_functions | IGNORED_FUNCTIONS | UNSUPPORTED_FUNCTIONS)


@pytest.mark.xfail(reason="coverage not completely set up yet")
def test_testing_completeness():
    assert not tested_functions.intersection(untested_functions)
    assert all_wrapped == (tested_functions | untested_functions)


class TestFunctionHelpersCompleteness:
    @pytest.mark.parametrize(
        "one, two",
        itertools.combinations(
            (
                MASKED_SAFE_FUNCTIONS,
                UNSUPPORTED_FUNCTIONS,
                set(APPLY_TO_BOTH_FUNCTIONS.keys()),
                set(DISPATCHED_FUNCTIONS.keys()),
            ),
            2,
        ),
    )
    def test_no_duplicates(self, one, two):
        assert not one.intersection(two)

    def test_all_included(self):
        included_in_helpers = (
            MASKED_SAFE_FUNCTIONS
            | UNSUPPORTED_FUNCTIONS
            | set(APPLY_TO_BOTH_FUNCTIONS.keys())
            | set(DISPATCHED_FUNCTIONS.keys())
        )
        assert all_wrapped == included_in_helpers

    @pytest.mark.xfail(reason="coverage not completely set up yet")
    def test_ignored_are_untested(self):
        assert IGNORED_FUNCTIONS == untested_functions

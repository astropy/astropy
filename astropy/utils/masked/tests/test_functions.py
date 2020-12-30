# Licensed under a 3-clause BSD style license - see LICENSE.rst
import itertools

import pytest
import numpy as np
from numpy.testing import assert_array_equal

from astropy.utils.compat import NUMPY_LT_1_18, NUMPY_LT_1_20
from astropy.units.tests.test_quantity_non_ufuncs import (
    get_wrapped_functions, CoverageMeta)

from ..core import Masked, MaskedNDArray
from ..function_helpers import (MASKED_SAFE_FUNCTIONS,
                                APPLY_TO_BOTH_FUNCTIONS,
                                DISPATCHED_FUNCTIONS,
                                IGNORED_FUNCTIONS,
                                UNSUPPORTED_FUNCTIONS)

from .test_masked import (MaskedArraySetup, QuantitySetup, LongitudeSetup,
                          assert_masked_equal)


all_wrapped_functions = get_wrapped_functions(np)
all_wrapped = set(all_wrapped_functions.values())


CoverageMeta.covered = set()


# We first do some explicit, more difficult comparisons, and then
# try to run through all numpy functions using just simple tests.

class TestMaskedArrayConcatenation(MaskedArraySetup, metaclass=CoverageMeta):
    def test_concatenate(self):
        mb = self.mb[np.newaxis]
        concat_a_b = np.concatenate((self.ma, mb), axis=0)
        expected_data = np.concatenate((self.a, self.b[np.newaxis]), axis=0)
        expected_mask = np.concatenate((self.mask_a, self.mask_b[np.newaxis]),
                                       axis=0)
        assert_array_equal(concat_a_b.unmasked, expected_data)
        assert_array_equal(concat_a_b.mask, expected_mask)

    def test_concatenate_not_all_masked(self):
        mb = self.mb[np.newaxis]
        concat_a_b = np.concatenate((self.a, mb), axis=0)
        expected_data = np.concatenate((self.a, self.b[np.newaxis]), axis=0)
        expected_mask = np.concatenate((np.zeros(self.a.shape, bool),
                                        self.mask_b[np.newaxis]), axis=0)
        assert_array_equal(concat_a_b.unmasked, expected_data)
        assert_array_equal(concat_a_b.mask, expected_mask)

    @pytest.mark.parametrize('obj', (1, slice(2, 3)))
    def test_insert(self, obj):
        mc_in_a = np.insert(self.ma, obj, self.mc, axis=-1)
        expected = Masked(np.insert(self.a, obj, self.c, axis=-1),
                          np.insert(self.mask_a, obj, self.mask_c, axis=-1))
        assert_masked_equal(mc_in_a, expected)

    def test_append(self):
        mc_to_a = np.append(self.ma, self.mc, axis=-1)
        expected = Masked(np.append(self.a, self.c, axis=-1),
                          np.append(self.mask_a, self.mask_c, axis=-1))
        assert_masked_equal(mc_to_a, expected)


class TestMaskedQuantityConcatenation(TestMaskedArrayConcatenation,
                                      QuantitySetup):
    pass


class TestMaskedLongitudeConcatenation(TestMaskedArrayConcatenation,
                                       LongitudeSetup):
    pass


class TestMaskedArrayBroadcast(MaskedArraySetup):
    def test_broadcast_to(self):
        shape = self.ma.shape
        ba = np.broadcast_to(self.mb, shape, subok=True)
        assert ba.shape == shape
        assert ba.mask.shape == shape
        expected = Masked(np.broadcast_to(self.mb.unmasked, shape, subok=True),
                          np.broadcast_to(self.mb.mask, shape, subok=True))
        assert_masked_equal(ba, expected)

    def test_broadcast_arrays(self):
        mb = np.broadcast_arrays(self.ma, self.mb, self.mc, subok=True)
        b = np.broadcast_arrays(self.a, self.b, self.c, subok=True)
        bm = np.broadcast_arrays(self.mask_a, self.mask_b, self.mask_c)
        for mb_, b_, bm_ in zip(mb, b, bm):
            assert_array_equal(mb_.unmasked, b_)
            assert_array_equal(mb_.mask, bm_)

    def test_broadcast_arrays_not_all_masked(self):
        mb = np.broadcast_arrays(self.a, self.mb, self.c, subok=True)
        assert_array_equal(mb[0], self.a)
        expected1 = np.broadcast_to(self.mb, self.a.shape, subok=True)
        assert_masked_equal(mb[1], expected1)
        expected2 = np.broadcast_to(self.c, self.a.shape, subok=True)
        assert_array_equal(mb[2], expected2)

    def test_broadcast_arrays_subok_false(self):
        # subok affects ndarray subclasses but not masking itself.
        mb = np.broadcast_arrays(self.ma, self.mb, self.mc, subok=False)
        assert all(type(mb_.unmasked) is np.ndarray for mb_ in mb)
        b = np.broadcast_arrays(self.a, self.b, self.c, subok=False)
        mask_b = np.broadcast_arrays(self.mask_a, self.mask_b,
                                     self.mask_c, subok=False)
        for mb_, b_, mask_ in zip(mb, b, mask_b):
            assert_array_equal(mb_.unmasked, b_)
            assert_array_equal(mb_.mask, mask_)


class TestMaskedQuantityBroadcast(TestMaskedArrayBroadcast, QuantitySetup):
    pass


class TestMaskedLongitudeBroadcast(TestMaskedArrayBroadcast, LongitudeSetup):
    pass


class TestMaskedArrayCalculation(MaskedArraySetup):
    @pytest.mark.parametrize('n,axis', [(1, -1), (2, -1), (1, 0)])
    def test_diff(self, n, axis):
        mda = np.diff(self.ma, n=n, axis=axis)
        expected_data = np.diff(self.a, n, axis)
        nan_mask = np.zeros_like(self.a)
        nan_mask[self.ma.mask] = np.nan
        expected_mask = np.isnan(np.diff(nan_mask, n=n, axis=axis))
        assert_array_equal(mda.unmasked, expected_data)
        assert_array_equal(mda.mask, expected_mask)

    def test_diff_explicit(self):
        ma = Masked(np.arange(8.),
                    [True, False, False, False, False, True, False, False])
        mda = np.diff(ma)
        assert np.all(mda.unmasked == 1.)
        assert np.all(mda.mask ==
                      [True, False, False, False, True, True, False])
        mda = np.diff(ma, n=2)
        assert np.all(mda.unmasked == 0.)
        assert np.all(mda.mask == [True, False, False, True, True, True])


class TestMaskedQuantityCalculation(TestMaskedArrayCalculation, QuantitySetup):
    pass


class TestMaskedLongitudeCalculation(TestMaskedArrayCalculation,
                                     LongitudeSetup):
    pass


class TestMaskedArraySorting(MaskedArraySetup):
    @pytest.mark.parametrize('axis', [-1, 0])
    def test_lexsort1(self, axis):
        ma_lexsort = np.lexsort((self.ma,), axis=axis)
        filled = self.a.copy()
        filled[self.mask_a] = 9e9
        expected_data = filled.argsort(axis)
        assert_array_equal(ma_lexsort, expected_data)

    @pytest.mark.parametrize('axis', [-1, 0])
    def test_lexsort2(self, axis):
        mb = np.broadcast_to(-self.mb, self.ma.shape).copy()
        mamb_lexsort = np.lexsort((self.ma, mb), axis=axis)
        filled_a = self.ma.unmask(9e9)
        filled_b = mb.unmask(9e9)
        expected_ab = np.lexsort((filled_a, filled_b), axis=axis)
        assert_array_equal(mamb_lexsort, expected_ab)
        mbma_lexsort = np.lexsort((mb, self.ma), axis=axis)
        expected_ba = np.lexsort((filled_b, filled_a), axis=axis)
        assert_array_equal(mbma_lexsort, expected_ba)
        mbma_lexsort2 = np.lexsort(np.stack([mb, self.ma], axis=0), axis=axis)
        assert_array_equal(mbma_lexsort2, expected_ba)

    @pytest.mark.parametrize('axis', [-1, 0])
    def test_lexsort_mix(self, axis):
        mb = np.broadcast_to(-self.mb, self.ma.shape).copy()
        mamb_lexsort = np.lexsort((self.a, mb), axis=axis)
        filled_b = mb.unmask(9e9)
        expected_ab = np.lexsort((self.a, filled_b), axis=axis)
        assert_array_equal(mamb_lexsort, expected_ab)
        mbma_lexsort = np.lexsort((mb, self.a), axis=axis)
        expected_ba = np.lexsort((filled_b, self.a), axis=axis)
        assert_array_equal(mbma_lexsort, expected_ba)
        mbma_lexsort2 = np.lexsort(np.stack([mb, self.a], axis=0), axis=axis)
        assert_array_equal(mbma_lexsort2, expected_ba)


# Now run through all functions, with simple tests.

class BasicTestSetup(MaskedArraySetup, metaclass=CoverageMeta):
    def check(self, func, *args, **kwargs):
        o = func(self.ma, *args, **kwargs)
        expected = Masked(func(self.a, *args, **kwargs),
                          mask=func(self.mask_a, *args, **kwargs))
        assert_array_equal(o.unmasked, expected.unmasked)
        assert_array_equal(o.mask, expected.mask)


class NoMaskTestSetup(MaskedArraySetup, metaclass=CoverageMeta):
    def check(self, func, *args, **kwargs):
        o = func(self.ma, *args, **kwargs)
        expected = func(self.a, *args, **kwargs)
        assert_array_equal(o, expected)


class InvariantMaskTestSetup(MaskedArraySetup, metaclass=CoverageMeta):
    def check(self, func, *args, **kwargs):
        o = func(self.ma, *args, **kwargs)
        expected = func(self.a, *args, **kwargs)
        assert_array_equal(o.unmasked, expected)
        assert_array_equal(o.mask, self.mask_a)


class TestShapeInformation(BasicTestSetup):
    # alen is deprecated in Numpy 1.8
    if NUMPY_LT_1_18:
        def test_alen(self):
            assert np.alen(self.ma) == 3

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


class TestArgFunctions(MaskedArraySetup, metaclass=CoverageMeta):
    def check(self, function, *args, fill_value=np.nan, **kwargs):
        o = function(self.ma, *args, **kwargs)
        a_filled = self.ma.unmask(fill_value=fill_value)
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
        self.check(np.nonzero, fill_value=0.)

    def test_argwhere(self):
        self.check(np.argwhere, fill_value=0.)

    @pytest.mark.xfail(reason='not implemented yet')
    def test_argpartition(self):
        self.check(np.argpartition, 2, fill_value=np.inf)

    def test_flatnonzero(self):
        self.check(np.flatnonzero, fill_value=0.)


class TestAlongAxis(BasicTestSetup):
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

    @pytest.mark.parametrize('axis', (0, 1))
    def test_apply_along_axis(self, axis):
        out = np.apply_along_axis(np.square, axis, self.ma)
        expected = np.apply_along_axis(np.square, axis, self.a)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, self.mask_a)

    @pytest.mark.parametrize('axes', ((1,), (0,), (0, 1)))
    def test_apply_over_axes(self, axes):
        def function(x, axis):
            return np.mean(np.square(x), axis)

        out = np.apply_over_axes(function, self.ma, axes)
        expected = self.ma
        for axis in axes:
            expected = (expected**2).mean(axis, keepdims=True)
        assert_array_equal(out.unmasked, expected.unmasked)
        assert_array_equal(out.mask, expected.mask)


class TestIndicesFrom(NoMaskTestSetup):
    def setup(self):
        self.a = np.arange(9).reshape(3, 3)
        self.mask_a = np.eye(3, dtype=bool)
        self.ma = Masked(self.a, self.mask_a)

    def test_diag_indices_from(self):
        self.check(np.diag_indices_from)

    def test_triu_indices_from(self):
        self.check(np.triu_indices_from)

    def test_tril_indices_from(self):
        self.check(np.tril_indices_from)


class TestRealImag(InvariantMaskTestSetup, metaclass=CoverageMeta):
    def setup(self):
        self.a = np.array([1+2j, 3+4j])
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

    @pytest.mark.parametrize('value', [0.5, Masked(0.5, mask=True), np.ma.masked])
    def test_full_like(self, value):
        o = np.full_like(self.ma, value)
        expected = Masked(np.zeros_like(self.a))
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

    @pytest.mark.xfail(reason='not implemented yet.')
    def test_trim_zeros(self):
        o = np.trim_zeros(self.ma.ravel())
        expected = 0
        assert np.all(o == expected)

    def test_roll(self):
        self.check(np.roll, 1)
        self.check(np.roll, 1, axis=0)

    def test_take(self):
        self.check(np.take, [0, 1], axis=1)
        self.check(np.take, 1)


class TestSettingParts(MaskedArraySetup, metaclass=CoverageMeta):
    def test_put(self):
        ma = self.ma.copy()
        np.put(ma, [0, 2], Masked([50, 150], [False, True]))
        expected = self.a.copy()
        np.put(expected, [0, 2], [50, 150])
        expected_mask = self.mask_a.copy()
        np.put(expected_mask, [0, 2], [False, True])
        assert_array_equal(ma.unmasked, expected)
        assert_array_equal(ma.mask, expected_mask)

    def test_putmask(self):
        ma = self.ma.flatten()
        mask = [True, False, False, False, True, False]
        values = Masked(np.arange(100, 650, 100),
                        mask=[False, True, True, True, False, False])
        np.putmask(ma, mask, values)
        expected = self.a.flatten()
        np.putmask(expected, mask, values.unmasked)
        expected_mask = self.mask_a.flatten()
        np.putmask(expected_mask, mask, values.mask)
        assert_array_equal(ma.unmasked, expected)
        assert_array_equal(ma.mask, expected_mask)

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

    def test_copyto(self):
        ma = self.ma.flatten()
        mask = [True, False, False, False, True, False]
        values = Masked(np.arange(100, 650, 100),
                        mask=[False, True, True, True, False, False])
        np.copyto(ma, values, where=mask)
        expected = self.a.flatten()
        np.copyto(expected, values.unmasked, where=mask)
        expected_mask = self.mask_a.flatten()
        np.copyto(expected_mask, values.mask, where=mask)
        assert_array_equal(ma.unmasked, expected)
        assert_array_equal(ma.mask, expected_mask)

    @pytest.mark.parametrize('value', [0.25, np.ma.masked])
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


class TestConcatenate(MaskedArraySetup, metaclass=CoverageMeta):
    # More tests at TestMaskedArrayConcatenation above.
    def check(self, func, *args, **kwargs):
        ma_list = kwargs.pop('ma_list', [self.ma, self.ma])
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

        out = Masked(np.empty((4, 3)))
        result = np.concatenate([self.ma, self.ma], out=out)
        assert out is result
        expected = np.concatenate([self.a, self.a])
        expected_mask = np.concatenate([self.mask_a, self.mask_a])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

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

        out = np.block([[0., Masked(1., True)],
                        [Masked(1, False), Masked(2, False)]])
        expected = np.array([[0, 1.], [1, 2]])
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
        out = np.insert(self.ma.flatten(), (1, 1),
                        Masked([50., 25.], mask=[True, False]))
        expected = np.insert(self.a.flatten(), (1, 1), [50., 25.])
        expected_mask = np.insert(self.mask_a.flatten(), (1, 1), [True, False])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)


class TestSplit(metaclass=CoverageMeta):
    def setup(self):
        self.ma = Masked(np.arange(54.).reshape(3, 3, 6))
        self.ma.mask[1, 1, 1] = True
        self.ma.mask[0, 1, 4] = True
        self.ma.mask[1, 2, 5] = True

    def check(self, func, *args, **kwargs):
        out = func(self.ma, *args, **kwargs)
        expected = func(self.ma.unmasked, *args, **kwargs)
        expected_mask = func(self.ma.mask, *args, **kwargs)
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


class TestMethodLikes(MaskedArraySetup, metaclass=CoverageMeta):
    def check(self, function, method=None, *args, **kwargs):
        if method is None:
            method = function.__name__

        o = function(self.ma, *args, **kwargs)
        x = getattr(self.ma, method)(*args, **kwargs)
        assert_array_equal(o.unmasked, x.unmasked)
        assert_array_equal(o.mask, x.mask)

    def test_amax(self):
        self.check(np.amax, method='max')

    def test_amin(self):
        self.check(np.amin, method='min')

    def test_sum(self):
        self.check(np.sum)

    @pytest.mark.xfail(reason='need to implement accumulate')
    def test_cumsum(self):
        self.check(np.cumsum)

    def test_any(self):
        self.check(np.any)

    def test_all(self):
        self.check(np.all)

    def test_sometrue(self):
        self.check(np.sometrue, method='any')

    def test_alltrue(self):
        self.check(np.alltrue, method='all')

    def test_prod(self):
        self.check(np.prod)

    def test_product(self):
        self.check(np.product, method='prod')

    @pytest.mark.xfail(reason='need to implement accumulate')
    def test_cumprod(self):
        self.check(np.cumprod)

    @pytest.mark.xfail(reason='need to implement accumulate')
    def test_cumproduct(self):
        self.check(np.cumproduct, method='cumprod')

    def test_ptp(self):
        self.check(np.ptp)
        self.check(np.ptp, axis=0)

    def test_round_(self):
        self.check(np.round_, method='round')

    def test_around(self):
        self.check(np.around, method='round')

    @pytest.mark.xfail(reason='need to implement clip properly')
    def test_clip(self):
        self.check(np.clip, 2., 4.)
        self.check(np.clip, self.mb, self.mc)


class TestUfuncLike(InvariantMaskTestSetup):
    def test_fix(self):
        self.check(np.fix)

    def test_angle(self):
        a = np.array([1+0j, 0+1j, 1+1j, 0+0j])
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
        out = np.where([True, False, True], self.ma, 1000.)
        expected = np.where([True, False, True], self.a, 1000.)
        expected_mask = np.where([True, False, True], self.mask_a, False)
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_choose(self):
        a = np.array([0, 1]).reshape((2, 1))
        out = np.choose(a, (self.ma, self.mb))
        expected = np.choose(a, (self.a, self.b))
        expected_mask = np.choose(a, (self.mask_a, self.mask_b))
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_choose_masked(self):
        ma = Masked(np.array([0, 1]), mask=[True, False]).reshape((2, 1))
        out = ma.choose((self.ma, self.mb))
        expected = np.choose(ma.unmasked, (self.a, self.b))
        expected_mask = np.choose(ma.unmasked, (self.mask_a, self.mask_b)) | ma.mask
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    @pytest.mark.parametrize('default', [-1., np.ma.masked, Masked(-1, mask=True)])
    def test_select(self, default):
        a, mask_a, ma = self.a, self.mask_a, self.ma
        out = np.select([a < 1.5, a > 3.5], [ma, ma+1], default=default)
        expected = np.select([a < 1.5, a > 3.5], [a, a+1],
                             default=-1 if default is not np.ma.masked else 0)
        expected_mask = np.select([a < 1.5, a > 3.5], [mask_a, mask_a],
                                  default=getattr(default, 'mask', False))
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_real_if_close(self):
        a = np.array([1+0j, 0+1j, 1+1j, 0+0j])
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


class TestUfuncLikeTests(metaclass=CoverageMeta):
    def setup(self):
        self.a = np.array([[-np.inf, +np.inf, np.nan, 3., 4.]]*2)
        self.mask_a = np.array([[False]*5, [True]*4+[False]])
        self.ma = Masked(self.a, mask=self.mask_a)
        self.b = np.array([[3.0001], [3.9999]])
        self.mask_b = np.array([[True], [False]])
        self.mb = Masked(self.b, mask=self.mask_b)

    def check(self, func):
        out = func(self.ma)
        expected = func(self.a)
        assert type(out) is MaskedNDArray
        assert out.dtype.kind == 'b'
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, self.mask_a)
        assert not np.may_share_memory(out.mask, self.mask_a)

    def test_isposinf(self):
        self.check(np.isposinf)

    def test_isneginf(self):
        self.check(np.isneginf)

    def test_isreal(self):
        self.check(np.isreal)
        o = np.isreal(Masked([1. + 1j], mask=False))
        assert not o.unmasked and not o.mask
        o = np.isreal(Masked([1. + 1j], mask=True))
        assert not o.unmasked and o.mask

    def test_iscomplex(self):
        self.check(np.iscomplex)
        o = np.iscomplex(Masked([1. + 1j], mask=False))
        assert o.unmasked and not o.mask
        o = np.iscomplex(Masked([1. + 1j], mask=True))
        assert o.unmasked and o.mask

    def test_isclose(self):
        out = np.isclose(self.ma, self.mb, atol=0.01)
        expected = np.isclose(self.ma, self.mb, atol=0.01)
        expected_mask = self.mask_a | self.mask_b
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_allclose(self):
        out = np.allclose(self.ma, self.mb, atol=0.01)
        expected = np.isclose(self.ma, self.mb,
                              atol=0.01)[self.mask_a | self.mask_b].all()
        assert_array_equal(out, expected)


class TestSpaceFunctions(metaclass=CoverageMeta):
    def setup_class(self):
        self.a = np.arange(1., 7.).reshape(2, 3)
        self.mask_a = np.array([[True, False, False],
                                [False, True, False]])
        self.ma = Masked(self.a, mask=self.mask_a)
        self.b = np.array([2.5, 10., 3.])
        self.mask_b = np.array([False, True, False])
        self.mb = Masked(self.b, mask=self.mask_b)

    def check(self, function, *args, **kwargs):
        out = function(self.ma, self.mb, 5)
        expected = function(self.a, self.b, 5)
        expected_mask = np.broadcast_to(self.mask_a | self.mask_b,
                                        expected.shape).copy()
        # TODO: make implementation that also ensures start point mask
        # is determined just by start point?
        expected_mask[-1] = self.mask_b
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)

    def test_linspace(self):
        self.check(np.linspace, 5)

    def test_logspace(self):
        self.check(np.logspace, 10)

    def test_geomspace(self):
        self.check(np.geomspace, 5)


class TestInterpolationFunctions(MaskedArraySetup, metaclass=CoverageMeta):
    def test_interp(self):
        xp = np.arange(5.)
        fp = np.array([1., 5., 6., 19., 20.])
        mask_fp = np.array([False, False, False, True, False])
        mfp = Masked(fp, mask=mask_fp)
        x = np.array([1.5, 17.])
        mask_x = np.array([False, True])
        mx = Masked(x, mask=mask_x)
        out = np.interp(mx, xp, mfp)
        expected = np.interp(x, xp[mask_fp], fp[mask_fp])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, mask_x)

    def test_piecewise(self):
        condlist = [self.a < 1, self.a >= 1]
        out = np.piecewise(self.ma, condlist, [Masked(-1, mask=True), 1.])
        expected = np.piecewise(self.a, condlist, [-1, 1.])
        expected_mask = np.piecewise(self.mask_a, condlist, [True, False])
        assert_array_equal(out.unmasked, expected)
        assert_array_equal(out.mask, expected_mask)
        condlist2 = [self.a < 1, self.a >= 3]
        out2 = np.piecewise(self.ma, condlist2,
                            [Masked(-1, True), 1, lambda x: Masked(np.full_like(x, 2.),
                                                                   mask=~x.mask)])
        expected = np.piecewise(self.a, condlist2, [-1, 1, 2])
        expected_mask = np.piecewise(self.mask_a, condlist2,
                                     [True, False, lambda x: ~x])
        assert_array_equal(out2.unmasked, expected)
        assert_array_equal(out2.mask, expected_mask)


class TestMemoryFunctions(MaskedArraySetup, metaclass=CoverageMeta):
    def test_shares_memory(self):
        assert np.shares_memory(self.ma, self.ma.unmasked)
        assert not np.shares_memory(self.ma, self.ma.mask)

    def test_may_share_memory(self):
        assert np.may_share_memory(self.ma, self.ma.unmasked)
        assert not np.may_share_memory(self.ma, self.ma.mask)


untested_functions = set()
if NUMPY_LT_1_20:
    financial_functions = {f for f in all_wrapped_functions.values()
                           if f in np.lib.financial.__dict__.values()}
    untested_functions |= financial_functions

deprecated_functions = {
    np.asscalar
    }
if NUMPY_LT_1_18:
    deprecated_functions |= {np.rank}
else:
    deprecated_functions |= {np.alen}

untested_functions |= deprecated_functions
io_functions = {np.save, np.savez, np.savetxt, np.savez_compressed}
untested_functions |= io_functions

poly_functions = {
    np.poly, np.polyadd, np.polyder, np.polydiv, np.polyfit, np.polyint,
    np.polymul, np.polysub, np.polyval, np.roots, np.vander
    }
untested_functions |= poly_functions


@pytest.mark.xfail(reason='not yet complete')
def test_testing_completeness():
    assert not CoverageMeta.covered.intersection(untested_functions)
    assert all_wrapped == (CoverageMeta.covered | untested_functions)


class TestFunctionHelpersCompleteness:
    @pytest.mark.parametrize('one, two', itertools.combinations(
        (MASKED_SAFE_FUNCTIONS,
         UNSUPPORTED_FUNCTIONS,
         set(APPLY_TO_BOTH_FUNCTIONS.keys()),
         set(DISPATCHED_FUNCTIONS.keys())), 2))
    def test_no_duplicates(self, one, two):
        assert not one.intersection(two)

    @pytest.mark.xfail(reason='not yet complete')
    def test_all_included(self):
        included_in_helpers = (MASKED_SAFE_FUNCTIONS |
                               UNSUPPORTED_FUNCTIONS |
                               set(APPLY_TO_BOTH_FUNCTIONS.keys()) |
                               set(DISPATCHED_FUNCTIONS.keys()))
        assert all_wrapped == included_in_helpers

    # untested_function is created using all_wrapped_functions
    def test_ignored_are_untested(self):
        assert IGNORED_FUNCTIONS == untested_functions

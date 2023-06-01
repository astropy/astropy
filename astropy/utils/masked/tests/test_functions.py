# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test numpy functions and ufuncs on Masked arrays and quantities.

The tests here are fairly detailed but do not aim for complete
coverage.  Complete coverage of all numpy functions is done
with less detailed tests in test_function_helpers.
"""
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.units import Quantity
from astropy.utils.compat.numpycompat import NUMPY_LT_1_25
from astropy.utils.masked.core import Masked

from .test_masked import (
    LongitudeSetup,
    MaskedArraySetup,
    QuantitySetup,
    assert_masked_equal,
)


class MaskedUfuncTests(MaskedArraySetup):
    @pytest.mark.parametrize(
        "ufunc", (np.add, np.subtract, np.divide, np.arctan2, np.minimum)
    )
    @pytest.mark.parametrize("a, b", [("ma", "mb"), ("ma", "b"), ("a", "mb")])
    def test_2op_ufunc(self, ufunc, a, b):
        a, b = getattr(self, a), getattr(self, b)
        mask_a = getattr(a, "mask", np.zeros(a.shape, bool))
        mask_b = getattr(b, "mask", np.zeros(b.shape, bool))
        result = ufunc(a, b)
        expected_data = ufunc(self.a, self.b)
        expected_mask = mask_a | mask_b
        # Note: assert_array_equal also checks type, i.e., that, e.g.,
        # Longitude decays into an Angle.
        assert_array_equal(result.unmasked, expected_data)
        assert_array_equal(result.mask, expected_mask)

        out = Masked(np.zeros_like(result.unmasked))
        result2 = ufunc(a, b, out=out)
        assert result2 is out
        assert_masked_equal(result2, result)

    @pytest.mark.parametrize("base_mask", [True, False])
    def test_ufunc_inplace_where(self, base_mask):
        # Construct base filled with -9 and base_mask (copying to get unit/class).
        base = self.ma.copy()
        base.unmasked.view(np.ndarray)[...] = -9.0
        base._mask[...] = base_mask
        out = base.copy()
        where = np.array([[True, False, False], [False, True, False]])
        result = np.add(self.ma, self.mb, out=out, where=where)
        # Direct checks.
        assert np.all(result.unmasked[~where] == base.unmasked[0, 0])
        assert np.all(result.unmasked[where] == (self.a + self.b)[where])
        # Full comparison.
        expected = base.unmasked.copy()
        np.add(self.a, self.b, out=expected, where=where)
        expected_mask = base.mask.copy()
        np.logical_or(self.mask_a, self.mask_b, out=expected_mask, where=where)
        assert_array_equal(result.unmasked, expected)
        assert_array_equal(result.mask, expected_mask)

    @pytest.mark.parametrize("base_mask", [True, False])
    def test_ufunc_inplace_masked_where(self, base_mask):
        base = self.ma.copy()
        base.unmasked.view(np.ndarray)[...] = -9.0
        base._mask[...] = base_mask
        out = base.copy()
        where = Masked(
            [[True, False, True], [False, False, True]],
            mask=[[True, False, False], [True, False, True]],
        )
        result = np.add(self.ma, self.mb, out=out, where=where)
        # Direct checks.
        assert np.all(result.unmasked[~where.unmasked] == base.unmasked[0, 0])
        assert np.all(
            result.unmasked[where.unmasked] == (self.a + self.b)[where.unmasked]
        )
        assert np.all(result.mask[where.mask])
        assert np.all(result.mask[~where.mask & ~where.unmasked] == base.mask[0, 0])
        assert np.all(
            result.mask[~where.mask & where.unmasked]
            == (self.mask_a | self.mask_b)[~where.mask & where.unmasked]
        )
        # Full comparison.
        expected = base.unmasked.copy()
        np.add(self.a, self.b, out=expected, where=where.unmasked)
        expected_mask = base.mask.copy()
        np.logical_or(self.mask_a, self.mask_b, out=expected_mask, where=where.unmasked)
        expected_mask |= where.mask
        assert_array_equal(result.unmasked, expected)
        assert_array_equal(result.mask, expected_mask)

    def test_ufunc_inplace_no_masked_input(self):
        a_b = np.add(self.a, self.b)
        out = Masked(np.zeros_like(a_b))
        result = np.add(self.a, self.b, out=out)
        assert result is out
        assert_array_equal(result.unmasked, a_b)
        assert_array_equal(result.mask, np.zeros(a_b.shape, bool))

    def test_ufunc_inplace_error(self):
        # Output is not masked.
        out = np.zeros(self.ma.shape)
        with pytest.raises(TypeError):
            np.add(self.ma, self.mb, out=out)

    @pytest.mark.xfail(NUMPY_LT_1_25, reason="masked where not supported in numpy<1.25")
    def test_ufunc_inplace_error_masked_where(self):
        # Input and output are not masked, but where is.
        # Note: prior to numpy 1.25, we cannot control this.
        out = self.a.copy()
        with pytest.raises(TypeError):
            np.add(self.a, self.b, out=out, where=Masked(True, mask=True))

    @pytest.mark.parametrize("ufunc", (np.add.outer, np.minimum.outer))
    @pytest.mark.parametrize("a, b", [("ma", "mb"), ("ma", "b"), ("a", "mb")])
    def test_2op_ufunc_outer(self, ufunc, a, b):
        a, b = getattr(self, a), getattr(self, b)
        mask_a = getattr(a, "mask", np.zeros(a.shape, bool))
        mask_b = getattr(b, "mask", np.zeros(b.shape, bool))
        result = ufunc(a, b)
        expected_data = ufunc(self.a, self.b)
        expected_mask = np.logical_or.outer(mask_a, mask_b)
        # Note: assert_array_equal also checks type, i.e., that, e.g.,
        # Longitude decays into an Angle.
        assert_array_equal(result.unmasked, expected_data)
        assert_array_equal(result.mask, expected_mask)

        out = Masked(np.zeros_like(result.unmasked))
        result2 = ufunc(a, b, out=out)
        assert result2 is out
        assert_masked_equal(result2, result)

    @pytest.mark.parametrize("ufunc", (np.add.outer, np.minimum.outer))
    def test_2op_ufunc_outer_no_masked_input(self, ufunc):
        expected_data = ufunc(self.a, self.b)
        out = Masked(np.zeros_like(expected_data), True)
        result = ufunc(self.a, self.b, out=out)
        assert_array_equal(out.unmasked, expected_data)
        assert_array_equal(out.mask, np.zeros(out.shape, dtype=bool))

    def test_3op_ufunc(self):
        ma_mb = np.clip(self.ma, self.b, self.c)
        expected_data = np.clip(self.a, self.b, self.c)
        expected_mask = self.mask_a
        assert_array_equal(ma_mb.unmasked, expected_data)
        assert_array_equal(ma_mb.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_add_reduce(self, axis):
        ma_reduce = np.add.reduce(self.ma, axis=axis)
        expected_data = np.add.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)

        out = Masked(np.zeros_like(ma_reduce.unmasked), np.ones_like(ma_reduce.mask))
        ma_reduce2 = np.add.reduce(self.ma, axis=axis, out=out)
        assert ma_reduce2 is out
        assert_masked_equal(ma_reduce2, ma_reduce)

    def test_add_reduce_no_masked_input(self):
        a_reduce = np.add.reduce(self.a, axis=0)
        out = Masked(np.zeros_like(a_reduce), np.ones(a_reduce.shape, bool))
        result = np.add.reduce(self.a, axis=0, out=out)
        assert result is out
        assert_array_equal(out.unmasked, a_reduce)
        assert_array_equal(out.mask, np.zeros(a_reduce.shape, bool))

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_minimum_reduce(self, axis):
        ma_reduce = np.minimum.reduce(self.ma, axis=axis)
        expected_data = np.minimum.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)

    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_maximum_reduce(self, axis):
        ma_reduce = np.maximum.reduce(self.ma, axis=axis)
        expected_data = np.maximum.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)


class TestMaskedArrayUfuncs(MaskedUfuncTests):
    # multiply.reduce does not work with units, so test only for plain array.
    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_multiply_reduce(self, axis):
        ma_reduce = np.multiply.reduce(self.ma, axis=axis)
        expected_data = np.multiply.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)

    def test_ufunc_not_implemented_for_other(self):
        """
        If the unmasked operation returns NotImplemented, this
        should lead to a TypeError also for the masked version.
        """
        a = np.array([1, 2])
        b = 3 * u.m
        with pytest.raises(TypeError):
            a & b

        ma = Masked(a)
        with pytest.raises(TypeError):
            ma & b


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


class TestMaskedArrayConcatenation(MaskedArraySetup):
    def test_concatenate(self):
        mb = self.mb[np.newaxis]
        concat_a_b = np.concatenate((self.ma, mb), axis=0)
        expected_data = np.concatenate((self.a, self.b[np.newaxis]), axis=0)
        expected_mask = np.concatenate((self.mask_a, self.mask_b[np.newaxis]), axis=0)
        assert_array_equal(concat_a_b.unmasked, expected_data)
        assert_array_equal(concat_a_b.mask, expected_mask)

    def test_concatenate_not_all_masked(self):
        mb = self.mb[np.newaxis]
        concat_a_b = np.concatenate((self.a, mb), axis=0)
        expected_data = np.concatenate((self.a, self.b[np.newaxis]), axis=0)
        expected_mask = np.concatenate(
            (np.zeros(self.a.shape, bool), self.mask_b[np.newaxis]), axis=0
        )
        assert_array_equal(concat_a_b.unmasked, expected_data)
        assert_array_equal(concat_a_b.mask, expected_mask)

    @pytest.mark.parametrize("obj", (1, slice(2, 3)))
    def test_insert(self, obj):
        mc_in_a = np.insert(self.ma, obj, self.mc, axis=-1)
        expected = Masked(
            np.insert(self.a, obj, self.c, axis=-1),
            np.insert(self.mask_a, obj, self.mask_c, axis=-1),
        )
        assert_masked_equal(mc_in_a, expected)

    def test_insert_masked_obj(self):
        with pytest.raises(TypeError):
            np.insert(self.ma, Masked(1, mask=False), self.mc, axis=-1)

    def test_append(self):
        mc_to_a = np.append(self.ma, self.mc, axis=-1)
        expected = Masked(
            np.append(self.a, self.c, axis=-1),
            np.append(self.mask_a, self.mask_c, axis=-1),
        )
        assert_masked_equal(mc_to_a, expected)


class TestMaskedQuantityConcatenation(TestMaskedArrayConcatenation, QuantitySetup):
    pass


class TestMaskedLongitudeConcatenation(TestMaskedArrayConcatenation, LongitudeSetup):
    pass


class TestMaskedArrayBroadcast(MaskedArraySetup):
    def test_broadcast_to(self):
        shape = self.ma.shape
        ba = np.broadcast_to(self.mb, shape, subok=True)
        assert ba.shape == shape
        assert ba.mask.shape == shape
        expected = Masked(
            np.broadcast_to(self.mb.unmasked, shape, subok=True),
            np.broadcast_to(self.mb.mask, shape, subok=True),
        )
        assert_masked_equal(ba, expected)

    def test_broadcast_to_using_apply(self):
        # Partially just to ensure we cover the relevant part of _apply.
        shape = self.ma.shape
        ba = self.mb._apply(np.broadcast_to, shape=shape, subok=True)
        assert ba.shape == shape
        assert ba.mask.shape == shape
        expected = Masked(
            np.broadcast_to(self.mb.unmasked, shape, subok=True),
            np.broadcast_to(self.mb.mask, shape, subok=True),
        )
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
        mask_b = np.broadcast_arrays(self.mask_a, self.mask_b, self.mask_c, subok=False)
        for mb_, b_, mask_ in zip(mb, b, mask_b):
            assert_array_equal(mb_.unmasked, b_)
            assert_array_equal(mb_.mask, mask_)


class TestMaskedQuantityBroadcast(TestMaskedArrayBroadcast, QuantitySetup):
    pass


class TestMaskedLongitudeBroadcast(TestMaskedArrayBroadcast, LongitudeSetup):
    pass


class TestMaskedArrayCalculation(MaskedArraySetup):
    @pytest.mark.parametrize("n,axis", [(1, -1), (2, -1), (1, 0)])
    def test_diff(self, n, axis):
        mda = np.diff(self.ma, n=n, axis=axis)
        expected_data = np.diff(self.a, n, axis)
        nan_mask = np.zeros_like(self.a)
        nan_mask[self.ma.mask] = np.nan
        expected_mask = np.isnan(np.diff(nan_mask, n=n, axis=axis))
        assert_array_equal(mda.unmasked, expected_data)
        assert_array_equal(mda.mask, expected_mask)

    def test_diff_explicit(self):
        ma = Masked(
            np.arange(8.0), [True, False, False, False, False, True, False, False]
        )
        mda = np.diff(ma)
        assert np.all(mda.unmasked == 1.0)
        assert np.all(mda.mask == [True, False, False, False, True, True, False])
        mda = np.diff(ma, n=2)
        assert np.all(mda.unmasked == 0.0)
        assert np.all(mda.mask == [True, False, False, True, True, True])


class TestMaskedQuantityCalculation(TestMaskedArrayCalculation, QuantitySetup):
    pass


class TestMaskedLongitudeCalculation(TestMaskedArrayCalculation, LongitudeSetup):
    pass


class TestMaskedArraySorting(MaskedArraySetup):
    @pytest.mark.parametrize("axis", [-1, 0])
    def test_lexsort1(self, axis):
        ma_lexsort = np.lexsort((self.ma,), axis=axis)
        filled = self.a.copy()
        filled[self.mask_a] = 9e9
        expected_data = filled.argsort(axis)
        assert_array_equal(ma_lexsort, expected_data)

    @pytest.mark.parametrize("axis", [-1, 0])
    def test_lexsort2(self, axis):
        mb = np.broadcast_to(-self.mb, self.ma.shape).copy()
        mamb_lexsort = np.lexsort((self.ma, mb), axis=axis)
        filled_a = self.ma.filled(9e9)
        filled_b = mb.filled(9e9)
        expected_ab = np.lexsort((filled_a, filled_b), axis=axis)
        assert_array_equal(mamb_lexsort, expected_ab)
        mbma_lexsort = np.lexsort((mb, self.ma), axis=axis)
        expected_ba = np.lexsort((filled_b, filled_a), axis=axis)
        assert_array_equal(mbma_lexsort, expected_ba)
        mbma_lexsort2 = np.lexsort(np.stack([mb, self.ma], axis=0), axis=axis)
        assert_array_equal(mbma_lexsort2, expected_ba)

    @pytest.mark.parametrize("axis", [-1, 0])
    def test_lexsort_mix(self, axis):
        mb = np.broadcast_to(-self.mb, self.ma.shape).copy()
        mamb_lexsort = np.lexsort((self.a, mb), axis=axis)
        filled_b = mb.filled(9e9)
        expected_ab = np.lexsort((self.a, filled_b), axis=axis)
        assert_array_equal(mamb_lexsort, expected_ab)
        mbma_lexsort = np.lexsort((mb, self.a), axis=axis)
        expected_ba = np.lexsort((filled_b, self.a), axis=axis)
        assert_array_equal(mbma_lexsort, expected_ba)
        mbma_lexsort2 = np.lexsort(np.stack([mb, self.a], axis=0), axis=axis)
        assert_array_equal(mbma_lexsort2, expected_ba)

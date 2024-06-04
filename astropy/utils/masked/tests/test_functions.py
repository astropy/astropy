# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test numpy functions and ufuncs on Masked arrays and quantities.

The tests here are fairly detailed but do not aim for complete
coverage.  Complete coverage of all numpy functions is done
with less detailed tests in test_function_helpers.
"""

# We generally call the ufunc in the tests, since those can take
# all ufunc arguments (like axes), but also test whether we can
# mask the exceptions and warnings from the wrappers in erfa itself.
import erfa
import erfa.ufunc as erfa_ufunc
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from astropy import units as u
from astropy.units import Quantity
from astropy.utils import minversion
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

    def test_multi_op_ufunc(self):
        mask = [True, False, False]
        iy = Masked([2000, 2001, 2002], mask=mask)
        im = Masked([1, 2, 3], mask=mask)
        idy = Masked([10, 20, 25], mask=mask)
        ihr = Masked([11, 12, 13], mask=[False, False, True])
        imn = np.array([50, 51, 52])
        isc = np.array([12.5, 13.6, 14.7])
        result = erfa_ufunc.dtf2d("utc", iy, im, idy, ihr, imn, isc)
        # Also test scalar
        result0 = erfa_ufunc.dtf2d("utc", iy[0], im[0], idy[0], ihr[0], imn[0], isc[0])
        expected = erfa_ufunc.dtf2d(
            "utc", iy.unmasked, im.unmasked, idy.unmasked, ihr.unmasked, imn, isc
        )
        expected_mask = np.array([True, False, True])
        for res, res0, exp in zip(result, result0, expected):
            assert_array_equal(res.unmasked, exp)
            assert_array_equal(res.mask, expected_mask)
            assert res0.unmasked == exp[0]
            assert res0.mask == expected_mask[0]

    def test_erfa_pdp_with_out(self):
        p = np.arange(9.0).reshape(3, 3)
        mp = Masked(p)
        mp.mask[1, 2] = True
        out = Masked(np.empty(3), mask=True)
        result = erfa_ufunc.pdp(mp, mp, out=out)
        assert result is out
        assert_array_equal(result.unmasked, erfa_ufunc.pdp(p, p))
        assert_array_equal(result.mask, [False, True, False])
        # With axes just for the inputs.
        axes = [0, 1]
        result2 = erfa_ufunc.pdp(mp, mp, out=out, axes=axes)
        assert result2 is out
        assert_array_equal(result2.unmasked, erfa_ufunc.pdp(p, p, axes=axes))
        assert_array_equal(result2.mask, [False, True, True])

    @pytest.mark.parametrize("kwargs", [{}, dict(axis=0), dict(axes=[0])])
    def test_erfa_p2s_with_out(self, kwargs):
        p = np.arange(9.0).reshape(3, 3)
        mp = Masked(p)
        mp.mask[1, 2] = True
        # Outputs are theta, phi, r.
        outs = tuple(Masked(np.empty(3), mask=True) for _ in range(3))
        masks = tuple(out.mask for out in outs)
        results = erfa_ufunc.p2s(mp, out=outs, **kwargs)
        assert len(results) == 3
        expected = erfa_ufunc.p2s(mp.unmasked, **kwargs)
        expected_mask = mp.mask.any(0 if kwargs else -1)
        for a, b, m, x in zip(results, outs, masks, expected):
            assert a is b
            assert a.mask is m
            assert_array_equal(a.unmasked, x)
            assert_array_equal(a.mask, expected_mask)

    def test_erfa_rxp(self):
        # Regression tests for gh-16116
        m = Masked(np.eye(3))
        v = Masked(np.arange(6).reshape(2, 3))
        rxp1 = erfa_ufunc.rxp(m, v)
        exp = erfa_ufunc.rxp(m.unmasked, v.unmasked)
        assert_array_equal(rxp1.unmasked, exp)
        assert_array_equal(rxp1.mask, False)
        v.mask[0, 0] = True
        rxp2 = erfa_ufunc.rxp(m, v)
        assert_array_equal(rxp2.unmasked, exp)
        assert_array_equal(rxp2.mask, [[True] * 3, [False] * 3])
        m.mask[1, 1] = True
        v.mask[...] = False
        rxp3 = erfa_ufunc.rxp(m, v)
        assert_array_equal(rxp3.unmasked, exp)
        assert_array_equal(rxp3.mask, True)

    def test_erfa_rxr_axes(self):
        m1 = Masked(np.arange(27.0).reshape(3, 3, 3))
        m2 = Masked(np.arange(-27.0, 0.0).reshape(3, 3, 3))
        rxr1 = erfa_ufunc.rxr(m1, m2)
        exp = erfa_ufunc.rxr(m1.unmasked, m2.unmasked)
        assert_array_equal(rxr1.unmasked, exp)
        assert_array_equal(rxr1.mask, False)
        m1.mask[0, 1, 2] = True
        rxr2 = erfa_ufunc.rxr(m1, m2)
        assert_array_equal(rxr2.unmasked, exp)
        assert np.all(rxr2.mask == [[[True]], [[False]], [[False]]])
        axes = [(0, 2), (-2, -1), (0, 1)]
        rxr3 = erfa_ufunc.rxr(m1, m2, axes=axes)
        exp3 = erfa_ufunc.rxr(m1.unmasked, m2.unmasked, axes=axes)
        assert_array_equal(rxr3.unmasked, exp3)
        assert np.all(rxr3.mask == [False, True, False])

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
        assert_array_equal(out.mask, False)
        # Also try with where (which should have different path)
        where = np.array([[True, False, False], [True, True, False]])
        a_reduce2 = np.add.reduce(self.a, axis=0, where=where)
        result2 = np.add.reduce(self.a, axis=0, out=out, where=where)
        assert result2 is out
        assert_array_equal(out.unmasked, a_reduce2)
        assert_array_equal(out.mask, [False, False, True])

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
    # Similarly, modf only works for dimensionless, but we are using it just
    # to check multiple outputs get the right mask.
    @pytest.mark.parametrize("axis", (0, 1, None))
    def test_multiply_reduce(self, axis):
        ma_reduce = np.multiply.reduce(self.ma, axis=axis)
        expected_data = np.multiply.reduce(self.a, axis=axis)
        expected_mask = np.logical_or.reduce(self.ma.mask, axis=axis)
        assert_array_equal(ma_reduce.unmasked, expected_data)
        assert_array_equal(ma_reduce.mask, expected_mask)

    def test_ufunc_two_out(self):
        out0 = np.empty_like(self.ma)
        out0_mask = out0.mask
        out1 = np.empty_like(self.ma)
        out1_mask = out1.mask
        res0, res1 = np.modf(self.ma, out=(out0, out1))
        assert res0 is out0
        assert res1 is out1
        assert out0.mask is out0_mask
        assert out1.mask is out1_mask
        assert_array_equal(out0.mask, self.mask_a)
        assert_array_equal(out1.mask, self.mask_a)

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


class TestStructuredUfuncs:
    """Test with structure dtypes, using erfa ufuncs."""

    def test_erfa_d2tf_tf2d(self):
        mask = np.array([True, False, False])
        days = Masked([0.25, 0.875, 0.0625], mask=mask)
        sign, ihmsf = erfa_ufunc.d2tf(3, days)
        assert_array_equal(sign.mask["sign"], mask)
        sign = sign.view("S1")  # Like is done by the erfa wrapper.
        assert_array_equal(sign.mask, mask)
        for name in ihmsf.dtype.names:
            assert_array_equal(ihmsf[name].mask, mask)

        # Check roundtrip.
        check, stat = erfa_ufunc.tf2d(
            sign, ihmsf["h"], ihmsf["m"], ihmsf["s"] + ihmsf["f"]
        )
        assert_allclose(check.unmasked, days, atol=1e-3 / 24 / 3600)
        assert_array_equal(check.mask, mask)
        assert_array_equal(stat.unmasked, 0)
        assert_array_equal(stat.mask, mask)

    def test_erfa_astrom(self):
        mask = np.array([True, False, False])
        jd2 = Masked([0, 0.401182685, 0.5], mask=mask)
        astrom, eo = erfa_ufunc.apci13(2456165.5, jd2)
        assert_array_equal(eo.mask, mask)
        for n in astrom.dtype.names:
            # .T for multi-element fields.
            assert np.all(astrom[n].mask.T == mask)

        along = np.array([0.125, 0.25, 0.35])
        # Not going to worry about different masks for different elements.
        # In principle aper could propagate just what it needs.
        astrom["along"] = Masked([0.125, 0.25, 0.35], mask)
        astrom2 = erfa_ufunc.aper(Masked(np.ones(3), [False, True, False]), astrom)
        assert_array_equal(astrom2["eral"].unmasked, along + 1.0)
        mask2 = mask | [False, True, False]
        for n in astrom2.dtype.names:
            # .T for multi-element fields.
            assert np.all(astrom2[n].mask.T == mask2)

    def test_erfa_atioq(self):
        # Regression test for gh-16123, using test from erfa.
        astrom, _ = erfa_ufunc.apio13(
            2456384.5,
            0.969254051,
            0.1550675,
            -0.527800806,
            -1.2345856,
            2738.0,
            2.47230737e-7,
            1.82640464e-6,
            731.0,
            12.8,
            0.59,
            0.55,
        )
        astrom = Masked(astrom)
        ri = 2.710121572969038991
        di = 0.1729371367218230438
        aob, zob, hob, dob, rob = erfa_ufunc.atioq(ri, di, astrom)
        assert isinstance(aob, Masked)
        # Really should not need to check the values, since done
        # in units/tests/test_quantity_erfa_ufuncs, but why not...
        assert_allclose(aob, 0.9233952224895122499e-1, atol=1e-12, rtol=0)
        assert_allclose(zob, 1.407758704513549991, atol=1e-12, rtol=0)
        assert_allclose(hob, -0.9247619879881698140e-1, atol=1e-12, rtol=0)
        assert_allclose(dob, 0.1717653435756234676, atol=1e-12, rtol=0)
        assert_allclose(rob, 2.710085107988480746, atol=1e-12, rtol=0)


def test_erfa_no_warnings_on_masked_entries():
    # Erfa warns for invalid inputs for some routines.
    msg = 'ERFA function "tf2d" yielded {count} of "ihour outside range 0-23"'
    ihour1 = [25, 26, 10]
    with pytest.warns(erfa.ErfaWarning, match=msg.format(count=2)):
        res1 = erfa.tf2d("+", ihour1, 0, 0.0)
    # But will not if they are masked.
    mask = [True, True, False]
    ihour2 = Masked(ihour1, mask)
    res2 = erfa.tf2d("+", ihour2, 0, 0.0)
    assert_array_equal(res2.unmasked, res1)
    assert_array_equal(res2.mask, mask)
    # And will count correctly.
    mask = [True, False, False]
    ihour3 = Masked(ihour1, mask)
    count = 1 if minversion(erfa, "2.0.1.1", inclusive=False) else "—"
    with pytest.warns(erfa.ErfaWarning, match=msg.format(count=count)):
        res3 = erfa.tf2d("+", ihour3, 0, 0.0)
    assert_array_equal(res3.unmasked, res1)
    assert_array_equal(res3.mask, mask)


def test_erfa_no_exceptions_on_masked_entries():
    # Erfa raises exceptions for invalid inputs in some routines.
    iday1 = [25, 30]
    with pytest.raises(erfa.ErfaError, match="bad day"):
        erfa.dat(2000, 2, iday1, 0.0)
    # But will not if they are masked.
    mask = [False, True]
    iday2 = Masked(iday1, mask)
    res = erfa.dat(2000, 2, iday2, 0.0)
    exp1 = erfa.dat(2000, 2, iday1[0], 0.0)
    assert_array_equal(res.unmasked[0], exp1)
    assert_array_equal(res.mask, mask)

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import itertools

import pytest
import numpy as np
from numpy.testing import assert_array_equal

from astropy.utils.masked.core import Masked, MaskedNDArray

from .test_masked import (MaskedArraySetup, QuantitySetup, LongitudeSetup,
                          assert_masked_equal)


class TestMaskedArrayConcatenation(MaskedArraySetup):
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

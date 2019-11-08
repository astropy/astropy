# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest
import numpy as np
from numpy.testing import assert_array_equal

from ..core import Masked

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

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test numpy functions and ufuncs on Distribution.

The tests here are fairly detailed but do not aim for complete
coverage.  Complete coverage of all numpy functions will be done
in test_function_helpers.

TODO: start test_function_helpers once more functions are supported.
"""

import numpy as np
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.uncertainty import Distribution


class ArraySetup:
    @classmethod
    def setup_class(cls):
        cls.a = (
            np.array([[[0.0]], [[10.0]]])
            + np.array([[0.0], [1.0], [2.0]])
            + np.arange(4.0) / 10.0
        )
        cls.b = -(np.arange(3.0, 6.0)[:, np.newaxis] + np.arange(4.0) / 10.0)
        cls.da = Distribution(cls.a)
        cls.db = Distribution(cls.b)
        cls.c = np.array([[200.0], [300.0]])


class QuantitySetup(ArraySetup):
    @classmethod
    def setup_class(cls):
        super().setup_class()
        cls.a <<= u.m
        cls.b <<= u.km
        cls.da <<= u.m
        cls.db <<= u.km
        cls.c <<= u.Mm


class TestConcatenation(ArraySetup):
    def test_concatenate(self):
        # Concatenate needs consistent shapes.
        db = self.db[np.newaxis]
        concat_a_b = np.concatenate((self.da, db), axis=0)
        expected_distr = np.concatenate((self.a, self.b[np.newaxis]), axis=0)
        assert_array_equal(concat_a_b.distribution, expected_distr)

    def test_concatenate_not_all_distribution(self):
        concat_c_a = np.concatenate((self.c, self.da), axis=1)
        assert isinstance(concat_c_a, Distribution)
        c_bcst = np.broadcast_to(
            self.c[..., np.newaxis], self.c.shape + (self.da.n_samples,), subok=True
        )

        expected_distr = np.concatenate((c_bcst, self.a), axis=1)
        assert_array_equal(concat_c_a.distribution, expected_distr)


class TestQuantityDistributionConcatenation(TestConcatenation, QuantitySetup):
    pass


class TestBroadcast(ArraySetup):
    def test_broadcast_to(self):
        shape = self.da.shape
        ba = np.broadcast_to(self.db, shape, subok=True)
        assert ba.shape == shape
        expected_distr = np.broadcast_to(self.b, self.a.shape, subok=True)
        assert_array_equal(ba.distribution, expected_distr)

    def test_broadcast_arrays(self):
        bda, bdb, bdc = np.broadcast_arrays(self.da, self.db, self.c, subok=True)
        assert type(bda) is type(bdb) is type(self.da)
        assert type(bdc) is type(self.c)
        ba, bb = np.broadcast_arrays(self.a, self.b, subok=True)
        bc = np.broadcast_to(self.c, self.da.shape, subok=True)
        assert_array_equal(bda.distribution, ba)
        assert_array_equal(bdb.distribution, bb)
        assert_array_equal(bdc, bc)

    def test_broadcast_arrays_subok_false(self):
        # subok affects ndarray subclasses but not distribution itself.
        bda, bdb, bdc = np.broadcast_arrays(self.da, self.db, self.c, subok=False)
        assert type(bda.distribution) is type(bdb.distribution) is np.ndarray
        assert type(bdc) is np.ndarray
        ba, bb = np.broadcast_arrays(self.a, self.b, subok=False)
        bc = np.broadcast_to(self.c, self.da.shape, subok=False)
        assert_array_equal(bda.distribution, ba)
        assert_array_equal(bdb.distribution, bb)
        assert_array_equal(bdc, bc)


class TestQuantityDistributionBroadcast(TestBroadcast, QuantitySetup):
    pass

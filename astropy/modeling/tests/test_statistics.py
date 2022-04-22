# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Module to test statistic functions
"""
import numpy as np
# pylint: disable=invalid-name
import pytest
from numpy.testing import assert_almost_equal

from astropy.modeling.models import Identity, Mapping
from astropy.modeling.statistic import leastsquare, leastsquare_1d, leastsquare_2d, leastsquare_3d


class TestLeastSquare_XD:
    """Tests for leastsquare with pre-specified number of dimensions."""

    @classmethod
    def setup_class(cls):
        cls.model1D = Identity(n_inputs=1)
        cls.model2D = Identity(n_inputs=2) | Mapping((0,), n_inputs=2)
        cls.model3D = Identity(n_inputs=3) | Mapping((0,), n_inputs=3)

        cls.data = cls.x = cls.y = cls.z = np.linspace(0, 10, num=100)
        cls.lsq_exp = 0

    def test_1d_no_weights(self):
        lsq = leastsquare_1d(self.data, self.model1D, None, self.x)
        assert_almost_equal(lsq, self.lsq_exp)

    def test_1d_with_weights(self):
        lsq = leastsquare_1d(self.data, self.model1D, np.ones(100), self.x)
        assert_almost_equal(lsq, self.lsq_exp)

    def test_2d_no_weights(self):
        lsq = leastsquare_2d(self.data, self.model2D, None, self.x, self.y)
        assert_almost_equal(lsq, self.lsq_exp)

    def test_2d_with_weights(self):
        lsq = leastsquare_2d(
            self.data, self.model2D, np.ones(100), self.x, self.y
        )
        assert_almost_equal(lsq, self.lsq_exp)

    def test_3d_no_weights(self):
        lsq = leastsquare_3d(
            self.data, self.model3D, None, self.x, self.y, self.z
        )
        assert_almost_equal(lsq, self.lsq_exp)

    def test_3d_with_weights(self):
        lsq = leastsquare_3d(
            self.data, self.model3D, np.ones(100), self.x, self.y, self.z
        )
        assert_almost_equal(lsq, self.lsq_exp)


class TestLeastSquare_ND:
    """Tests for leastsquare."""

    @classmethod
    def setup_class(cls):
        cls.model1D = Identity(n_inputs=1)
        cls.model3D = Identity(n_inputs=3) | Mapping((0,), n_inputs=3)

        cls.data = cls.x = cls.y = cls.z = np.linspace(0, 10, num=100)
        cls.lsq_exp = 0

    def test_1d_no_weights(self):
        lsq = leastsquare(self.data, self.model1D, None, self.x)
        assert_almost_equal(lsq, self.lsq_exp)

    def test_1d_with_weights(self):
        lsq = leastsquare(self.data, self.model1D, np.ones(100), self.x)
        assert_almost_equal(lsq, self.lsq_exp)

    def test_3d_no_weights(self):
        lsq = leastsquare(
            self.data, self.model3D, None, self.x, self.y, self.z
        )
        assert_almost_equal(lsq, self.lsq_exp)

    def test_3d_with_weights(self):
        lsq = leastsquare(
            self.data, self.model3D, np.ones(100), self.x, self.y, self.z
        )
        assert_almost_equal(lsq, self.lsq_exp)

    def test_shape_mismatch(self):
        with pytest.raises(ValueError):
            leastsquare(0, self.model1D, None, self.x)

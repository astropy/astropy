# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module tests fitting and model evaluation with various inputs
"""
from __future__ import division
import numpy as np
from .. import models
from .. import fitting
from numpy.testing import utils
from ...tests.helper import pytest

try:
    from scipy import optimize  # pylint: disable=W0611
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


model1d_params = [
    (models.Polynomial1D, [2]),
    (models.Legendre1D, [2]),
    (models.Chebyshev1D, [2]),
    (models.Shift, [2]),
    (models.Scale, [2])
]

model2d_params = [
    (models.Polynomial2D, [2]),
    (models.Legendre2D, [1, 2]),
    (models.Chebyshev2D, [1, 2])
]


class TestInputType(object):

    """
    This class tests that models accept numbers, lists and arrays.

    Add new models to one of the lists above to test for this.
    """
    def setup_class(self):
        self.x = 5.3
        self.y = 6.7
        self.x1 = np.arange(1, 10, .1)
        self.y1 = np.arange(1, 10, .1)
        self.x2, self.y2 = np.mgrid[:10, :8]

    @pytest.mark.parametrize(('model', 'params'), model1d_params)
    def test_input1D(self, model, params):
        m = model(*params)
        m(self.x)
        m(self.x1)
        m(self.x2)

    @pytest.mark.parametrize(('model', 'params'), model2d_params)
    def test_input2D(self, model, params):
        m = model(*params)
        m(self.x, self.y)
        m(self.x1, self.y1)
        m(self.x2, self.y2)


class TestFitting(object):

    """
    test various input options to fitting routines
    """
    def setup_class(self):
        self.x1 = np.arange(10)
        self.x, self.y = np.mgrid[:10, :10]

    def test_linear_fitter_1set(self):
        """
        1 set 1D x, 1pset
        """
        expected = np.array([0, 1, 1, 1])
        p1 = models.Polynomial1D(3)
        p1.parameters = [0, 1, 1, 1]
        y1 = p1(self.x1)
        pfit = fitting.LinearLSQFitter()
        model = pfit(p1, self.x1, y1)
        utils.assert_allclose(model.parameters, expected, atol=10 ** (-7))

    def test_linear_fitter_Nset(self):
        """
        1 set 1D x, 2 sets 1D y, 2 param_sets
        """
        expected = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
        p1 = models.Polynomial1D(3, param_dim=2)
        p1.parameters = [0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0]
        params = {}
        for i in range(4):
            params[p1.param_names[i]] = [i, i]
        p1 = models.Polynomial1D(3, param_dim=2, **params)
        y1 = p1(self.x1)
        pfit = fitting.LinearLSQFitter()
        model = pfit(p1, self.x1, y1)
        utils.assert_allclose(model.param_sets, expected, atol=10 ** (-7))

    def test_linear_fitter_1dcheb(self):
        """
        1 pset, 1 set 1D x, 1 set 1D y, Chebyshev 1D polynomial
        """
        expected = np.array(
            [[2817.2499999999995,
              4226.6249999999991,
              1680.7500000000009,
              273.37499999999926]]).T
        ch1 = models.Chebyshev1D(3)
        ch1.parameters = [0, 1, 2, 3]
        y1 = ch1(self.x1)
        pfit = fitting.LinearLSQFitter()
        model = pfit(ch1, self.x1, y1)
        utils.assert_allclose(model.param_sets, expected, atol=10 ** (-2))

    def test_linear_fitter_1dlegend(self):
        """
        1 pset, 1 set 1D x, 1 set 1D y, Legendre 1D polynomial
        """
        expected = np.array(
            [[1925.5000000000011,
              3444.7500000000005,
              1883.2500000000014,
              364.4999999999996]]).T
        leg1 = models.Legendre1D(3)
        leg1.parameters = [1, 2, 3, 4]
        y1 = leg1(self.x1)
        pfit = fitting.LinearLSQFitter()
        model = pfit(leg1, self.x1, y1)
        utils.assert_allclose(model.param_sets, expected, atol=10 ** (-12))

    def test_linear_fitter_1set2d(self):
        p2 = models.Polynomial2D(2)
        p2.parameters = [0, 1, 2, 3, 4, 5]
        expected = [0, 1, 2, 3, 4, 5]
        z = p2(self.x, self.y)
        pfit = fitting.LinearLSQFitter()
        model = pfit(p2, self.x, self.y, z)
        utils.assert_allclose(model.parameters, expected, atol=10 ** (-12))
        utils.assert_allclose(model(self.x, self.y), z, atol=10 ** (-12))

    def test_wrong_numpset(self):
        """
        A ValueError is raised if a 1 data set (1d x, 1d y) is fit
        with a model with multiple parameter sets.
        """
        with pytest.raises(ValueError):
            p1 = models.Polynomial1D(5)
            y1 = p1(self.x1)
            p1 = models.Polynomial1D(5, param_dim=2)
            pfit = fitting.LinearLSQFitter()
            model = pfit(p1, self.x1, y1)

    def test_wrong_pset(self):
        """
        A case of 1 set of x and multiple sets of y and parameters
        """
        expected = np.array([[1., 0],
                             [1, 1],
                             [1, 2],
                             [1, 3],
                             [1, 4],
                             [1, 5]])
        p1 = models.Polynomial1D(5, param_dim=2)
        params = {}
        for i in range(6):
            params[p1.param_names[i]] = [1, i]
        p1 = models.Polynomial1D(5, param_dim=2, **params)
        y1 = p1(self.x1)
        pfit = fitting.LinearLSQFitter()
        model = pfit(p1, self.x1, y1)
        utils.assert_allclose(model.param_sets, expected, atol=10 ** (-7))

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_nonlinear_lsqt_1set_1d(self):
        """
        1 set 1D x, 1 set 1D y, 1 pset NonLinearFitter
        """
        g1 = models.Gaussian1D(10, mean=3, stddev=.2)
        y1 = g1(self.x1)
        gfit = fitting.NonLinearLSQFitter()
        model = gfit(g1, self.x1, y1)
        utils.assert_allclose(model.parameters, [10, 3, .2])

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_nonlinear_lsqt_Nset_1d(self):
        """
        1 set 1D x, 1 set 1D y, 2 param_sets, NonLinearFitter
        """
        with pytest.raises(ValueError):
            g1 = models.Gaussian1D([10.2, 10], mean=[3, 3.2], stddev=[.23, .2])
            y1 = g1(self.x1)
            gfit = fitting.NonLinearLSQFitter()
            model = gfit(g1, self.x1, y1)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_nonlinear_lsqt_1set_2d(self):
        """
        1 set 2d x, 1set 2D y, 1 pset, NonLinearFitter
        """
        g2 = models.Gaussian2D(10, x_mean=3, y_mean=4, x_stddev=.3, y_stddev=.2, theta=0)
        z = g2(self.x, self.y)
        gfit = fitting.NonLinearLSQFitter()
        model = gfit(g2, self.x, self.y, z)
        utils.assert_allclose(model.parameters, [10, 3, 4, .3, .2, 0])

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_nonlinear_lsqt_Nset_2d(self):
        """
         1 set 2d x, 1set 2D y, 2 param_sets, NonLinearFitter
        """
        with pytest.raises(ValueError):
            g2 = models.Gaussian2D([10, 10], [3, 3], [4, 4], x_stddev=[.3, .3],
                                   y_stddev=[.2, .2], theta=[0, 0])
            z = g2(self.x.flatten(), self.y.flatten())
            gfit = fitting.NonLinearLSQFitter()
            model = gfit(g2, self.x, self.y, z)


class TestEvaluation(object):

    """
    test various input options to model evaluation

    TestFitting actually covers evaluation of polynomials
    """
    def setup_class(self):
        self.x1 = np.arange(20)
        self.x, self.y = np.mgrid[:10, :10]

    def test_non_linear_NYset(self):
        """
        This case covers:
            N param sets , 1 set 1D x --> N 1D y data
        """
        g1 = models.Gaussian1D([10, 10], [3, 3], [.2, .2])
        y1 = g1(self.x1)
        utils.assert_equal((y1[:, 0] - y1[:, 1]).nonzero(), (np.array([]),))

    def test_non_linear_NXYset(self):
        """
        This case covers: N param sets , N sets 1D x --> N N sets 1D y data
        """
        g1 = models.Gaussian1D([10, 10], [3, 3], [.2, .2])
        xx = np.array([self.x1, self.x1])
        y1 = g1(xx.T)
        utils.assert_allclose(y1[:, 0], y1[:, 1], atol=10 ** (-12))

    def test_p1_1set_1pset(self):
        """
        1 data set, 1 pset, Polynomial1D
        """
        p1 = models.Polynomial1D(4)
        y1 = p1(self.x1)
        assert y1.shape == (20,)

    def test_p1_nset_npset(self):
        """
        N data sets, N param_sets, Polynomial1D
        """
        p1 = models.Polynomial1D(4, param_dim=2)
        y1 = p1(np.array([self.x1, self.x1]).T)
        assert y1.shape == (20, 2)
        utils.assert_allclose(y1[:, 0], y1[:, 1], atol=10 ** (-12))

    def test_p2_1set_1pset(self):
        """
        1 pset, 1 2D data set, Polynomial2D
        """
        p2 = models.Polynomial2D(5)
        z = p2(self.x, self.y)
        assert z.shape == (10, 10)

    def test_p2_nset_npset(self):
        """
        N param_sets, N 2D data sets, Poly2d
        """
        p2 = models.Polynomial2D(5, param_dim=2)
        xx = np.array([self.x, self.x])
        yy = np.array([self.y, self.y])
        z = p2(xx, yy)
        assert z.shape == (2, 10, 10)

    def test_nset_domain(self):
        """
        Polynomial evaluation of multiple data sets with different domain
        """
        xx = np.array([self.x1, self.x1]).T
        xx[0, 0] = 100
        xx[1, 0] = 100
        xx[2, 0] = 99
        p1 = models.Polynomial1D(5, param_dim=2)
        yy = p1(xx)
        x1 = xx[:, 0]
        x2 = xx[:, 1]
        p1 = models.Polynomial1D(5)
        utils.assert_allclose(p1(x1), yy[:, 0], atol=10 ** (-12))
        p1 = models.Polynomial1D(5)
        utils.assert_allclose(p1(x2), yy[:, 1], atol=10 ** (-12))

    def test_evaluate_gauss2d(self):
        cov = np.array([[1., 0.8], [0.8, 3]])
        g = models.Gaussian2D(1., 5., 4., cov_matrix=cov)
        x, y = np.mgrid[:10, :10]
        g(x, y)

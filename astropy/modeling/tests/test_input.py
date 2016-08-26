# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module tests fitting and model evaluation with various inputs
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from numpy.testing.utils import assert_allclose

from .. import models
from .. import fitting
from ..core import Model, FittableModel, Fittable1DModel
from ..parameters import Parameter
from ...tests.helper import pytest
from ...extern.six.moves import range

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
        self.y2, self.x2 = np.mgrid[:10, :8]

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
    """Test various input options to fitting routines."""

    def setup_class(self):
        self.x1 = np.arange(10)
        self.y, self.x = np.mgrid[:10, :10]

    def test_linear_fitter_1set(self):
        """1 set 1D x, 1pset"""

        expected = np.array([0, 1, 1, 1])
        p1 = models.Polynomial1D(3)
        p1.parameters = [0, 1, 1, 1]
        y1 = p1(self.x1)
        pfit = fitting.LinearLSQFitter()
        model = pfit(p1, self.x1, y1)
        assert_allclose(model.parameters, expected, atol=10 ** (-7))

    def test_linear_fitter_Nset(self):
        """1 set 1D x, 2 sets 1D y, 2 param_sets"""

        expected = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
        p1 = models.Polynomial1D(3, n_models=2)
        p1.parameters = [0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0]
        params = {}
        for i in range(4):
            params[p1.param_names[i]] = [i, i]
        p1 = models.Polynomial1D(3, model_set_axis=0, **params)
        y1 = p1(self.x1, model_set_axis=False)
        pfit = fitting.LinearLSQFitter()
        model = pfit(p1, self.x1, y1)
        assert_allclose(model.param_sets, expected, atol=10 ** (-7))

    def test_linear_fitter_1dcheb(self):
        """1 pset, 1 set 1D x, 1 set 1D y, Chebyshev 1D polynomial"""

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
        assert_allclose(model.param_sets, expected, atol=10 ** (-2))

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
        assert_allclose(model.param_sets, expected, atol=10 ** (-12))

    def test_linear_fitter_1set2d(self):
        p2 = models.Polynomial2D(2)
        p2.parameters = [0, 1, 2, 3, 4, 5]
        expected = [0, 1, 2, 3, 4, 5]
        z = p2(self.x, self.y)
        pfit = fitting.LinearLSQFitter()
        model = pfit(p2, self.x, self.y, z)
        assert_allclose(model.parameters, expected, atol=10 ** (-12))
        assert_allclose(model(self.x, self.y), z, atol=10 ** (-12))

    def test_wrong_numpset(self):
        """
        A ValueError is raised if a 1 data set (1d x, 1d y) is fit
        with a model with multiple parameter sets.
        """

        with pytest.raises(ValueError):
            p1 = models.Polynomial1D(5)
            y1 = p1(self.x1)
            p1 = models.Polynomial1D(5, n_models=2)
            pfit = fitting.LinearLSQFitter()
            model = pfit(p1, self.x1, y1)

    def test_wrong_pset(self):
        """A case of 1 set of x and multiple sets of y and parameters."""

        expected = np.array([[1., 0],
                             [1, 1],
                             [1, 2],
                             [1, 3],
                             [1, 4],
                             [1, 5]])
        p1 = models.Polynomial1D(5, n_models=2)
        params = {}
        for i in range(6):
            params[p1.param_names[i]] = [1, i]
        p1 = models.Polynomial1D(5, model_set_axis=0, **params)
        y1 = p1(self.x1, model_set_axis=False)
        pfit = fitting.LinearLSQFitter()
        model = pfit(p1, self.x1, y1)
        assert_allclose(model.param_sets, expected, atol=10 ** (-7))

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_nonlinear_lsqt_1set_1d(self):
        """1 set 1D x, 1 set 1D y, 1 pset NonLinearFitter"""

        g1 = models.Gaussian1D(10, mean=3, stddev=.2)
        y1 = g1(self.x1)
        gfit = fitting.LevMarLSQFitter()
        model = gfit(g1, self.x1, y1)
        assert_allclose(model.parameters, [10, 3, .2])

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_nonlinear_lsqt_Nset_1d(self):
        """1 set 1D x, 1 set 1D y, 2 param_sets, NonLinearFitter"""

        with pytest.raises(ValueError):
            g1 = models.Gaussian1D([10.2, 10], mean=[3, 3.2], stddev=[.23, .2],
                                   n_models=2)
            y1 = g1(self.x1, model_set_axis=False)
            gfit = fitting.LevMarLSQFitter()
            model = gfit(g1, self.x1, y1)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_nonlinear_lsqt_1set_2d(self):
        """1 set 2d x, 1set 2D y, 1 pset, NonLinearFitter"""

        g2 = models.Gaussian2D(10, x_mean=3, y_mean=4, x_stddev=.3,
                               y_stddev=.2, theta=0)
        z = g2(self.x, self.y)
        gfit = fitting.LevMarLSQFitter()
        model = gfit(g2, self.x, self.y, z)
        assert_allclose(model.parameters, [10, 3, 4, .3, .2, 0])

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_nonlinear_lsqt_Nset_2d(self):
        """1 set 2d x, 1set 2D y, 2 param_sets, NonLinearFitter"""

        with pytest.raises(ValueError):
            g2 = models.Gaussian2D([10, 10], [3, 3], [4, 4], x_stddev=[.3, .3],
                                   y_stddev=[.2, .2], theta=[0, 0], n_models=2)
            z = g2(self.x.flatten(), self.y.flatten())
            gfit = fitting.LevMarLSQFitter()
            model = gfit(g2, self.x, self.y, z)


class TestEvaluation(object):
    """
    Test various input options to model evaluation

    TestFitting actually covers evaluation of polynomials
    """

    def setup_class(self):
        self.x1 = np.arange(20)
        self.y, self.x = np.mgrid[:10, :10]

    def test_non_linear_NYset(self):
        """
        This case covers:
            N param sets , 1 set 1D x --> N 1D y data
        """

        g1 = models.Gaussian1D([10, 10], [3, 3], [.2, .2], n_models=2)
        y1 = g1(self.x1, model_set_axis=False)
        assert np.all((y1[0, :] - y1[1, :]).nonzero() == np.array([]))

    def test_non_linear_NXYset(self):
        """
        This case covers: N param sets , N sets 1D x --> N N sets 1D y data
        """

        g1 = models.Gaussian1D([10, 10], [3, 3], [.2, .2], n_models=2)
        xx = np.array([self.x1, self.x1])
        y1 = g1(xx)
        assert_allclose(y1[:, 0], y1[:, 1], atol=10 ** (-12))

    def test_p1_1set_1pset(self):
        """1 data set, 1 pset, Polynomial1D"""

        p1 = models.Polynomial1D(4)
        y1 = p1(self.x1)
        assert y1.shape == (20,)

    def test_p1_nset_npset(self):
        """N data sets, N param_sets, Polynomial1D"""

        p1 = models.Polynomial1D(4, n_models=2)
        y1 = p1(np.array([self.x1, self.x1]).T, model_set_axis=-1)
        assert y1.shape == (2, 20)
        assert_allclose(y1[0, :], y1[1, :], atol=10 ** (-12))

    def test_p2_1set_1pset(self):
        """1 pset, 1 2D data set, Polynomial2D"""

        p2 = models.Polynomial2D(5)
        z = p2(self.x, self.y)
        assert z.shape == (10, 10)

    def test_p2_nset_npset(self):
        """N param_sets, N 2D data sets, Poly2d"""

        p2 = models.Polynomial2D(5, n_models=2)
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
        p1 = models.Polynomial1D(5, n_models=2)
        yy = p1(xx, model_set_axis=-1)
        x1 = xx[:, 0]
        x2 = xx[:, 1]
        p1 = models.Polynomial1D(5)
        assert_allclose(p1(x1), yy[0, :], atol=10 ** (-12))
        p1 = models.Polynomial1D(5)
        assert_allclose(p1(x2), yy[1, :], atol=10 ** (-12))

    def test_evaluate_gauss2d(self):
        cov = np.array([[1., 0.8], [0.8, 3]])
        g = models.Gaussian2D(1., 5., 4., cov_matrix=cov)
        y, x = np.mgrid[:10, :10]
        g(x, y)


class TModel_1_1(Fittable1DModel):
    p1 = Parameter()
    p2 = Parameter()

    @staticmethod
    def evaluate(x, p1, p2):
        return x + p1 + p2


class TestSingleInputSingleOutputSingleModel(object):
    """
    A suite of tests to check various cases of parameter and input combinations
    on models with n_input = n_output = 1 on a toy model with n_models=1.

    Many of these tests mirror test cases in
    ``astropy.modeling.tests.test_parameters.TestParameterInitialization``,
    except that this tests how different parameter arrangements interact with
    different types of model inputs.
    """

    def test_scalar_parameters_scalar_input(self):
        """
        Scalar parameters with a scalar input should return a scalar.
        """

        t = TModel_1_1(1, 10)
        y = t(100)
        assert isinstance(y, float)
        assert np.ndim(y) == 0
        assert y == 111

    def test_scalar_parameters_1d_array_input(self):
        """
        Scalar parameters should broadcast with an array input to result in an
        array output of the same shape as the input.
        """

        t = TModel_1_1(1, 10)
        y = t(np.arange(5) * 100)
        assert isinstance(y, np.ndarray)
        assert np.shape(y) == (5,)
        assert np.all(y == [11, 111, 211, 311, 411])

    def test_scalar_parameters_2d_array_input(self):
        """
        Scalar parameters should broadcast with an array input to result in an
        array output of the same shape as the input.
        """

        t = TModel_1_1(1, 10)
        y = t(np.arange(6).reshape(2, 3) * 100)
        assert isinstance(y, np.ndarray)
        assert np.shape(y) == (2, 3)
        assert np.all(y == [[11, 111, 211],
                            [311, 411, 511]])

    def test_scalar_parameters_3d_array_input(self):
        """
        Scalar parameters should broadcast with an array input to result in an
        array output of the same shape as the input.
        """

        t = TModel_1_1(1, 10)
        y = t(np.arange(12).reshape(2, 3, 2) * 100)
        assert isinstance(y, np.ndarray)
        assert np.shape(y) == (2, 3, 2)
        assert np.all(y == [[[11, 111], [211, 311], [411, 511]],
                            [[611, 711], [811, 911], [1011, 1111]]])

    def test_1d_array_parameters_scalar_input(self):
        """
        Array parameters should all be broadcastable with each other, and with
        a scalar input the output should be broadcast to the maximum dimensions
        of the parameters.
        """

        t = TModel_1_1([1, 2], [10, 20])
        y = t(100)
        assert isinstance(y, np.ndarray)
        assert np.shape(y) == (2,)
        assert np.all(y == [111, 122])

    def test_1d_array_parameters_1d_array_input(self):
        """
        When given an array input it must be broadcastable with all the
        parameters.
        """

        t = TModel_1_1([1, 2], [10, 20])
        y1 = t([100, 200])
        assert np.shape(y1) == (2,)
        assert np.all(y1 == [111, 222])

        y2 = t([[100], [200]])
        assert np.shape(y2) == (2, 2)
        assert np.all(y2 == [[111, 122], [211, 222]])

        with pytest.raises(ValueError):
            # Doesn't broadcast
            y3 = t([100, 200, 300])

    def test_2d_array_parameters_2d_array_input(self):
        """
        When given an array input it must be broadcastable with all the
        parameters.
        """

        t = TModel_1_1([[1, 2], [3, 4]], [[10, 20], [30, 40]])

        y1 = t([[100, 200], [300, 400]])
        assert np.shape(y1) == (2, 2)
        assert np.all(y1 == [[111, 222], [333, 444]])

        y2 = t([[[[100]], [[200]]], [[[300]], [[400]]]])
        assert np.shape(y2) == (2, 2, 2, 2)
        assert np.all(y2 == [[[[111, 122], [133, 144]],
                              [[211, 222], [233, 244]]],
                             [[[311, 322], [333, 344]],
                              [[411, 422], [433, 444]]]])

        with pytest.raises(ValueError):
            # Doesn't broadcast
            y3 = t([[100, 200, 300], [400, 500, 600]])

    def test_mixed_array_parameters_1d_array_input(self):
        """
        When given an array input it must be broadcastable with all the
        parameters.
        """

        t = TModel_1_1([[[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]],
                           [[0.07, 0.08, 0.09], [0.10, 0.11, 0.12]]],
                          [1, 2, 3])

        y1 = t([10, 20, 30])
        assert np.shape(y1) == (2, 2, 3)
        assert_allclose(y1, [[[11.01, 22.02, 33.03], [11.04, 22.05, 33.06]],
                             [[11.07, 22.08, 33.09], [11.10, 22.11, 33.12]]])

        y2 = t([[[[10]]], [[[20]]], [[[30]]]])
        assert np.shape(y2) == (3, 2, 2, 3)
        assert_allclose(y2, [[[[11.01, 12.02, 13.03],
                               [11.04, 12.05, 13.06]],
                              [[11.07, 12.08, 13.09],
                               [11.10, 12.11, 13.12]]],
                             [[[21.01, 22.02, 23.03],
                               [21.04, 22.05, 23.06]],
                              [[21.07, 22.08, 23.09],
                               [21.10, 22.11, 23.12]]],
                             [[[31.01, 32.02, 33.03],
                               [31.04, 32.05, 33.06]],
                              [[31.07, 32.08, 33.09],
                               [31.10, 32.11, 33.12]]]])


class TestSingleInputSingleOutputTwoModel(object):
    """
    A suite of tests to check various cases of parameter and input combinations
    on models with n_input = n_output = 1 on a toy model with n_models=2.

    Many of these tests mirror test cases in
    ``astropy.modeling.tests.test_parameters.TestParameterInitialization``,
    except that this tests how different parameter arrangements interact with
    different types of model inputs.

    With n_models=2 all outputs should have a first dimension of size 2 (unless
    defined with model_set_axis != 0).
    """

    def test_scalar_parameters_scalar_input(self):
        """
        Scalar parameters with a scalar input should return a 1-D array with
        size equal to the number of models.
        """

        t = TModel_1_1([1, 2], [10, 20], n_models=2)

        y = t(100)
        assert np.shape(y) == (2,)
        assert np.all(y == [111, 122])

    def test_scalar_parameters_1d_array_input(self):
        """
        The dimension of the input should match the number of models unless
        model_set_axis=False is given, in which case the input is copied across
        all models.
        """

        t = TModel_1_1([1, 2], [10, 20], n_models=2)

        with pytest.raises(ValueError):
            y = t(np.arange(5) * 100)

        y1 = t([100, 200])
        assert np.shape(y1) == (2,)
        assert np.all(y1 == [111, 222])

        y2 = t([100, 200], model_set_axis=False)
        # In this case the value [100, 200, 300] should be evaluated on each
        # model rather than evaluating the first model with 100 and the second
        # model  with 200
        assert np.shape(y2) == (2, 2)
        assert np.all(y2 == [[111, 211], [122, 222]])

        y3 = t([100, 200, 300], model_set_axis=False)
        assert np.shape(y3) == (2, 3)
        assert np.all(y3 == [[111, 211, 311], [122, 222, 322]])

    def test_scalar_parameters_2d_array_input(self):
        """
        The dimension of the input should match the number of models unless
        model_set_axis=False is given, in which case the input is copied across
        all models.
        """

        t = TModel_1_1([1, 2], [10, 20], n_models=2)

        y1 = t(np.arange(6).reshape(2, 3) * 100)
        assert np.shape(y1) == (2, 3)
        assert np.all(y1 == [[11, 111, 211],
                            [322, 422, 522]])

        y2 = t(np.arange(6).reshape(2, 3) * 100, model_set_axis=False)
        assert np.shape(y2) == (2, 2, 3)
        assert np.all(y2 == [[[11, 111, 211], [311, 411, 511]],
                             [[22, 122, 222], [322, 422, 522]]])

    def test_scalar_parameters_3d_array_input(self):
        """
        The dimension of the input should match the number of models unless
        model_set_axis=False is given, in which case the input is copied across
        all models.
        """

        t = TModel_1_1([1, 2], [10, 20], n_models=2)
        data = np.arange(12).reshape(2, 3, 2) * 100

        y1 = t(data)
        assert np.shape(y1) == (2, 3, 2)
        assert np.all(y1 == [[[11, 111], [211, 311], [411, 511]],
                             [[622, 722], [822, 922], [1022, 1122]]])

        y2 = t(data, model_set_axis=False)
        assert np.shape(y2) == (2, 2, 3, 2)
        assert np.all(y2 == np.array([data + 11, data + 22]))

    def test_1d_array_parameters_scalar_input(self):
        """
        Array parameters should all be broadcastable with each other, and with
        a scalar input the output should be broadcast to the maximum dimensions
        of the parameters.
        """

        t = TModel_1_1([[1, 2, 3], [4, 5, 6]],
                          [[10, 20, 30], [40, 50, 60]], n_models=2)

        y = t(100)
        assert np.shape(y) == (2, 3)
        assert np.all(y == [[111, 122, 133], [144, 155, 166]])

    def test_1d_array_parameters_1d_array_input(self):
        """
        When the input is an array, if model_set_axis=False then it must
        broadcast with the shapes of the parameters (excluding the
        model_set_axis).

        Otherwise all dimensions must be broadcastable.
        """

        t = TModel_1_1([[1, 2, 3], [4, 5, 6]],
                          [[10, 20, 30], [40, 50, 60]], n_models=2)

        with pytest.raises(ValueError):
            y1 = t([100, 200, 300])

        y1 = t([100, 200])
        assert np.shape(y1) == (2, 3)
        assert np.all(y1 == [[111, 122, 133], [244, 255, 266]])

        with pytest.raises(ValueError):
            # Doesn't broadcast with the shape of the parameters, (3,)
            y2 = t([100, 200], model_set_axis=False)

        y2 = t([100, 200, 300], model_set_axis=False)
        assert np.shape(y2) == (2, 3)
        assert np.all(y2 == [[111, 222, 333],
                             [144, 255, 366]])

    def test_2d_array_parameters_2d_array_input(self):
        t = TModel_1_1([[[1, 2], [3, 4]], [[5, 6], [7, 8]]],
                          [[[10, 20], [30, 40]], [[50, 60], [70, 80]]],
                          n_models=2)
        y1 = t([[100, 200], [300, 400]])
        assert np.shape(y1) == (2, 2, 2)
        assert np.all(y1 == [[[111, 222], [133, 244]],
                             [[355, 466], [377, 488]]])

        with pytest.raises(ValueError):
            y2 = t([[100, 200, 300], [400, 500, 600]])

        y2 = t([[[100, 200], [300, 400]], [[500, 600], [700, 800]]])
        assert np.shape(y2) == (2, 2, 2)
        assert np.all(y2 == [[[111, 222], [333, 444]],
                             [[555, 666], [777, 888]]])

    def test_mixed_array_parameters_1d_array_input(self):
        t = TModel_1_1([[[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]],
                           [[0.07, 0.08, 0.09], [0.10, 0.11, 0.12]]],
                          [[1, 2, 3], [4, 5, 6]], n_models=2)

        with pytest.raises(ValueError):
            y = t([10, 20, 30])

        y = t([10, 20, 30], model_set_axis=False)
        assert np.shape(y) == (2, 2, 3)
        assert_allclose(y, [[[11.01, 22.02, 33.03], [11.04, 22.05, 33.06]],
                            [[14.07, 25.08, 36.09], [14.10, 25.11, 36.12]]])


class TModel_1_2(FittableModel):
    inputs = ('x',)
    outputs = ('y', 'z')

    p1 = Parameter()
    p2 = Parameter()
    p3 = Parameter()

    @staticmethod
    def evaluate(x, p1, p2, p3):
        return (x + p1 + p2, x + p1 + p2 + p3)


class TestSingleInputDoubleOutputSingleModel(object):
    """
    A suite of tests to check various cases of parameter and input combinations
    on models with n_input = 1 but n_output = 2 on a toy model with n_models=1.

    As of writing there are not enough controls to adjust how outputs from such
    a model should be formatted (currently the shapes of outputs are assumed to
    be directly associated with the shapes of corresponding inputs when
    n_inputs == n_outputs).  For now, the approach taken for cases like this is
    to assume all outputs should have the same format.
    """

    def test_scalar_parameters_scalar_input(self):
        """
        Scalar parameters with a scalar input should return a scalar.
        """

        t = TModel_1_2(1, 10, 1000)
        y, z = t(100)
        assert isinstance(y, float)
        assert isinstance(z, float)
        assert np.ndim(y) == np.ndim(z) == 0
        assert y == 111
        assert z == 1111

    def test_scalar_parameters_1d_array_input(self):
        """
        Scalar parameters should broadcast with an array input to result in an
        array output of the same shape as the input.
        """

        t = TModel_1_2(1, 10, 1000)
        y, z = t(np.arange(5) * 100)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)
        assert np.shape(y) == np.shape(z) == (5,)
        assert np.all(y == [11, 111, 211, 311, 411])
        assert np.all(z == (y + 1000))

    def test_scalar_parameters_2d_array_input(self):
        """
        Scalar parameters should broadcast with an array input to result in an
        array output of the same shape as the input.
        """

        t = TModel_1_2(1, 10, 1000)
        y, z = t(np.arange(6).reshape(2, 3) * 100)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)
        assert np.shape(y) == np.shape(z) == (2, 3)
        assert np.all(y == [[11, 111, 211],
                            [311, 411, 511]])
        assert np.all(z == (y + 1000))

    def test_scalar_parameters_3d_array_input(self):
        """
        Scalar parameters should broadcast with an array input to result in an
        array output of the same shape as the input.
        """

        t = TModel_1_2(1, 10, 1000)
        y, z = t(np.arange(12).reshape(2, 3, 2) * 100)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)
        assert np.shape(y) == np.shape(z) == (2, 3, 2)
        assert np.all(y == [[[11, 111], [211, 311], [411, 511]],
                            [[611, 711], [811, 911], [1011, 1111]]])
        assert np.all(z == (y + 1000))

    def test_1d_array_parameters_scalar_input(self):
        """
        Array parameters should all be broadcastable with each other, and with
        a scalar input the output should be broadcast to the maximum dimensions
        of the parameters.
        """

        t = TModel_1_2([1, 2], [10, 20], [1000, 2000])
        y, z = t(100)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)
        assert np.shape(y) == np.shape(z) == (2,)
        assert np.all(y == [111, 122])
        assert np.all(z == [1111, 2122])

    def test_1d_array_parameters_1d_array_input(self):
        """
        When given an array input it must be broadcastable with all the
        parameters.
        """

        t = TModel_1_2([1, 2], [10, 20], [1000, 2000])
        y1, z1 = t([100, 200])
        assert np.shape(y1) == np.shape(z1) == (2,)
        assert np.all(y1 == [111, 222])
        assert np.all(z1 == [1111, 2222])

        y2, z2 = t([[100], [200]])
        assert np.shape(y2) == np.shape(z2) == (2, 2)
        assert np.all(y2 == [[111, 122], [211, 222]])
        assert np.all(z2 == [[1111, 2122], [1211, 2222]])

        with pytest.raises(ValueError):
            # Doesn't broadcast
            y3, z3 = t([100, 200, 300])

    def test_2d_array_parameters_2d_array_input(self):
        """
        When given an array input it must be broadcastable with all the
        parameters.
        """

        t = TModel_1_2([[1, 2], [3, 4]], [[10, 20], [30, 40]],
                          [[1000, 2000], [3000, 4000]])

        y1, z1 = t([[100, 200], [300, 400]])
        assert np.shape(y1) == np.shape(z1) == (2, 2)
        assert np.all(y1 == [[111, 222], [333, 444]])
        assert np.all(z1 == [[1111, 2222], [3333, 4444]])

        y2, z2 = t([[[[100]], [[200]]], [[[300]], [[400]]]])
        assert np.shape(y2) == np.shape(z2) == (2, 2, 2, 2)
        assert np.all(y2 == [[[[111, 122], [133, 144]],
                              [[211, 222], [233, 244]]],
                             [[[311, 322], [333, 344]],
                              [[411, 422], [433, 444]]]])
        assert np.all(z2 == [[[[1111, 2122], [3133, 4144]],
                              [[1211, 2222], [3233, 4244]]],
                             [[[1311, 2322], [3333, 4344]],
                              [[1411, 2422], [3433, 4444]]]])

        with pytest.raises(ValueError):
            # Doesn't broadcast
            y3, z3 = t([[100, 200, 300], [400, 500, 600]])

    def test_mixed_array_parameters_1d_array_input(self):
        """
        When given an array input it must be broadcastable with all the
        parameters.
        """

        t = TModel_1_2([[[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]],
                           [[0.07, 0.08, 0.09], [0.10, 0.11, 0.12]]],
                          [1, 2, 3], [100, 200, 300])

        y1, z1 = t([10, 20, 30])
        assert np.shape(y1) == np.shape(z1) == (2, 2, 3)
        assert_allclose(y1, [[[11.01, 22.02, 33.03], [11.04, 22.05, 33.06]],
                             [[11.07, 22.08, 33.09], [11.10, 22.11, 33.12]]])
        assert_allclose(z1, [[[111.01, 222.02, 333.03],
                              [111.04, 222.05, 333.06]],
                             [[111.07, 222.08, 333.09],
                              [111.10, 222.11, 333.12]]])


        y2, z2 = t([[[[10]]], [[[20]]], [[[30]]]])
        assert np.shape(y2) == np.shape(z2) == (3, 2, 2, 3)
        assert_allclose(y2, [[[[11.01, 12.02, 13.03],
                               [11.04, 12.05, 13.06]],
                              [[11.07, 12.08, 13.09],
                               [11.10, 12.11, 13.12]]],
                             [[[21.01, 22.02, 23.03],
                               [21.04, 22.05, 23.06]],
                              [[21.07, 22.08, 23.09],
                               [21.10, 22.11, 23.12]]],
                             [[[31.01, 32.02, 33.03],
                               [31.04, 32.05, 33.06]],
                              [[31.07, 32.08, 33.09],
                               [31.10, 32.11, 33.12]]]])

        assert_allclose(z2, [[[[111.01, 212.02, 313.03],
                               [111.04, 212.05, 313.06]],
                              [[111.07, 212.08, 313.09],
                               [111.10, 212.11, 313.12]]],
                             [[[121.01, 222.02, 323.03],
                               [121.04, 222.05, 323.06]],
                              [[121.07, 222.08, 323.09],
                               [121.10, 222.11, 323.12]]],
                             [[[131.01, 232.02, 333.03],
                               [131.04, 232.05, 333.06]],
                              [[131.07, 232.08, 333.09],
                               [131.10, 232.11, 333.12]]]])


class TInputFormatter(Model):
    """
    A toy model to test input/output formatting.
    """

    inputs = ('x', 'y')
    outputs = ('x', 'y')

    @staticmethod
    def evaluate(x, y):
        return x, y


def test_format_input_scalars():
    model = TInputFormatter()
    result = model(1, 2)
    assert result == (1, 2)


def test_format_input_arrays():
    model = TInputFormatter()
    result = model([1, 1], [2, 2])
    assert_allclose(result, (np.array([1, 1]), np.array([2, 2])))


def test_format_input_arrays_transposed():
    model = TInputFormatter()
    input = np.array([[1, 1]]).T, np.array([[2, 2]]).T
    result = model(*input)
    assert_allclose(result, input)

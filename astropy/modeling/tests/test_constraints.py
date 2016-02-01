# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import types

from ..core import Fittable1DModel
from ..parameters import Parameter
from .. import models
from .. import fitting
import numpy as np
from numpy.testing import utils
from numpy.random import RandomState
from ...tests.helper import pytest
from .utils import ignore_non_integer_warning

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class TestNonLinearConstraints(object):

    def setup_class(self):
        self.g1 = models.Gaussian1D(10, 14.9, stddev=.3)
        self.g2 = models.Gaussian1D(10, 13, stddev=.4)
        self.x = np.arange(10, 20, .1)
        self.y1 = self.g1(self.x)
        self.y2 = self.g2(self.x)
        rsn = RandomState(1234567890)
        self.n = rsn.randn(100)
        self.ny1 = self.y1 + 2 * self.n
        self.ny2 = self.y2 + 2 * self.n

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_fixed_par(self):
        g1 = models.Gaussian1D(10, mean=14.9, stddev=.3,
                               fixed={'amplitude': True})
        fitter = fitting.LevMarLSQFitter()
        model = fitter(g1, self.x, self.ny1)
        assert model.amplitude.value == 10

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_tied_par(self):

        def tied(model):
            mean = 50 * model.stddev
            return mean
        g1 = models.Gaussian1D(10, mean=14.9, stddev=.3, tied={'mean': tied})
        fitter = fitting.LevMarLSQFitter()
        model = fitter(g1, self.x, self.ny1)
        utils.assert_allclose(model.mean.value, 50 * model.stddev,
                              rtol=10 ** (-5))

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_joint_fitter(self):
        g1 = models.Gaussian1D(10, 14.9, stddev=.3)
        g2 = models.Gaussian1D(10, 13, stddev=.4)
        jf = fitting.JointFitter([g1, g2], {g1: ['amplitude'],
                                            g2: ['amplitude']}, [9.8])
        x = np.arange(10, 20, .1)
        y1 = g1(x)
        y2 = g2(x)
        n = np.random.randn(100)
        ny1 = y1 + 2 * n
        ny2 = y2 + 2 * n
        jf(x, ny1, x, ny2)
        p1 = [14.9, .3]
        p2 = [13, .4]
        A = 9.8
        p = np.r_[A, p1, p2]

        def compmodel(A, p, x):
            return A * np.exp(-0.5 / p[1] ** 2 * (x - p[0]) ** 2)

        def errf(p, x1, y1, x2, y2):
            return np.ravel(
                np.r_[compmodel(p[0], p[1:3], x1) - y1,
                      compmodel(p[0], p[3:], x2) - y2])

        fitparams, _ = optimize.leastsq(errf, p, args=(x, ny1, x, ny2))
        utils.assert_allclose(jf.fitparams, fitparams, rtol=10 ** (-5))
        utils.assert_allclose(g1.amplitude.value, g2.amplitude.value)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_no_constraints(self):
        g1 = models.Gaussian1D(9.9, 14.5, stddev=.3)

        def func(p, x):
            return p[0] * np.exp(-0.5 / p[2] ** 2 * (x - p[1]) ** 2)

        def errf(p, x, y):
            return func(p, x) - y

        p0 = [9.9, 14.5, 0.3]
        y = g1(self.x)
        n = np.random.randn(100)
        ny = y + n
        fitpar, s = optimize.leastsq(errf, p0, args=(self.x, ny))
        fitter = fitting.LevMarLSQFitter()
        model = fitter(g1, self.x, ny)
        utils.assert_allclose(model.parameters, fitpar, rtol=5 * 10 ** (-3))


@pytest.mark.skipif('not HAS_SCIPY')
class TestBounds(object):

    def setup_class(self):
        A = -2.0
        B = 0.5
        self.x = np.linspace(-1.0, 1.0, 100)
        self.y = A * self.x + B + np.random.normal(scale=0.1, size=100)
        data = np.array([505.0, 556.0, 630.0, 595.0, 561.0, 553.0, 543.0, 496.0, 460.0, 469.0,
                         426.0, 518.0, 684.0, 798.0, 830.0, 794.0, 649.0, 706.0, 671.0, 545.0,
                         479.0, 454.0, 505.0, 700.0, 1058.0, 1231.0, 1325.0, 997.0, 1036.0, 884.0,
                         610.0, 487.0, 453.0, 527.0, 780.0, 1094.0, 1983.0, 1993.0, 1809.0, 1525.0,
                         1056.0, 895.0, 604.0, 466.0, 510.0, 678.0, 1130.0, 1986.0, 2670.0, 2535.0,
                         1878.0, 1450.0, 1200.0, 663.0, 511.0, 474.0, 569.0, 848.0, 1670.0, 2611.0,
                         3129.0, 2507.0, 1782.0, 1211.0, 723.0, 541.0, 511.0, 518.0, 597.0, 1137.0,
                         1993.0, 2925.0, 2438.0, 1910.0, 1230.0, 738.0, 506.0, 461.0, 486.0, 597.0,
                         733.0, 1262.0, 1896.0, 2342.0, 1792.0, 1180.0, 667.0, 482.0, 454.0, 482.0,
                         504.0, 566.0, 789.0, 1194.0, 1545.0, 1361.0, 933.0, 562.0, 418.0, 463.0,
                         435.0, 466.0, 528.0, 487.0, 664.0, 799.0, 746.0, 550.0, 478.0, 535.0, 443.0,
                         416.0, 439.0, 472.0, 472.0, 492.0, 523.0, 569.0, 487.0, 441.0, 428.0])
        self.data = data.reshape(11, 11)

    def test_bounds_lsq(self):
        guess_slope = 1.1
        guess_intercept = 0.0
        bounds = {'slope': (-1.5, 5.0), 'intercept': (-1.0, 1.0)}
        line_model = models.Linear1D(guess_slope, guess_intercept,
                                     bounds=bounds)
        fitter = fitting.LevMarLSQFitter()
        model = fitter(line_model, self.x, self.y)
        slope = model.slope.value
        intercept = model.intercept.value
        assert slope + 10 ** -5 >= bounds['slope'][0]
        assert slope - 10 ** -5 <= bounds['slope'][1]
        assert intercept + 10 ** -5 >= bounds['intercept'][0]
        assert intercept - 10 ** -5 <= bounds['intercept'][1]

    def test_bounds_slsqp(self):
        guess_slope = 1.1
        guess_intercept = 0.0
        bounds = {'slope': (-1.5, 5.0), 'intercept': (-1.0, 1.0)}
        line_model = models.Linear1D(guess_slope, guess_intercept,
                                     bounds=bounds)
        fitter = fitting.SLSQPLSQFitter()

        with ignore_non_integer_warning():
            model = fitter(line_model, self.x, self.y)

        slope = model.slope.value
        intercept = model.intercept.value
        assert slope + 10 ** -5 >= bounds['slope'][0]
        assert slope - 10 ** -5 <= bounds['slope'][1]
        assert intercept + 10 ** -5 >= bounds['intercept'][0]
        assert intercept - 10 ** -5 <= bounds['intercept'][1]

    def test_bounds_gauss2d_lsq(self):
        X, Y = np.meshgrid(np.arange(11), np.arange(11))
        bounds = {"x_mean": [0., 11.],
                  "y_mean": [0., 11.],
                  "x_stddev": [1., 4],
                  "y_stddev": [1., 4]}
        gauss = models.Gaussian2D(amplitude=10., x_mean=5., y_mean=5.,
                                  x_stddev=4., y_stddev=4., theta=0.5,
                                  bounds=bounds)
        gauss_fit = fitting.LevMarLSQFitter()
        model = gauss_fit(gauss, X, Y, self.data)
        x_mean = model.x_mean.value
        y_mean = model.y_mean.value
        x_stddev = model.x_stddev.value
        y_stddev = model.y_stddev.value
        assert x_mean + 10 ** -5 >= bounds['x_mean'][0]
        assert x_mean - 10 ** -5 <= bounds['x_mean'][1]
        assert y_mean + 10 ** -5 >= bounds['y_mean'][0]
        assert y_mean - 10 ** -5 <= bounds['y_mean'][1]
        assert x_stddev + 10 ** -5 >= bounds['x_stddev'][0]
        assert x_stddev - 10 ** -5 <= bounds['x_stddev'][1]
        assert y_stddev + 10 ** -5 >= bounds['y_stddev'][0]
        assert y_stddev - 10 ** -5 <= bounds['y_stddev'][1]

    def test_bounds_gauss2d_slsqp(self):
        X, Y = np.meshgrid(np.arange(11), np.arange(11))
        bounds = {"x_mean": [0., 11.],
                  "y_mean": [0., 11.],
                  "x_stddev": [1., 4],
                  "y_stddev": [1., 4]}
        gauss = models.Gaussian2D(amplitude=10., x_mean=5., y_mean=5.,
                                  x_stddev=4., y_stddev=4., theta=0.5,
                                  bounds=bounds)
        gauss_fit = fitting.SLSQPLSQFitter()
        with ignore_non_integer_warning():
            model = gauss_fit(gauss, X, Y, self.data)
        x_mean = model.x_mean.value
        y_mean = model.y_mean.value
        x_stddev = model.x_stddev.value
        y_stddev = model.y_stddev.value
        assert x_mean + 10 ** -5 >= bounds['x_mean'][0]
        assert x_mean - 10 ** -5 <= bounds['x_mean'][1]
        assert y_mean + 10 ** -5 >= bounds['y_mean'][0]
        assert y_mean - 10 ** -5 <= bounds['y_mean'][1]
        assert x_stddev + 10 ** -5 >= bounds['x_stddev'][0]
        assert x_stddev - 10 ** -5 <= bounds['x_stddev'][1]
        assert y_stddev + 10 ** -5 >= bounds['y_stddev'][0]
        assert y_stddev - 10 ** -5 <= bounds['y_stddev'][1]


class TestLinearConstraints(object):

    def setup_class(self):
        self.p1 = models.Polynomial1D(4)
        self.p1.c0 = 0
        self.p1.c1 = 0
        self.p1.window = [0., 9.]
        self.x = np.arange(10)
        self.y = self.p1(self.x)
        rsn = RandomState(1234567890)
        self.n = rsn.randn(10)
        self.ny = self.y + self.n

    def test(self):
        self.p1.c0.fixed = True
        self.p1.c1.fixed = True
        pfit = fitting.LinearLSQFitter()
        model = pfit(self.p1, self.x, self.y)
        utils.assert_allclose(self.y, model(self.x))

# Test constraints as parameter properties


def test_set_fixed_1():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1)
    gauss.mean.fixed = True
    assert gauss.fixed == {'amplitude': False, 'mean': True, 'stddev': False}


def test_set_fixed_2():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1,
                              fixed={'mean': True})
    assert gauss.mean.fixed is True


def test_set_tied_1():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1)
    gauss.amplitude.tied = tie_amplitude
    assert gauss.amplitude.tied is not False
    assert isinstance(gauss.tied['amplitude'], types.FunctionType)


def test_set_tied_2():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1,
                              tied={'amplitude': tie_amplitude})
    assert gauss.amplitude.tied != False


def test_unset_fixed():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1,
                              fixed={'mean': True})
    gauss.mean.fixed = False
    assert gauss.fixed == {'amplitude': False, 'mean': False, 'stddev': False}


def test_unset_tied():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1,
                              tied={'amplitude': tie_amplitude})
    gauss.amplitude.tied = False
    assert gauss.tied == {'amplitude': False, 'mean': False, 'stddev': False}


def test_set_bounds_1():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1,
                              bounds={'stddev': (0, None)})
    assert gauss.bounds == {'amplitude': (None, None),
                            'mean': (None, None),
                            'stddev': (0.0, None)}


def test_set_bounds_2():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1)
    gauss.stddev.min = 0.
    assert gauss.bounds == {'amplitude': (None, None),
                            'mean': (None, None),
                            'stddev': (0.0, None)}


def test_unset_bounds():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1,
                              bounds={'stddev': (0, 2)})
    gauss.stddev.min = None
    gauss.stddev.max = None
    assert gauss.bounds == {'amplitude': (None, None),
                            'mean': (None, None),
                            'stddev': (None, None)}


def test_default_constraints():
    """Regression test for https://github.com/astropy/astropy/issues/2396

    Ensure that default constraints defined on parameters are carried through
    to instances of the models those parameters are defined for.
    """

    class MyModel(Fittable1DModel):
        a = Parameter(default=1)
        b = Parameter(default=0, min=0, fixed=True)

        @staticmethod
        def evaluate(x, a, b):
            return x * a + b

    assert MyModel.a.default == 1
    assert MyModel.b.default == 0
    assert MyModel.b.min == 0
    assert MyModel.b.bounds == (0, None)
    assert MyModel.b.fixed is True

    m = MyModel()
    assert m.a.value == 1
    assert m.b.value == 0
    assert m.b.min == 0
    assert m.b.bounds == (0, None)
    assert m.b.fixed is True
    assert m.bounds == {'a': (None, None), 'b': (0, None)}
    assert m.fixed == {'a': False, 'b': True}

    # Make a model instance that overrides the default constraints and values
    m = MyModel(3, 4, bounds={'a': (1, None), 'b': (2, None)},
                fixed={'a': True, 'b': False})
    assert m.a.value == 3
    assert m.b.value == 4
    assert m.a.min == 1
    assert m.b.min == 2
    assert m.a.bounds == (1, None)
    assert m.b.bounds == (2, None)
    assert m.a.fixed is True
    assert m.b.fixed is False
    assert m.bounds == {'a': (1, None), 'b': (2, None)}
    assert m.fixed == {'a': True, 'b': False}


@pytest.mark.skipif('not HAS_SCIPY')
def test_fit_with_fixed_and_bound_constraints():
    """
    Regression test for https://github.com/astropy/astropy/issues/2235

    Currently doesn't test that the fit is any *good*--just that parameters
    stay within their given constraints.
    """

    m = models.Gaussian1D(amplitude=3, mean=4, stddev=1,
                          bounds={'mean': (4, 5)},
                          fixed={'amplitude': True})
    x = np.linspace(0, 10, 10)
    y = np.exp(-x ** 2 / 2)

    f = fitting.LevMarLSQFitter()
    fitted_1 = f(m, x, y)
    assert fitted_1.mean >= 4
    assert fitted_1.mean <= 5
    assert fitted_1.amplitude == 3.0

    m.amplitude.fixed = False
    fitted_2 = f(m, x, y)
    # It doesn't matter anymore what the amplitude ends up as so long as the
    # bounds constraint was still obeyed
    assert fitted_1.mean >= 4
    assert fitted_1.mean <= 5


@pytest.mark.skipif('not HAS_SCIPY')
def test_fit_with_bound_constraints_estimate_jacobian():
    """
    Regression test for https://github.com/astropy/astropy/issues/2400

    Checks that bounds constraints are obeyed on a custom model that does not
    define fit_deriv (and thus its Jacobian must be estimated for non-linear
    fitting).
    """

    class MyModel(Fittable1DModel):
        a = Parameter(default=1)
        b = Parameter(default=2)

        @staticmethod
        def evaluate(x, a, b):
            return a * x + b

    m_real = MyModel(a=1.5, b=-3)
    x = np.arange(100)
    y = m_real(x)

    m = MyModel()
    f = fitting.LevMarLSQFitter()
    fitted_1 = f(m, x, y)

    # This fit should be trivial so even without constraints on the bounds it
    # should be right
    assert np.allclose(fitted_1.a, 1.5)
    assert np.allclose(fitted_1.b, -3)

    m2 = MyModel()
    m2.a.bounds = (-2, 2)
    f2 = fitting.LevMarLSQFitter()
    fitted_2 = f2(m2, x, y)
    assert np.allclose(fitted_1.a, 1.5)
    assert np.allclose(fitted_1.b, -3)

    # Check that the estimated Jacobian was computed (it doesn't matter what
    # the values are so long as they're not all zero.
    assert np.any(f2.fit_info['fjac'] != 0)

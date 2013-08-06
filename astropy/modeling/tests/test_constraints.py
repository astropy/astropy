# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import types
from .. import models
from .. import fitting
import numpy as np
from numpy.testing import utils
from numpy.random import RandomState
from ...tests.helper import pytest

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class TestNonLinearConstraints(object):

    def setup_class(self):
        self.g1 = models.Gaussian1DModel(10, 14.9, stddev=.3)
        self.g2 = models.Gaussian1DModel(10, 13, stddev=.4)
        self.x = np.arange(10, 20, .1)
        self.y1 = self.g1(self.x)
        self.y2 = self.g2(self.x)
        rsn = RandomState(1234567890)
        self.n = rsn.randn(100)
        self.ny1 = self.y1 + 2 * self.n
        self.ny2 = self.y2 + 2 * self.n

    @pytest.mark.skipif('not HAS_SCIPY')
    def testFixedPar(self):
        g1 = models.Gaussian1DModel(10, mean=14.9, stddev=.3, fixed={'amplitude': True})
        fitter = fitting.NonLinearLSQFitter(g1)
        fitter(self.x, self.ny1)
        assert g1.amplitude.value == 10

    @pytest.mark.skipif('not HAS_SCIPY')
    def testTiedPar(self):

        def tied(model):
            mean = 50 * model.stddev
            return mean
        g1 = models.Gaussian1DModel(10, mean=14.9, stddev=.3, tied={'mean': tied})
        fitter = fitting.NonLinearLSQFitter(g1)
        fitter(self.x, self.ny1)
        utils.assert_allclose(g1.mean.value, 50 * g1.stddev, rtol=10 ** (-5))

    @pytest.mark.skipif('not HAS_SCIPY')
    def testJointFitter(self):
        g1 = models.Gaussian1DModel(10, 14.9, stddev=.3)
        g2 = models.Gaussian1DModel(10, 13, stddev=.4)
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
        compmodel = lambda A, p, x: A * np.exp(-0.5 / p[1] ** 2 * (x - p[0]) ** 2)
        errf = lambda p, x1, y1, x2, y2: np.ravel(
            np.r_[compmodel(p[0], p[1:3], x1) - y1,
                  compmodel(p[0], p[3:], x2) - y2])
        fitpars, _ = optimize.leastsq(errf, p, args=(x, ny1, x, ny2))
        utils.assert_allclose(jf.fitpars, fitpars, rtol=10 ** (-5))
        utils.assert_allclose(g1.amplitude.value, g2.amplitude.value)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_no_constraints(self):
        g1 = models.Gaussian1DModel(9.9, 14.5, stddev=.3)
        func = lambda p, x: p[0] * np.exp(-0.5 / p[2] ** 2 * (x - p[1]) ** 2)
        errf = lambda p, x, y: func(p, x) - y
        p0 = [9.9, 14.5, 0.3]
        y = g1(self.x)
        n = np.random.randn(100)
        ny = y + n
        fitpar, s = optimize.leastsq(errf, p0, args=(self.x, n + y))
        fitter = fitting.NonLinearLSQFitter(g1)
        fitter(self.x, n + y)
        utils.assert_allclose(g1.parameters, fitpar, rtol=5 * 10 ** (-3))


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
        line_model = models.Linear1DModel(guess_slope, guess_intercept, bounds=bounds)
        fitter = fitting.NonLinearLSQFitter(line_model)
        fitter(self.x, self.y)
        slope = line_model.slope.value
        intercept = line_model.intercept.value
        assert slope + 10 ** -5 >= bounds['slope'][0]
        assert slope - 10 ** -5 <= bounds['slope'][1]
        assert intercept + 10 ** -5 >= bounds['intercept'][0]
        assert intercept - 10 ** -5 <= bounds['intercept'][1]

    def test_bounds_slsqp(self):
        guess_slope = 1.1
        guess_intercept = 0.0
        bounds = {'slope': (-1.5, 5.0), 'intercept': (-1.0, 1.0)}
        line_model = models.Linear1DModel(guess_slope, guess_intercept, bounds=bounds)
        fitter = fitting.SLSQPFitter(line_model)
        fitter(self.x, self.y)
        slope = line_model.slope.value
        intercept = line_model.intercept.value
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
        gauss = models.Gaussian2DModel(amplitude=10., x_mean=5., y_mean=5.,
                                       x_stddev=4., y_stddev=4., theta=0.5,
                                       bounds=bounds)
        gauss_fit = fitting.NonLinearLSQFitter(gauss)
        gauss_fit(X, Y, self.data)
        x_mean = gauss.x_mean.value
        y_mean = gauss.y_mean.value
        x_stddev = gauss.x_stddev.value
        y_stddev = gauss.y_stddev.value
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
        gauss = models.Gaussian2DModel(amplitude=10., x_mean=5., y_mean=5.,
                                       x_stddev=4., y_stddev=4., theta=0.5,
                                       bounds=bounds)
        gauss_fit = fitting.SLSQPFitter(gauss)
        gauss_fit(X, Y, self.data)
        x_mean = gauss.x_mean.value
        y_mean = gauss.y_mean.value
        x_stddev = gauss.x_stddev.value
        y_stddev = gauss.y_stddev.value
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
        self.p1 = models.Poly1DModel(4)
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
        pfit = fitting.LinearLSQFitter(self.p1)
        pfit(self.x, self.y)
        utils.assert_allclose(self.y, self.p1(self.x))

# Test constraints as parameter properties


def test_set_fixed_1():
    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1)
    gauss.mean.fixed = True
    assert gauss.fixed == {'amplitude': False, 'mean': True, 'stddev': False}
    gfit = fitting.NonLinearLSQFitter(gauss)
    assert gfit.fixed == [False, True, False]


def test_set_fixed_2():
    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1, fixed={'mean': True})
    assert gauss.mean.fixed is True
    assert gauss.fixed == {'amplitude': False, 'mean': True, 'stddev': False}


def test_set_tied_1():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1)
    gauss.amplitude.tied = tie_amplitude
    assert gauss.amplitude.tied is not False
    assert type(gauss.tied['amplitude']) == types.FunctionType


def test_set_tied_2():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1,
                                   tied={'amplitude': tie_amplitude})
    assert gauss.amplitude.tied is not False


def test_unset_fixed():
    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1, fixed={'mean': True})
    gfit = fitting.NonLinearLSQFitter(gauss)
    gauss.mean.fixed = False
    assert gauss.fixed == {'amplitude': False, 'mean': False, 'stddev': False}
    gfit._update_constraints()
    assert gfit.fixed == [False, False, False]


def test_unset_tied():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1,
                                   tied={'amplitude': tie_amplitude})
    gfit = fitting.NonLinearLSQFitter(gauss)
    gauss.amplitude.tied = False
    assert gauss.tied == {'amplitude': False, 'mean': False, 'stddev': False}
    gfit._update_constraints()
    assert gfit.fixed == [False, False, False]


def test_set_bounds_1():
    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1,
                                   bounds={'stddev': (0, None)})
    gfit = fitting.NonLinearLSQFitter(gauss)
    assert gauss.bounds == {'amplitude': (None, None),
                            'mean': (None, None),
                            'stddev': (0.0, None)}
    assert gfit.bounds == [(-10 ** 12, 10 ** 12), (-10 ** 12, 10 ** 12), (0.0, 10 ** 12)]


def test_set_bounds_2():
    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1)
    gfit = fitting.NonLinearLSQFitter(gauss)
    gauss.stddev.min = 0.
    assert gauss.bounds == {'amplitude': (None, None),
                            'mean': (None, None),
                            'stddev': (0.0, None)}
    gfit._update_constraints()
    assert gfit.bounds == [(-10 ** 12, 10 ** 12), (-10 ** 12, 10 ** 12), (0.0, 10 ** 12)]


def test_unset_bounds():
    gauss = models.Gaussian1DModel(amplitude=20, mean=2, stddev=1,
                                   bounds={'stddev': (0, 2)})
    gauss.stddev.min = None
    gauss.stddev.max = None
    assert gauss.bounds == {'amplitude': (None, None),
                            'mean': (None, None),
                            'stddev': (None, None)}
    gfit = fitting.NonLinearLSQFitter(gauss)
    assert gfit.bounds == [(-10 ** 12, 10 ** 12), (-10 ** 12, 10 ** 12), (-10 ** 12, 10 ** 12)]

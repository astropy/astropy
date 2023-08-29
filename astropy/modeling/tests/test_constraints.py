# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name

import platform
import types
import warnings
from contextlib import nullcontext

import numpy as np
import pytest
from numpy.random import default_rng
from numpy.testing import assert_allclose

from astropy.modeling import fitting, models
from astropy.modeling.core import Fittable1DModel
from astropy.modeling.parameters import Parameter
from astropy.utils import minversion
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning

SCIPY_LT_1_11_2 = not minversion("scipy", "1.11.2") if HAS_SCIPY else True

fitters = [
    fitting.LevMarLSQFitter,
    fitting.TRFLSQFitter,
    fitting.LevMarLSQFitter,
    fitting.DogBoxLSQFitter,
]


class TestNonLinearConstraints:
    def setup_class(self):
        self.g1 = models.Gaussian1D(10, 14.9, stddev=0.3)
        self.g2 = models.Gaussian1D(10, 13, stddev=0.4)
        self.x = np.arange(10, 20, 0.1)
        self.y1 = self.g1(self.x)
        self.y2 = self.g2(self.x)
        rsn = default_rng(1234567890)
        self.n = rsn.standard_normal(100)
        self.ny1 = self.y1 + 2 * self.n
        self.ny2 = self.y2 + 2 * self.n

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    @pytest.mark.parametrize("fitter", fitters)
    def test_fixed_par(self, fitter):
        fitter = fitter()

        g1 = models.Gaussian1D(10, mean=14.9, stddev=0.3, fixed={"amplitude": True})
        model = fitter(g1, self.x, self.ny1)
        assert model.amplitude.value == 10

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    @pytest.mark.parametrize("fitter", fitters)
    def test_tied_par(self, fitter):
        fitter = fitter()

        def tied(model):
            mean = 50 * model.stddev
            return mean

        g1 = models.Gaussian1D(10, mean=14.9, stddev=0.3, tied={"mean": tied})
        model = fitter(g1, self.x, self.ny1)
        assert_allclose(model.mean.value, 50 * model.stddev, rtol=10 ** (-5))

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    def test_joint_fitter(self):
        from scipy import optimize

        g1 = models.Gaussian1D(10, 14.9, stddev=0.3)
        g2 = models.Gaussian1D(10, 13, stddev=0.4)
        jf = fitting.JointFitter(
            [g1, g2], {g1: ["amplitude"], g2: ["amplitude"]}, [9.8]
        )
        x = np.arange(10, 20, 0.1)
        y1 = g1(x)
        y2 = g2(x)
        n = np.random.randn(100)
        ny1 = y1 + 2 * n
        ny2 = y2 + 2 * n
        jf(x, ny1, x, ny2)
        p1 = [14.9, 0.3]
        p2 = [13, 0.4]
        A = 9.8
        p = np.r_[A, p1, p2]

        def compmodel(A, p, x):
            return A * np.exp(-0.5 / p[1] ** 2 * (x - p[0]) ** 2)

        def errf(p, x1, y1, x2, y2):
            return np.ravel(
                np.r_[compmodel(p[0], p[1:3], x1) - y1, compmodel(p[0], p[3:], x2) - y2]
            )

        fitparams, _ = optimize.leastsq(errf, p, args=(x, ny1, x, ny2))
        assert_allclose(jf.fitparams, fitparams, rtol=10 ** (-5))
        assert_allclose(g1.amplitude.value, g2.amplitude.value)

    @pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
    @pytest.mark.parametrize("fitter", fitters)
    def test_no_constraints(self, fitter):
        from scipy import optimize

        fitter = fitter()

        g1 = models.Gaussian1D(9.9, 14.5, stddev=0.3)

        def func(p, x):
            return p[0] * np.exp(-0.5 / p[2] ** 2 * (x - p[1]) ** 2)

        def errf(p, x, y):
            return func(p, x) - y

        p0 = [9.9, 14.5, 0.3]
        y = g1(self.x)
        n = np.random.randn(100)
        ny = y + n
        fitpar, s = optimize.leastsq(errf, p0, args=(self.x, ny))
        model = fitter(g1, self.x, ny)
        assert_allclose(model.parameters, fitpar, rtol=5 * 10 ** (-3))


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
class TestBounds:
    def setup_class(self):
        A = -2.0
        B = 0.5
        self.x = np.linspace(-1.0, 1.0, 100)
        self.y = A * self.x + B + np.random.normal(scale=0.1, size=100)
        # fmt: off
        data = np.array(
            [
                505.0,  556.0,  630.0,  595.0,  561.0,  553.0,  543.0,  496.0,  460.0,  469.0,
                426.0,  518.0,  684.0,  798.0,  830.0,  794.0,  649.0,  706.0,  671.0,  545.0,
                479.0,  454.0,  505.0,  700.0,  1058.0, 1231.0, 1325.0, 997.0,  1036.0, 884.0,
                610.0,  487.0,  453.0,  527.0,  780.0,  1094.0, 1983.0, 1993.0, 1809.0, 1525.0,
                1056.0, 895.0,  604.0,  466.0,  510.0,  678.0,  1130.0, 1986.0, 2670.0, 2535.0,
                1878.0, 1450.0, 1200.0, 663.0,  511.0,  474.0,  569.0,  848.0,  1670.0, 2611.0,
                3129.0, 2507.0, 1782.0, 1211.0, 723.0,  541.0,  511.0,  518.0,  597.0,  1137.0,
                1993.0, 2925.0, 2438.0, 1910.0, 1230.0, 738.0,  506.0,  461.0,  486.0,  597.0,
                733.0,  1262.0, 1896.0, 2342.0, 1792.0, 1180.0, 667.0,  482.0,  454.0,  482.0,
                504.0,  566.0,  789.0,  1194.0, 1545.0, 1361.0, 933.0,  562.0,  418.0,  463.0,
                435.0,  466.0,  528.0,  487.0,  664.0,  799.0,  746.0,  550.0,  478.0,  535.0,
                443.0,  416.0,  439.0,  472.0,  472.0,  492.0,  523.0,  569.0,  487.0,  441.0,
                428.0
            ]
        )
        # fmt: on

        self.data = data.reshape(11, 11)

    @pytest.mark.parametrize("fitter", fitters)
    def test_bounds_lsq(self, fitter):
        fitter = fitter()

        guess_slope = 1.1
        guess_intercept = 0.0
        bounds = {"slope": (-1.5, 5.0), "intercept": (-1.0, 1.0)}
        line_model = models.Linear1D(guess_slope, guess_intercept, bounds=bounds)
        with pytest.warns(AstropyUserWarning, match=r"Model is linear in parameters"):
            model = fitter(line_model, self.x, self.y)
        slope = model.slope.value
        intercept = model.intercept.value
        assert slope + 10**-5 >= bounds["slope"][0]
        assert slope - 10**-5 <= bounds["slope"][1]
        assert intercept + 10**-5 >= bounds["intercept"][0]
        assert intercept - 10**-5 <= bounds["intercept"][1]

    def test_bounds_slsqp(self):
        guess_slope = 1.1
        guess_intercept = 0.0
        bounds = {"slope": (-1.5, 5.0), "intercept": (-1.0, 1.0)}
        line_model = models.Linear1D(guess_slope, guess_intercept, bounds=bounds)
        fitter = fitting.SLSQPLSQFitter()
        with pytest.warns(
            AstropyUserWarning, match="consider using linear fitting methods"
        ):
            model = fitter(line_model, self.x, self.y)

        slope = model.slope.value
        intercept = model.intercept.value
        assert slope + 10**-5 >= bounds["slope"][0]
        assert slope - 10**-5 <= bounds["slope"][1]
        assert intercept + 10**-5 >= bounds["intercept"][0]
        assert intercept - 10**-5 <= bounds["intercept"][1]

    @pytest.mark.filterwarnings("ignore:The fit may be unsuccessful")
    @pytest.mark.parametrize("fitter", fitters)
    def test_bounds_gauss2d_lsq(self, fitter):
        fitter = fitter()

        X, Y = np.meshgrid(np.arange(11), np.arange(11))
        bounds = {
            "x_mean": [0.0, 11.0],
            "y_mean": [0.0, 11.0],
            "x_stddev": [1.0, 4],
            "y_stddev": [1.0, 4],
        }
        gauss = models.Gaussian2D(
            amplitude=10.0,
            x_mean=5.0,
            y_mean=5.0,
            x_stddev=4.0,
            y_stddev=4.0,
            theta=0.5,
            bounds=bounds,
        )
        if isinstance(fitter, fitting.TRFLSQFitter):
            ctx = np.errstate(invalid="ignore", divide="ignore")
        else:
            ctx = nullcontext()
        with ctx:
            model = fitter(gauss, X, Y, self.data)
        x_mean = model.x_mean.value
        y_mean = model.y_mean.value
        x_stddev = model.x_stddev.value
        y_stddev = model.y_stddev.value
        assert x_mean + 10**-5 >= bounds["x_mean"][0]
        assert x_mean - 10**-5 <= bounds["x_mean"][1]
        assert y_mean + 10**-5 >= bounds["y_mean"][0]
        assert y_mean - 10**-5 <= bounds["y_mean"][1]
        assert x_stddev + 10**-5 >= bounds["x_stddev"][0]
        assert x_stddev - 10**-5 <= bounds["x_stddev"][1]
        assert y_stddev + 10**-5 >= bounds["y_stddev"][0]
        assert y_stddev - 10**-5 <= bounds["y_stddev"][1]

    def test_bounds_gauss2d_slsqp(self):
        X, Y = np.meshgrid(np.arange(11), np.arange(11))
        bounds = {
            "x_mean": [0.0, 11.0],
            "y_mean": [0.0, 11.0],
            "x_stddev": [1.0, 4],
            "y_stddev": [1.0, 4],
        }
        gauss = models.Gaussian2D(
            amplitude=10.0,
            x_mean=5.0,
            y_mean=5.0,
            x_stddev=4.0,
            y_stddev=4.0,
            theta=0.5,
            bounds=bounds,
        )
        gauss_fit = fitting.SLSQPLSQFitter()
        # Warning does not appear in all the CI jobs.
        # TODO: Rewrite the test for more consistent warning behavior.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r".*The fit may be unsuccessful.*",
                category=AstropyUserWarning,
            )
            model = gauss_fit(gauss, X, Y, self.data)
        x_mean = model.x_mean.value
        y_mean = model.y_mean.value
        x_stddev = model.x_stddev.value
        y_stddev = model.y_stddev.value
        assert x_mean + 10**-5 >= bounds["x_mean"][0]
        assert x_mean - 10**-5 <= bounds["x_mean"][1]
        assert y_mean + 10**-5 >= bounds["y_mean"][0]
        assert y_mean - 10**-5 <= bounds["y_mean"][1]
        assert x_stddev + 10**-5 >= bounds["x_stddev"][0]
        assert x_stddev - 10**-5 <= bounds["x_stddev"][1]
        assert y_stddev + 10**-5 >= bounds["y_stddev"][0]
        assert y_stddev - 10**-5 <= bounds["y_stddev"][1]


class TestLinearConstraints:
    def setup_class(self):
        self.p1 = models.Polynomial1D(4)
        self.p1.c0 = 0
        self.p1.c1 = 0
        self.p1.window = [0.0, 9.0]
        self.x = np.arange(10)
        self.y = self.p1(self.x)
        rsn = default_rng(1234567890)
        self.n = rsn.standard_normal(10)
        self.ny = self.y + self.n

    def test(self):
        self.p1.c0.fixed = True
        self.p1.c1.fixed = True
        pfit = fitting.LinearLSQFitter()
        model = pfit(self.p1, self.x, self.y)
        assert_allclose(self.y, model(self.x))


# Test constraints as parameter properties


def test_set_fixed_1():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1)
    gauss.mean.fixed = True
    assert gauss.fixed == {"amplitude": False, "mean": True, "stddev": False}


def test_set_fixed_2():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1, fixed={"mean": True})
    assert gauss.mean.fixed is True


def test_set_tied_1():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1)
    gauss.amplitude.tied = tie_amplitude
    assert gauss.amplitude.tied is not False
    assert isinstance(gauss.tied["amplitude"], types.FunctionType)


def test_set_tied_2():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1D(
        amplitude=20, mean=2, stddev=1, tied={"amplitude": tie_amplitude}
    )
    assert gauss.amplitude.tied


def test_unset_fixed():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1, fixed={"mean": True})
    gauss.mean.fixed = False
    assert gauss.fixed == {"amplitude": False, "mean": False, "stddev": False}


def test_unset_tied():
    def tie_amplitude(model):
        return 50 * model.stddev

    gauss = models.Gaussian1D(
        amplitude=20, mean=2, stddev=1, tied={"amplitude": tie_amplitude}
    )
    gauss.amplitude.tied = False
    assert gauss.tied == {"amplitude": False, "mean": False, "stddev": False}


def test_set_bounds_1():
    gauss = models.Gaussian1D(
        amplitude=20, mean=2, stddev=1, bounds={"stddev": (0, None)}
    )
    assert gauss.bounds == {
        "amplitude": (None, None),
        "mean": (None, None),
        "stddev": (0.0, None),
    }


def test_set_bounds_2():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1)
    gauss.stddev.min = 0.0
    assert gauss.bounds == {
        "amplitude": (None, None),
        "mean": (None, None),
        "stddev": (0.0, None),
    }


def test_unset_bounds():
    gauss = models.Gaussian1D(amplitude=20, mean=2, stddev=1, bounds={"stddev": (0, 2)})
    gauss.stddev.min = None
    gauss.stddev.max = None
    assert gauss.bounds == {
        "amplitude": (None, None),
        "mean": (None, None),
        "stddev": (None, None),
    }


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
    assert m.bounds == {"a": (None, None), "b": (0, None)}
    assert m.fixed == {"a": False, "b": True}

    # Make a model instance that overrides the default constraints and values
    m = MyModel(
        3, 4, bounds={"a": (1, None), "b": (2, None)}, fixed={"a": True, "b": False}
    )
    assert m.a.value == 3
    assert m.b.value == 4
    assert m.a.min == 1
    assert m.b.min == 2
    assert m.a.bounds == (1, None)
    assert m.b.bounds == (2, None)
    assert m.a.fixed is True
    assert m.b.fixed is False
    assert m.bounds == {"a": (1, None), "b": (2, None)}
    assert m.fixed == {"a": True, "b": False}


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.filterwarnings(r"ignore:divide by zero encountered.*")
@pytest.mark.parametrize("fitter", fitters)
def test_fit_with_fixed_and_bound_constraints(fitter):
    """
    Regression test for https://github.com/astropy/astropy/issues/2235

    Currently doesn't test that the fit is any *good*--just that parameters
    stay within their given constraints.
    """

    # DogBoxLSQFitter causes failure on s390x, aremel possibly others (not x86_64 or arm64)
    if fitter == fitting.DogBoxLSQFitter and (
        platform.machine() not in ("x86_64", "arm64")
    ):
        pytest.xfail(
            "DogBoxLSQFitter can to be unstable on non-standard platforms leading to "
            "random test failures"
        )

    fitter = fitter()

    m = models.Gaussian1D(
        amplitude=3,
        mean=4,
        stddev=1,
        bounds={"mean": (4, 5)},
        fixed={"amplitude": True},
    )
    x = np.linspace(0, 10, 10)
    y = np.exp(-(x**2) / 2)

    if isinstance(fitter, fitting.TRFLSQFitter):
        ctx = np.errstate(invalid="ignore", divide="ignore")
    else:
        ctx = nullcontext()

    with ctx:
        fitted_1 = fitter(m, x, y)
        assert fitted_1.mean >= 4
        assert fitted_1.mean <= 5
        assert fitted_1.amplitude == 3.0

        m.amplitude.fixed = False
        # Cannot enter np.errstate twice, so we need to indent everything in between.
        _ = fitter(m, x, y)
    # It doesn't matter anymore what the amplitude ends up as so long as the
    # bounds constraint was still obeyed
    assert fitted_1.mean >= 4
    assert fitted_1.mean <= 5


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("fitter", fitters)
def test_fit_with_bound_constraints_estimate_jacobian(fitter):
    """
    Regression test for https://github.com/astropy/astropy/issues/2400

    Checks that bounds constraints are obeyed on a custom model that does not
    define fit_deriv (and thus its Jacobian must be estimated for non-linear
    fitting).
    """
    fitter = fitter()

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
    fitted_1 = fitter(m, x, y)

    # This fit should be trivial so even without constraints on the bounds it
    # should be right
    assert np.allclose(fitted_1.a, 1.5)
    assert np.allclose(fitted_1.b, -3)

    m2 = MyModel()
    m2.a.bounds = (-2, 2)
    _ = fitter(m2, x, y)
    assert np.allclose(fitted_1.a, 1.5)
    assert np.allclose(fitted_1.b, -3)

    # Check that the estimated Jacobian was computed (it doesn't matter what
    # the values are so long as they're not all zero.
    if fitter == fitting.LevMarLSQFitter:
        assert np.any(fitter.fit_info["fjac"] != 0)


# https://github.com/astropy/astropy/issues/6014
@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("fitter", fitters)
def test_gaussian2d_positive_stddev(fitter):
    # This is 2D Gaussian with noise to be fitted, as provided by @ysBach
    fitter = fitter()

    # fmt: off
    test = [
        [-54.33, 13.81, -34.55, 8.95, -143.71, -0.81, 59.25, -14.78, -204.9,
         -30.87, -124.39, 123.53, 70.81, -109.48, -106.77, 35.64, 18.29],
        [-126.19, -89.13, 63.13, 50.74, 61.83, 19.06, 65.7, 77.94, 117.14,
         139.37, 52.57, 236.04, 100.56, 242.28, -180.62, 154.02, -8.03],
        [91.43, 96.45, -118.59, -174.58, -116.49, 80.11, -86.81, 14.62, 79.26,
         7.56, 54.99, 260.13, -136.42, -20.77, -77.55, 174.52, 134.41],
        [33.88, 7.63, 43.54, 70.99, 69.87, 33.97, 273.75, 176.66, 201.94,
         336.34, 340.54, 163.77, -156.22, 21.49, -148.41, 94.88, 42.55],
        [82.28, 177.67, 26.81, 17.66, 47.81, -31.18, 353.23, 589.11, 553.27,
         242.35, 444.12, 186.02, 140.73, 75.2, -87.98, -18.23, 166.74],
        [113.09, -37.01, 134.23, 71.89, 107.88, 198.69, 273.88, 626.63, 551.8,
         547.61, 580.35, 337.8, 139.8, 157.64, -1.67, -26.99, 37.35],
        [106.47, 31.97, 84.99, -125.79, 195.0, 493.65, 861.89, 908.31, 803.9,
         781.01, 532.59, 404.67, 115.18, 111.11, 28.08, 122.05, -58.36],
        [183.62, 45.22, 40.89, 111.58, 425.81, 321.53, 545.09, 866.02, 784.78,
         731.35, 609.01, 405.41, -19.65, 71.2, -140.5, 144.07, 25.24],
        [137.13, -86.95, 15.39, 180.14, 353.23, 699.01, 1033.8, 1014.49,
         814.11, 647.68, 461.03, 249.76, 94.8, 41.17, -1.16, 183.76, 188.19],
        [35.39, 26.92, 198.53, -37.78, 638.93, 624.41, 816.04, 867.28, 697.0,
         491.56, 378.21, -18.46, -65.76, 98.1, 12.41, -102.18, 119.05],
        [190.73, 125.82, 311.45, 369.34, 554.39, 454.37, 755.7, 736.61, 542.43,
         188.24, 214.86, 217.91, 7.91, 27.46, -172.14, -82.36, -80.31],
        [-55.39, 80.18, 267.19, 274.2, 169.53, 327.04, 488.15, 437.53, 225.38,
         220.94, 4.01, -92.07, 39.68, 57.22, 144.66, 100.06, 34.96],
        [130.47, -4.23, 46.3, 101.49, 115.01, 217.38, 249.83, 115.9, 87.36,
         105.81, -47.86, -9.94, -82.28, 144.45, 83.44, 23.49, 183.9],
        [-110.38, -115.98, 245.46, 103.51, 255.43, 163.47, 56.52, 33.82,
         -33.26, -111.29, 88.08, 193.2, -100.68, 15.44, 86.32, -26.44, -194.1],
        [109.36, 96.01, -124.89, -16.4, 84.37, 114.87, -65.65, -58.52, -23.22,
         42.61, 144.91, -209.84, 110.29, 66.37, -117.85, -147.73, -122.51],
        [10.94, 45.98, 118.12, -46.53, -72.14, -74.22, 21.22, 0.39, 86.03,
         23.97, -45.42, 12.05, -168.61, 27.79, 61.81, 84.07, 28.79],
        [46.61, -104.11, 56.71, -90.85, -16.51, -66.45, -141.34, 0.96, 58.08,
         285.29, -61.41, -9.01, -323.38, 58.35, 80.14, -101.22, 145.65]
    ]
    # fmt: on

    g_init = models.Gaussian2D(x_mean=8, y_mean=8)
    if isinstance(fitter, (fitting.TRFLSQFitter, fitting.DogBoxLSQFitter)):
        pytest.xfail("TRFLSQFitter seems to be broken for this test.")

    y, x = np.mgrid[:17, :17]
    g_fit = fitter(g_init, x, y, test)

    # Compare with @ysBach original result:
    # - x_stddev was negative, so its abs value is used for comparison here.
    # - theta is beyond (-90, 90) deg, which doesn't make sense, so ignored.
    assert_allclose(
        [g_fit.amplitude.value, g_fit.y_stddev.value],
        [984.7694929790363, 3.1840618351417307],
        rtol=1.5e-6,
    )
    assert_allclose(g_fit.x_mean.value, 7.198391516587464)
    assert_allclose(g_fit.y_mean.value, 7.49720660088511, rtol=5e-7)
    assert_allclose(g_fit.x_stddev.value, 1.9840185107597297, rtol=2e-6)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.filterwarnings(r"ignore:Model is linear in parameters.*")
@pytest.mark.parametrize("fitter", fitters)
def test_2d_model(fitter):
    """Issue #6403"""
    from astropy.utils import NumpyRNGContext

    fitter = fitter()

    # 2D model with LevMarLSQFitter
    gauss2d = models.Gaussian2D(10.2, 4.3, 5, 2, 1.2, 1.4)
    X = np.linspace(-1, 7, 200)
    Y = np.linspace(-1, 7, 200)
    x, y = np.meshgrid(X, Y)
    z = gauss2d(x, y)
    w = np.ones(x.size)
    w.shape = x.shape

    with NumpyRNGContext(1234567890):
        n = np.random.randn(x.size)
        n.shape = x.shape
        m = fitter(gauss2d, x, y, z + 2 * n, weights=w)
        assert_allclose(m.parameters, gauss2d.parameters, rtol=0.05)
        m = fitter(gauss2d, x, y, z + 2 * n, weights=None)
        assert_allclose(m.parameters, gauss2d.parameters, rtol=0.05)

        # 2D model with LevMarLSQFitter, fixed constraint
        gauss2d.x_stddev.fixed = True
        m = fitter(gauss2d, x, y, z + 2 * n, weights=w)
        assert_allclose(m.parameters, gauss2d.parameters, rtol=0.05)
        m = fitter(gauss2d, x, y, z + 2 * n, weights=None)
        assert_allclose(m.parameters, gauss2d.parameters, rtol=0.05)

        # Polynomial2D, col_fit_deriv=False
        p2 = models.Polynomial2D(1, c0_0=1, c1_0=1.2, c0_1=3.2)
        z = p2(x, y)
        m = fitter(p2, x, y, z + 2 * n, weights=None)
        assert_allclose(m.parameters, p2.parameters, rtol=0.05)
        m = fitter(p2, x, y, z + 2 * n, weights=w)
        assert_allclose(m.parameters, p2.parameters, rtol=0.05)

        # Polynomial2D, col_fit_deriv=False, fixed constraint
        p2.c1_0.fixed = True
        m = fitter(p2, x, y, z + 2 * n, weights=w)
        assert_allclose(m.parameters, p2.parameters, rtol=0.05)
        m = fitter(p2, x, y, z + 2 * n, weights=None)
        assert_allclose(m.parameters, p2.parameters, rtol=0.05)


def test_set_prior_posterior():
    model = models.Polynomial1D(1)
    model.c0.prior = models.Gaussian1D(2.3, 2, 0.1)
    assert model.c0.prior(2) == 2.3

    model.c0.posterior = models.Linear1D(1, 0.2)
    assert model.c0.posterior(1) == 1.2


def test_set_constraints():
    g = models.Gaussian1D()
    p = models.Polynomial1D(1)

    # Set bounds before model combination
    g.stddev.bounds = (0, 3)
    m = g + p
    assert m.bounds == {
        "amplitude_0": (None, None),
        "mean_0": (None, None),
        "stddev_0": (0.0, 3.0),
        "c0_1": (None, None),
        "c1_1": (None, None),
    }

    # Set bounds on the compound model
    m.stddev_0.bounds = (1, 3)
    assert m.bounds == {
        "amplitude_0": (None, None),
        "mean_0": (None, None),
        "stddev_0": (1.0, 3.0),
        "c0_1": (None, None),
        "c1_1": (None, None),
    }

    # Set the bounds of a Parameter directly in the bounds dict
    m.bounds["stddev_0"] = (4, 5)
    assert m.bounds == {
        "amplitude_0": (None, None),
        "mean_0": (None, None),
        "stddev_0": (4, 5),
        "c0_1": (None, None),
        "c1_1": (None, None),
    }

    # Set the bounds of a Parameter on the child model bounds dict
    g.bounds["stddev"] = (1, 5)
    m = g + p
    assert m.bounds == {
        "amplitude_0": (None, None),
        "mean_0": (None, None),
        "stddev_0": (1, 5),
        "c0_1": (None, None),
        "c1_1": (None, None),
    }

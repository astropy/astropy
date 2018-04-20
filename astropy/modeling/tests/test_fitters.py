# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Module to test fitting routines
"""

import os.path

import pytest
import numpy as np
from numpy import linalg
from numpy.testing import assert_allclose, assert_almost_equal
from unittest import mock

from . import irafutil
from .. import models
from ..core import Fittable2DModel, Parameter
from ..fitting import *
from ...utils import NumpyRNGContext
from ...utils.data import get_pkg_data_filename
from .utils import ignore_non_integer_warning
from ...stats import sigma_clip

from ...utils.exceptions import AstropyUserWarning
from ..fitting import populate_entry_points
import warnings

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    from pkg_resources import EntryPoint
    HAS_PKG = True
except ImportError:
    HAS_PKG = False


fitters = [SimplexLSQFitter, SLSQPLSQFitter]

_RANDOM_SEED = 0x1337


class TestPolynomial2D:
    """Tests for 2D polynomail fitting."""

    def setup_class(self):
        self.model = models.Polynomial2D(2)
        self.y, self.x = np.mgrid[:5, :5]

        def poly2(x, y):
            return 1 + 2 * x + 3 * x ** 2 + 4 * y + 5 * y ** 2 + 6 * x * y
        self.z = poly2(self.x, self.y)
        self.fitter = LinearLSQFitter()

    def test_poly2D_fitting(self):
        v = self.model.fit_deriv(x=self.x, y=self.y)
        p = linalg.lstsq(v, self.z.flatten())[0]
        new_model = self.fitter(self.model, self.x, self.y, self.z)
        assert_allclose(new_model.parameters, p)

    def test_eval(self):
        new_model = self.fitter(self.model, self.x, self.y, self.z)
        assert_allclose(new_model(self.x, self.y), self.z)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_polynomial2D_nonlinear_fitting(self):
        self.model.parameters = [.6, 1.8, 2.9, 3.7, 4.9, 6.7]
        nlfitter = LevMarLSQFitter()
        new_model = nlfitter(self.model, self.x, self.y, self.z)
        assert_allclose(new_model.parameters, [1, 2, 3, 4, 5, 6])


class TestICheb2D:
    """
    Tests 2D Chebyshev polynomial fitting

    Create a 2D polynomial (z) using Polynomial2DModel and default coefficients
    Fit z using a ICheb2D model
    Evaluate the ICheb2D polynomial and compare with the initial z
    """

    def setup_class(self):
        self.pmodel = models.Polynomial2D(2)
        self.y, self.x = np.mgrid[:5, :5]
        self.z = self.pmodel(self.x, self.y)
        self.cheb2 = models.Chebyshev2D(2, 2)
        self.fitter = LinearLSQFitter()

    def test_default_params(self):
        self.cheb2.parameters = np.arange(9)
        p = np.array([1344., 1772., 400., 1860., 2448., 552., 432., 568.,
                      128.])
        z = self.cheb2(self.x, self.y)
        model = self.fitter(self.cheb2, self.x, self.y, z)
        assert_almost_equal(model.parameters, p)

    def test_poly2D_cheb2D(self):
        model = self.fitter(self.cheb2, self.x, self.y, self.z)
        z1 = model(self.x, self.y)
        assert_almost_equal(self.z, z1)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_chebyshev2D_nonlinear_fitting(self):
        cheb2d = models.Chebyshev2D(2, 2)
        cheb2d.parameters = np.arange(9)
        z = cheb2d(self.x, self.y)
        cheb2d.parameters = [0.1, .6, 1.8, 2.9, 3.7, 4.9, 6.7, 7.5, 8.9]
        nlfitter = LevMarLSQFitter()
        model = nlfitter(cheb2d, self.x, self.y, z)
        assert_allclose(model.parameters, [0, 1, 2, 3, 4, 5, 6, 7, 8],
                        atol=10**-9)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_chebyshev2D_nonlinear_fitting_with_weights(self):
        cheb2d = models.Chebyshev2D(2, 2)
        cheb2d.parameters = np.arange(9)
        z = cheb2d(self.x, self.y)
        cheb2d.parameters = [0.1, .6, 1.8, 2.9, 3.7, 4.9, 6.7, 7.5, 8.9]
        nlfitter = LevMarLSQFitter()
        weights = np.ones_like(self.y)
        model = nlfitter(cheb2d, self.x, self.y, z, weights=weights)
        assert_allclose(model.parameters, [0, 1, 2, 3, 4, 5, 6, 7, 8],
                        atol=10**-9)


@pytest.mark.skipif('not HAS_SCIPY')
class TestJointFitter:

    """
    Tests the joint fitting routine using 2 gaussian models
    """

    def setup_class(self):
        """
        Create 2 gaussian models and some data with noise.
        Create a fitter for the two models keeping the amplitude parameter
        common for the two models.
        """
        self.g1 = models.Gaussian1D(10, mean=14.9, stddev=.3)
        self.g2 = models.Gaussian1D(10, mean=13, stddev=.4)
        self.jf = JointFitter([self.g1, self.g2],
                                      {self.g1: ['amplitude'],
                                       self.g2: ['amplitude']}, [9.8])
        self.x = np.arange(10, 20, .1)
        y1 = self.g1(self.x)
        y2 = self.g2(self.x)

        with NumpyRNGContext(_RANDOM_SEED):
            n = np.random.randn(100)

        self.ny1 = y1 + 2 * n
        self.ny2 = y2 + 2 * n
        self.jf(self.x, self.ny1, self.x, self.ny2)

    def test_joint_parameter(self):
        """
        Tests that the amplitude of the two models is the same
        """
        assert_allclose(self.jf.fitparams[0], self.g1.parameters[0])
        assert_allclose(self.jf.fitparams[0], self.g2.parameters[0])

    def test_joint_fitter(self):
        """
        Tests the fitting routine with similar procedure.
        Compares the fitted parameters.
        """
        p1 = [14.9, .3]
        p2 = [13, .4]
        A = 9.8
        p = np.r_[A, p1, p2]

        def model(A, p, x):
            return A * np.exp(-0.5 / p[1] ** 2 * (x - p[0]) ** 2)

        def errfunc(p, x1, y1, x2, y2):
            return np.ravel(np.r_[model(p[0], p[1:3], x1) - y1,
                            model(p[0], p[3:], x2) - y2])

        coeff, _ = optimize.leastsq(errfunc, p,
                                    args=(self.x, self.ny1, self.x, self.ny2))
        assert_allclose(coeff, self.jf.fitparams, rtol=10 ** (-2))

class TestLinearLSQFitter:
    def test_compound_model_raises_error(self):
        """Test that if an user tries to use a compound model, raises an error"""
        with pytest.raises(ValueError) as excinfo:
            init_model1 = models.Polynomial1D(degree=2, c0=[1, 1], n_models=2)
            init_model2 = models.Polynomial1D(degree=2, c0=[1, 1], n_models=2)
            init_model_comp = init_model1 + init_model2
            x = np.arange(10)
            y = init_model_comp(x, model_set_axis=False)
            fitter = LinearLSQFitter()
            fitted_model = fitter(init_model_comp, x, y)
        assert "Model must be simple, not compound" in str(excinfo.value)

    def test_chebyshev1D(self):
        """Tests fitting a 1D Chebyshev polynomial to some real world data."""

        test_file = get_pkg_data_filename(os.path.join('data',
                                                       'idcompspec.fits'))
        with open(test_file) as f:
            lines = f.read()
            reclist = lines.split('begin')

        record = irafutil.IdentifyRecord(reclist[1])
        coeffs = record.coeff
        order = int(record.fields['order'])

        initial_model = models.Chebyshev1D(order - 1,
                                           domain=record.get_range())
        fitter = LinearLSQFitter()

        fitted_model = fitter(initial_model, record.x, record.z)
        assert_allclose(fitted_model.parameters, np.array(coeffs),
                        rtol=10e-2)

    def test_linear_fit_model_set(self):
        """Tests fitting multiple models simultaneously."""

        init_model = models.Polynomial1D(degree=2, c0=[1, 1], n_models=2)
        x = np.arange(10)
        y_expected = init_model(x, model_set_axis=False)
        assert y_expected.shape == (2, 10)

        # Add a bit of random noise
        with NumpyRNGContext(_RANDOM_SEED):
            y = y_expected + np.random.normal(0, 0.01, size=y_expected.shape)

        fitter = LinearLSQFitter()
        fitted_model = fitter(init_model, x, y)
        assert_allclose(fitted_model(x, model_set_axis=False), y_expected,
                        rtol=1e-1)

    def test_linear_fit_2d_model_set(self):
        """Tests fitted multiple 2-D models simultaneously."""

        init_model = models.Polynomial2D(degree=2, c0_0=[1, 1], n_models=2)
        x = np.arange(10)
        y = np.arange(10)
        z_expected = init_model(x, y, model_set_axis=False)
        assert z_expected.shape == (2, 10)

        # Add a bit of random noise
        with NumpyRNGContext(_RANDOM_SEED):
            z = z_expected + np.random.normal(0, 0.01, size=z_expected.shape)

        fitter = LinearLSQFitter()
        fitted_model = fitter(init_model, x, y, z)
        assert_allclose(fitted_model(x, y, model_set_axis=False), z_expected,
                        rtol=1e-1)

    def test_linear_fit_fixed_parameter(self):
        """
        Tests fitting a polynomial model with a fixed parameter (issue #6135).
        """
        init_model = models.Polynomial1D(degree=2, c1=1)
        init_model.c1.fixed = True

        x = np.arange(10)
        y = 2 + x + 0.5*x*x

        fitter = LinearLSQFitter()
        fitted_model = fitter(init_model, x, y)
        assert_allclose(fitted_model.parameters, [2., 1., 0.5], atol=1e-14)

    def test_linear_fit_model_set_fixed_parameter(self):
        """
        Tests fitting a polynomial model set with a fixed parameter (#6135).
        """
        init_model = models.Polynomial1D(degree=2, c1=[1, -2], n_models=2)
        init_model.c1.fixed = True

        x = np.arange(10)
        yy = np.array([2 + x + 0.5*x*x, -2*x])

        fitter = LinearLSQFitter()
        fitted_model = fitter(init_model, x, yy)

        assert_allclose(fitted_model.c0, [2., 0.], atol=1e-14)
        assert_allclose(fitted_model.c1, [1., -2.], atol=1e-14)
        assert_allclose(fitted_model.c2, [0.5, 0.], atol=1e-14)

    def test_linear_fit_2d_model_set_fixed_parameters(self):
        """
        Tests fitting a 2d polynomial model set with fixed parameters (#6135).
        """
        init_model = models.Polynomial2D(degree=2, c1_0=[1, 2], c0_1=[-0.5, 1],
                                         n_models=2,
                                         fixed={'c1_0': True, 'c0_1': True})

        x, y = np.mgrid[0:5, 0:5]
        zz = np.array([1+x-0.5*y+0.1*x*x, 2*x+y-0.2*y*y])

        fitter = LinearLSQFitter()
        fitted_model = fitter(init_model, x, y, zz)

        assert_allclose(fitted_model(x, y, model_set_axis=False), zz,
                        atol=1e-14)

    def test_linear_fit_model_set_masked_values(self):
        """
        Tests model set fitting with masked value(s) (#4824, #6819).
        """
        # NB. For single models, there is an equivalent doctest.

        init_model = models.Polynomial1D(degree=1, n_models=2)
        x = np.arange(10)
        y = np.ma.masked_array([2*x+1, x-2], mask=np.zeros_like([x, x]))

        y[0, 7] = 100.  # throw off fit coefficients if unmasked
        y.mask[0, 7] = True
        y[1, 1:3] = -100.
        y.mask[1, 1:3] = True

        fitter = LinearLSQFitter()
        fitted_model = fitter(init_model, x, y)

        assert_allclose(fitted_model.c0, [1., -2.], atol=1e-14)
        assert_allclose(fitted_model.c1, [2., 1.], atol=1e-14)

    def test_linear_fit_2d_model_set_masked_values(self):
        """
        Tests 2D model set fitting with masked value(s) (#4824, #6819).
        """
        init_model = models.Polynomial2D(1, n_models=2)
        x, y = np.mgrid[0:5, 0:5]
        z = np.ma.masked_array([2*x+3*y+1, x-0.5*y-2],
                               mask=np.zeros_like([x, x]))

        z[0, 3, 1] = -1000.  # throw off fit coefficients if unmasked
        z.mask[0, 3, 1] = True

        fitter = LinearLSQFitter()
        fitted_model = fitter(init_model, x, y, z)

        assert_allclose(fitted_model.c0_0, [1., -2.], atol=1e-14)
        assert_allclose(fitted_model.c1_0, [2., 1.], atol=1e-14)
        assert_allclose(fitted_model.c0_1, [3., -0.5], atol=1e-14)


@pytest.mark.skipif('not HAS_SCIPY')
class TestNonLinearFitters:
    """Tests non-linear least squares fitting and the SLSQP algorithm."""

    def setup_class(self):
        self.initial_values = [100, 5, 1]

        self.xdata = np.arange(0, 10, 0.1)
        sigma = 4. * np.ones_like(self.xdata)

        with NumpyRNGContext(_RANDOM_SEED):
            yerror = np.random.normal(0, sigma)

        def func(p, x):
            return p[0] * np.exp(-0.5 / p[2] ** 2 * (x - p[1]) ** 2)

        self.ydata = func(self.initial_values, self.xdata) + yerror
        self.gauss = models.Gaussian1D(100, 5, stddev=1)

    def test_estimated_vs_analytic_deriv(self):
        """
        Runs `LevMarLSQFitter` with estimated and analytic derivatives of a
        `Gaussian1D`.
        """

        fitter = LevMarLSQFitter()
        model = fitter(self.gauss, self.xdata, self.ydata)
        g1e = models.Gaussian1D(100, 5.0, stddev=1)
        efitter = LevMarLSQFitter()
        emodel = efitter(g1e, self.xdata, self.ydata, estimate_jacobian=True)
        assert_allclose(model.parameters, emodel.parameters, rtol=10 ** (-3))

    def test_estimated_vs_analytic_deriv_with_weights(self):
        """
        Runs `LevMarLSQFitter` with estimated and analytic derivatives of a
        `Gaussian1D`.
        """

        weights = 1.0 / (self.ydata / 10.)

        fitter = LevMarLSQFitter()
        model = fitter(self.gauss, self.xdata, self.ydata, weights=weights)
        g1e = models.Gaussian1D(100, 5.0, stddev=1)
        efitter = LevMarLSQFitter()
        emodel = efitter(g1e, self.xdata, self.ydata, weights=weights, estimate_jacobian=True)
        assert_allclose(model.parameters, emodel.parameters, rtol=10 ** (-3))

    def test_with_optimize(self):
        """
        Tests results from `LevMarLSQFitter` against `scipy.optimize.leastsq`.
        """

        fitter = LevMarLSQFitter()
        model = fitter(self.gauss, self.xdata, self.ydata,
                       estimate_jacobian=True)

        def func(p, x):
            return p[0] * np.exp(-0.5 / p[2] ** 2 * (x - p[1]) ** 2)

        def errfunc(p, x, y):
            return func(p, x) - y

        result = optimize.leastsq(errfunc, self.initial_values,
                                  args=(self.xdata, self.ydata))
        assert_allclose(model.parameters, result[0], rtol=10 ** (-3))

    def test_with_weights(self):
        """
        Tests results from `LevMarLSQFitter` with weights.
        """
        # part 1: weights are equal to 1
        fitter = LevMarLSQFitter()
        model = fitter(self.gauss, self.xdata, self.ydata,
                       estimate_jacobian=True)
        withw = fitter(self.gauss, self.xdata, self.ydata,
                       estimate_jacobian=True, weights=np.ones_like(self.xdata))

        assert_allclose(model.parameters, withw.parameters, rtol=10 ** (-4))

        # part 2: weights are 0 or 1 (effectively, they are a mask)
        weights = np.zeros_like(self.xdata)
        weights[::2] = 1.
        mask = weights >= 1.

        model = fitter(self.gauss, self.xdata[mask], self.ydata[mask],
                       estimate_jacobian=True)
        withw = fitter(self.gauss, self.xdata, self.ydata,
                       estimate_jacobian=True, weights=weights)

        assert_allclose(model.parameters, withw.parameters, rtol=10 ** (-4))

    @pytest.mark.parametrize('fitter_class', fitters)
    def test_fitter_against_LevMar(self, fitter_class):
        """Tests results from non-linear fitters against `LevMarLSQFitter`."""

        levmar = LevMarLSQFitter()
        fitter = fitter_class()
        with ignore_non_integer_warning():
            new_model = fitter(self.gauss, self.xdata, self.ydata)
        model = levmar(self.gauss, self.xdata, self.ydata)
        assert_allclose(model.parameters, new_model.parameters,
                        rtol=10 ** (-4))

    def test_LSQ_SLSQP_with_constraints(self):
        """
        Runs `LevMarLSQFitter` and `SLSQPLSQFitter` on a model with
        constraints.
        """

        g1 = models.Gaussian1D(100, 5, stddev=1)
        g1.mean.fixed = True
        fitter = LevMarLSQFitter()
        fslsqp = SLSQPLSQFitter()
        with ignore_non_integer_warning():
            slsqp_model = fslsqp(g1, self.xdata, self.ydata)
        model = fitter(g1, self.xdata, self.ydata)
        assert_allclose(model.parameters, slsqp_model.parameters,
                        rtol=10 ** (-4))

    def test_simplex_lsq_fitter(self):
        """A basic test for the `SimplexLSQ` fitter."""

        class Rosenbrock(Fittable2DModel):
            a = Parameter()
            b = Parameter()

            @staticmethod
            def evaluate(x, y, a, b):
                return (a - x) ** 2 + b * (y - x ** 2) ** 2

        x = y = np.linspace(-3.0, 3.0, 100)
        with NumpyRNGContext(_RANDOM_SEED):
            z = Rosenbrock.evaluate(x, y, 1.0, 100.0)
            z += np.random.normal(0., 0.1, size=z.shape)

        fitter = SimplexLSQFitter()
        r_i = Rosenbrock(1, 100)
        r_f = fitter(r_i, x, y, z)

        assert_allclose(r_f.parameters, [1.0, 100.0], rtol=1e-2)

    def test_param_cov(self):
        """
        Tests that the 'param_cov' fit_info entry gets the right answer for
        *linear* least squares, where the answer is exact
        """

        a = 2
        b = 100

        with NumpyRNGContext(_RANDOM_SEED):
            x = np.linspace(0, 1, 100)
            # y scatter is amplitude ~1 to make sure covarience is
            # non-negligible
            y = x*a + b + np.random.randn(len(x))

        # first compute the ordinary least squares covariance matrix
        X = np.matrix(np.vstack([x, np.ones(len(x))]).T)
        beta = np.linalg.inv(X.T * X) * X.T * np.matrix(y).T
        s2 = np.sum((y - (X * beta).A.ravel())**2) / (len(y) - len(beta))
        olscov = np.linalg.inv(X.T * X) * s2

        # now do the non-linear least squares fit
        mod = models.Linear1D(a, b)
        fitter = LevMarLSQFitter()

        fmod = fitter(mod, x, y)

        assert_allclose(fmod.parameters, beta.A.ravel())
        assert_allclose(olscov, fitter.fit_info['param_cov'])


@pytest.mark.skipif('not HAS_PKG')
class TestEntryPoint:
    """Tests population of fitting with entry point fitters"""

    def setup_class(self):
        self.exception_not_thrown = Exception("The test should not have gotten here. There was no exception thrown")

    def successfulimport(self):
            # This should work
            class goodclass(Fitter):
                __name__ = "GoodClass"
            return goodclass

    def raiseimporterror(self):
        #  This should fail as it raises an Import Error
        raise ImportError

    def returnbadfunc(self):
        def badfunc():
            # This should import but it should fail type check
            pass
        return badfunc

    def returnbadclass(self):
        # This should import But it should fail subclass type check
        class badclass:
            pass
        return badclass

    def test_working(self):
        """This should work fine"""
        mock_entry_working = mock.create_autospec(EntryPoint)
        mock_entry_working.name = "Working"
        mock_entry_working.load = self.successfulimport
        populate_entry_points([mock_entry_working])

    def test_import_error(self):
        """This raises an import error on load to test that it is handled correctly"""
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                mock_entry_importerror = mock.create_autospec(EntryPoint)
                mock_entry_importerror.name = "IErr"
                mock_entry_importerror.load = self.raiseimporterror
                populate_entry_points([mock_entry_importerror])
            except AstropyUserWarning as w:
                if "ImportError" in w.args[0]:  # any error for this case should have this in it.
                    pass
                else:
                    raise w
            else:
                raise self.exception_not_thrown

    def test_bad_func(self):
        """This returns a function which fails the type check"""
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                mock_entry_badfunc = mock.create_autospec(EntryPoint)
                mock_entry_badfunc.name = "BadFunc"
                mock_entry_badfunc.load = self.returnbadfunc
                populate_entry_points([mock_entry_badfunc])
            except AstropyUserWarning as w:
                if "Class" in w.args[0]:  # any error for this case should have this in it.
                    pass
                else:
                    raise w
            else:
                raise self.exception_not_thrown

    def test_bad_class(self):
        """This returns a class which doesn't inherient from fitter """
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                mock_entry_badclass = mock.create_autospec(EntryPoint)
                mock_entry_badclass.name = "BadClass"
                mock_entry_badclass.load = self.returnbadclass
                populate_entry_points([mock_entry_badclass])
            except AstropyUserWarning as w:
                if 'modeling.Fitter' in w.args[0]:  # any error for this case should have this in it.
                    pass
                else:
                    raise w
            else:
                raise self.exception_not_thrown


@pytest.mark.skipif('not HAS_SCIPY')
class Test1DFittingWithOutlierRemoval:
    def setup_class(self):
        self.x = np.linspace(-5., 5., 200)
        self.model_params = (3.0, 1.3, 0.8)

        def func(p, x):
            return p[0]*np.exp(-0.5*(x - p[1])**2/p[2]**2)

        self.y = func(self.model_params, self.x)

    def test_with_fitters_and_sigma_clip(self):
        import scipy.stats as stats

        np.random.seed(0)
        c = stats.bernoulli.rvs(0.25, size=self.x.shape)
        self.y += (np.random.normal(0., 0.2, self.x.shape) +
                   c*np.random.normal(3.0, 5.0, self.x.shape))

        g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
        # test with Levenberg-Marquardt Least Squares fitter
        fit = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                        niter=3, sigma=3.0)
        _, fitted_model = fit(g_init, self.x, self.y)
        assert_allclose(fitted_model.parameters, self.model_params, rtol=1e-1)
        # test with Sequential Least Squares Programming fitter
        fit = FittingWithOutlierRemoval(SLSQPLSQFitter(), sigma_clip,
                                        niter=3, sigma=3.0)
        _, fitted_model = fit(g_init, self.x, self.y)
        assert_allclose(fitted_model.parameters, self.model_params, rtol=1e-1)
        # test with Simplex LSQ fitter
        fit = FittingWithOutlierRemoval(SimplexLSQFitter(), sigma_clip,
                                        niter=3, sigma=3.0)
        _, fitted_model = fit(g_init, self.x, self.y)
        assert_allclose(fitted_model.parameters, self.model_params, atol=1e-1)


@pytest.mark.skipif('not HAS_SCIPY')
class Test2DFittingWithOutlierRemoval:
    def setup_class(self):
        self.y, self.x = np.mgrid[-3:3:128j, -3:3:128j]
        self.model_params = (3.0, 1.0, 0.0, 0.8, 0.8)

        def Gaussian_2D(p, pos):
            return p[0]*np.exp(-0.5*(pos[0] - p[2])**2 / p[4]**2 -
                               0.5*(pos[1] - p[1])**2 / p[3]**2)

        self.z = Gaussian_2D(self.model_params, np.array([self.y, self.x]))

    def initial_guess(self, data, pos):
        y = pos[0]
        x = pos[1]

        """computes the centroid of the data as the initial guess for the
        center position"""

        wx = x * data
        wy = y * data
        total_intensity = np.sum(data)
        x_mean = np.sum(wx) / total_intensity
        y_mean = np.sum(wy) / total_intensity

        x_to_pixel = x[0].size / (x[x[0].size - 1][x[0].size - 1] - x[0][0])
        y_to_pixel = y[0].size / (y[y[0].size - 1][y[0].size - 1] - y[0][0])
        x_pos = np.around(x_mean * x_to_pixel + x[0].size / 2.).astype(int)
        y_pos = np.around(y_mean * y_to_pixel + y[0].size / 2.).astype(int)

        amplitude = data[y_pos][x_pos]

        return amplitude, x_mean, y_mean

    def test_with_fitters_and_sigma_clip(self):
        import scipy.stats as stats

        np.random.seed(0)
        c = stats.bernoulli.rvs(0.25, size=self.z.shape)
        self.z += (np.random.normal(0., 0.2, self.z.shape) +
                   c*np.random.normal(self.z, 2.0, self.z.shape))

        guess = self.initial_guess(self.z, np.array([self.y, self.x]))
        g2_init = models.Gaussian2D(amplitude=guess[0], x_mean=guess[1],
                                    y_mean=guess[2], x_stddev=0.75,
                                    y_stddev=1.25)

        # test with Levenberg-Marquardt Least Squares fitter
        fit = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                        niter=3, sigma=3.)
        _, fitted_model = fit(g2_init, self.x, self.y, self.z)
        assert_allclose(fitted_model.parameters[0:5], self.model_params,
                        atol=1e-1)
        # test with Sequential Least Squares Programming fitter
        fit = FittingWithOutlierRemoval(SLSQPLSQFitter(), sigma_clip, niter=3,
                                        sigma=3.)
        _, fitted_model = fit(g2_init, self.x, self.y, self.z)
        assert_allclose(fitted_model.parameters[0:5], self.model_params,
                        atol=1e-1)
        # test with Simplex LSQ fitter
        fit = FittingWithOutlierRemoval(SimplexLSQFitter(), sigma_clip,
                                        niter=3, sigma=3.)
        _, fitted_model = fit(g2_init, self.x, self.y, self.z)
        assert_allclose(fitted_model.parameters[0:5], self.model_params,
                        atol=1e-1)


def test_1d_set_fitting_with_outlier_removal():
    """Test model set fitting with outlier removal (issue #6819)"""

    poly_set = models.Polynomial1D(2, n_models=2)

    fitter = FittingWithOutlierRemoval(LinearLSQFitter(),
                                       sigma_clip, sigma=2.5, niter=3,
                                       cenfunc=np.ma.mean, stdfunc=np.ma.std)

    x = np.arange(10)
    y = np.array([2.5*x - 4, 2*x*x + x + 10])
    y[1,5] = -1000  # outlier

    filt_y, poly_set = fitter(poly_set, x, y)

    assert_allclose(poly_set.c0, [-4., 10.], atol=1e-14)
    assert_allclose(poly_set.c1, [2.5, 1.], atol=1e-14)
    assert_allclose(poly_set.c2, [0., 2.], atol=1e-14)


def test_2d_set_axis_2_fitting_with_outlier_removal():
    """Test fitting 2D model set (axis 2) with outlier removal (issue #6819)"""

    poly_set = models.Polynomial2D(1, n_models=2, model_set_axis=2)

    fitter = FittingWithOutlierRemoval(LinearLSQFitter(),
                                       sigma_clip, sigma=2.5, niter=3,
                                       cenfunc=np.ma.mean, stdfunc=np.ma.std)

    y, x = np.mgrid[0:5, 0:5]
    z = np.rollaxis(np.array([x+y, 1-0.1*x+0.2*y]), 0, 3)
    z[3,3:5,0] = 100.   # outliers

    filt_z, poly_set = fitter(poly_set, x, y, z)

    assert_allclose(poly_set.c0_0, [[[0., 1.]]], atol=1e-14)
    assert_allclose(poly_set.c1_0, [[[1., -0.1]]], atol=1e-14)
    assert_allclose(poly_set.c0_1, [[[1., 0.2]]], atol=1e-14)


@pytest.mark.skipif('not HAS_SCIPY')
class TestWeightedFittingWithOutlierRemoval:
    """Issue #7020 """

    def setup_class(self):
        # values of x,y not important as we fit y(x,y) = p0 model here
        self.y, self.x = np.mgrid[0:20, 0:20]
        self.z = np.mod(self.x + self.y, 2) * 2 - 1 # -1,1 chessboard
        self.weights = np.mod(self.x + self.y, 2) * 2 + 1 # 1,3 chessboard
        self.z[0,0] = 1000.0 # outlier
        self.z[0,1] = 1000.0 # outlier
        self.x1d = self.x.flatten()
        self.z1d = self.z.flatten()
        self.weights1d = self.weights.flatten()

    def test_1d_without_weights_without_sigma_clip(self):
        model = models.Polynomial1D(0)
        fitter = LinearLSQFitter()
        fit = fitter(model, self.x1d, self.z1d)
        assert_allclose(fit.parameters[0], self.z1d.mean(), atol=10**(-2))

    def test_1d_without_weights_with_sigma_clip(self):
        model = models.Polynomial1D(0)
        fitter = FittingWithOutlierRemoval(LinearLSQFitter(), sigma_clip,
                                           niter=3, sigma=3.)
        filtered, fit = fitter(model, self.x1d, self.z1d)
        assert(filtered.count() == self.z1d.size - 2)
        assert(filtered.mask[0] and filtered.mask[1])
        assert_allclose(fit.parameters[0], 0.0, atol=10**(-2)) # with removed outliers mean is 0.0

    def test_1d_with_weights_without_sigma_clip(self):
        model = models.Polynomial1D(0)
        fitter = LinearLSQFitter()
        fit = fitter(model, self.x1d, self.z1d, weights=self.weights1d)
        assert(fit.parameters[0] > 1.0)     # outliers pulled it high

    def test_1d_with_weights_with_sigma_clip(self):
        """smoke test for #7020 - fails without fitting.py patch because weights does not propagate"""
        model = models.Polynomial1D(0)
        fitter = FittingWithOutlierRemoval(LinearLSQFitter(), sigma_clip,
                                           niter=3, sigma=3.)
        filtered, fit = fitter(model, self.x1d, self.z1d, weights=self.weights1d)
        assert(fit.parameters[0] > 10**(-2))  # weights pulled it > 0
        assert(fit.parameters[0] < 1.0)       # outliers didn't pull it out of [-1:1] because they had been removed

    def test_1d_set_with_common_weights_with_sigma_clip(self):
        """added for #6819 (1D model set with weights in common)"""
        model = models.Polynomial1D(0, n_models=2)
        fitter = FittingWithOutlierRemoval(LinearLSQFitter(), sigma_clip,
                                           niter=3, sigma=3.)
        z1d = np.array([self.z1d, self.z1d])

        filtered, fit = fitter(model, self.x1d, z1d, weights=self.weights1d)
        assert_allclose(fit.parameters, [0.8, 0.8], atol=1e-14)

    def test_2d_without_weights_without_sigma_clip(self):
        model = models.Polynomial2D(0)
        fitter = LinearLSQFitter()
        fit = fitter(model, self.x, self.y, self.z)
        assert_allclose(fit.parameters[0], self.z.mean(), atol=10**(-2))

    def test_2d_without_weights_with_sigma_clip(self):
        model = models.Polynomial2D(0)
        fitter = FittingWithOutlierRemoval(LinearLSQFitter(), sigma_clip,
                                           niter=3, sigma=3.)
        filtered, fit = fitter(model, self.x, self.y, self.z)
        assert(filtered.count() == self.z.size - 2)
        assert(filtered.mask[0,0] and filtered.mask[0,1])
        assert_allclose(fit.parameters[0], 0.0, atol=10**(-2))

    def test_2d_with_weights_without_sigma_clip(self):
        model = models.Polynomial2D(0)
        fitter = LevMarLSQFitter() # LinearLSQFitter doesn't handle weights properly in 2D
        fit = fitter(model, self.x, self.y, self.z, weights=self.weights)
        assert(fit.parameters[0] > 1.0)     # outliers pulled it high

    def test_2d_with_weights_with_sigma_clip(self):
        """smoke test for #7020 - fails without fitting.py patch because weights does not propagate"""
        model = models.Polynomial2D(0)
        fitter = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                           niter=3, sigma=3.)
        filtered, fit = fitter(model, self.x, self.y, self.z, weights=self.weights)
        assert(fit.parameters[0] > 10**(-2))  # weights pulled it > 0
        assert(fit.parameters[0] < 1.0)       # outliers didn't pull it out of [-1:1] because they had been removed


@pytest.mark.skipif('not HAS_SCIPY')
def test_fitters_with_weights():
    """Issue #5737 """
    Xin, Yin = np.mgrid[0:21, 0:21]
    fitter = LevMarLSQFitter()

    with NumpyRNGContext(_RANDOM_SEED):
        zsig = np.random.normal(0, 0.01, size=Xin.shape)

    # Non-linear model
    g2 = models.Gaussian2D(10, 10, 9, 2, 3)
    z = g2(Xin, Yin)
    gmod = fitter(models.Gaussian2D(15, 7, 8, 1.3, 1.2), Xin, Yin, z + zsig)
    assert_allclose(gmod.parameters, g2.parameters, atol=10 ** (-2))

    # Linear model
    p2 = models.Polynomial2D(3)
    p2.parameters = np.arange(10)/1.2
    z = p2(Xin, Yin)
    pmod = fitter(models.Polynomial2D(3), Xin, Yin, z + zsig)
    assert_allclose(pmod.parameters, p2.parameters, atol=10 ** (-2))


@pytest.mark.skipif('not HAS_SCIPY')
def test_fitters_interface():
    """
    Test that **kwargs work with all optimizers.
    This is a basic smoke test.
    """
    levmar = LevMarLSQFitter()
    slsqp = SLSQPLSQFitter()
    simplex = SimplexLSQFitter()

    kwargs = {'maxiter': 77, 'verblevel': 1, 'epsilon': 1e-2, 'acc': 1e-6}
    simplex_kwargs = {'maxiter': 77, 'verblevel': 1, 'acc': 1e-6}
    model = models.Gaussian1D(10, 4, .3)
    x = np.arange(21)
    y = model(x)

    slsqp_model = slsqp(model, x, y, **kwargs)
    simplex_model = simplex(model, x, y, **simplex_kwargs)
    kwargs.pop('verblevel')
    lm_model = levmar(model, x, y, **kwargs)

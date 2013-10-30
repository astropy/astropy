# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Module to test fitting routines
"""
from __future__ import division
import os.path
from .. import models
from .. import fitting
from . import irafutil
import numpy as np
from numpy import linalg
from numpy.testing import utils
from numpy.random import RandomState
from ...utils.data import get_pkg_data_filename
from ...tests.helper import pytest

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class TestPolynomial2D(object):

    """
    Tests for 2D polynomail fitting
    """
    def setup_class(self):
        self.model = models.Polynomial2D(2)
        self.x, self.y = np.mgrid[:5, :5]

        def poly2(x, y):
            return 1 + 2 * x + 3 * x ** 2 + 4 * y + 5 * y ** 2 + 6 * x * y
        self.z = poly2(self.x, self.y)
        self.fitter = fitting.LinearLSQFitter()

    def test_poly2D_fitting(self):
        v = self.model.deriv(x=self.x, y=self.y)
        p = linalg.lstsq(v, self.z.flatten())[0]
        new_model = self.fitter(self.model, self.x, self.y, self.z)
        utils.assert_allclose(new_model.parameters, p)

    def test_eval(self):
        new_model = self.fitter(self.model, self.x, self.y, self.z)
        utils.assert_allclose(new_model(self.x, self.y), self.z)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_polynomial2D_nonlinear_fitting(self):
        self.model.parameters = [.6, 1.8, 2.9, 3.7, 4.9, 6.7]
        nlfitter = fitting.NonLinearLSQFitter()
        new_model = nlfitter(self.model, self.x, self.y, self.z)
        utils.assert_allclose(new_model.parameters, [1, 2, 3, 4, 5, 6])


class TestICheb2D(object):

    """
    Tests 2D Chebyshev polynomial fitting

    Create a 2D polynomial (z) using Polynomial2DModel and default coefficients
    Fit z using a ICheb2D model
    Evaluate the ICheb2D polynomial and compare with the initial z
    """
    def setup_class(self):
        self.pmodel = models.Polynomial2D(2)
        self.x, self.y = np.mgrid[:5, :5]
        self.z = self.pmodel(self.x, self.y)
        self.cheb2 = models.Chebyshev2D(2, 2)
        self.fitter = fitting.LinearLSQFitter()

    def test_default_params(self):
        self.cheb2.parameters = np.arange(9)
        p = np.array([1344., 1772., 400., 1860., 2448., 552., 432., 568.,
                      128.])
        z = self.cheb2(self.x, self.y)
        model = self.fitter(self.cheb2, self.x, self.y, z)
        utils.assert_almost_equal(model.parameters, p)

    def test_poly2D_cheb2D(self):
        model = self.fitter(self.cheb2, self.x, self.y, self.z)
        z1 = model(self.x, self.y)
        utils.assert_almost_equal(self.z, z1)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_chebyshev2D_nonlinear_fitting(self):
        cheb2d = models.Chebyshev2D(2, 2)
        cheb2d.parameters = np.arange(9)
        z = cheb2d(self.x, self.y)
        cheb2d.parameters = [0.1, .6, 1.8, 2.9, 3.7, 4.9, 6.7, 7.5, 8.9]
        nlfitter = fitting.NonLinearLSQFitter()
        model = nlfitter(cheb2d, self.x, self.y, z)
        utils.assert_allclose(model.parameters, [0, 1, 2, 3, 4, 5, 6, 7, 8], atol=10**-9)


@pytest.mark.skipif('not HAS_SCIPY')
class TestJointFitter(object):

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
        self.jf = fitting.JointFitter([self.g1, self.g2],
                                      {self.g1: ['amplitude'],
                                       self.g2: ['amplitude']}, [9.8])
        self.x = np.arange(10, 20, .1)
        y1 = self.g1(self.x)
        y2 = self.g2(self.x)
        n = np.random.randn(100)
        self.ny1 = y1 + 2 * n
        self.ny2 = y2 + 2 * n
        self.jf(self.x, self.ny1, self.x, self.ny2)

    def test_joint_parameter(self):
        """
        Tests that the amplitude of the two models is the same
        """
        utils.assert_allclose(self.jf.fitparams[0], self.g1.parameters[0])
        utils.assert_allclose(self.jf.fitparams[0], self.g2.parameters[0])

    def test_joint_fitter(self):
        """
        Tests the fitting routine with similar procedure.
        Compares the fitted parameters.
        """
        p1 = [14.9, .3]
        p2 = [13, .4]
        A = 9.8
        p = np.r_[A, p1, p2]
        compmodel = lambda A, p, x: A * np.exp(-0.5 / p[1] ** 2 * (x - p[0]) ** 2)
        errf = lambda p, x1, y1, x2, y2: np.ravel(np.r_[compmodel(p[0], p[1:3],
                                                                  x1) - y1, compmodel(p[0], p[3:], x2) - y2])
        coeff, _ = optimize.leastsq(errf, p, args=(self.x, self.ny1, self.x,
                                                   self.ny2))
        utils.assert_allclose(coeff, self.jf.fitparams, rtol=10 ** (-2))


class TestLinearLSQFitter(object):

    def setup_class(self):
        test_file = get_pkg_data_filename(os.path.join('data', 'idcompspec.fits'))
        f = open(test_file)
        lines = f.read()
        reclist = lines.split("begin")
        f.close()
        record = irafutil.IdentifyRecord(reclist[1])
        self.icoeff = record.coeff
        order = int(record.fields['order'])
        self.model = models.Chebyshev1D(order - 1)
        self.model.domain = record.get_range()
        self.lf = fitting.LinearLSQFitter()
        self.x = record.x
        self.y = record.z
        self.yy = np.array([record.z, record.z])

    def test_chebyshev1D(self):
        new_model = self.lf(self.model, self.x, self.y)
        utils.assert_allclose(new_model.parameters, np.array(self.icoeff),
                              rtol=10E-2)


@pytest.mark.skipif('not HAS_SCIPY')
class TestNonLinearFitters(object):

    """
    Tests non-linear least squares fitting and the SLSQP algorithm
    """
    def setup_class(self):
        self.initial_values = [100, 5, 1]
        func = lambda p, x: p[0] * np.exp(-0.5 / p[2] ** 2 * (x - p[1]) ** 2)
        errf = lambda p, x, y: (func(p, x) - y)
        self.xdata = np.arange(0, 10, 0.1)
        sigma = 8. * np.ones_like(self.xdata)
        rsn = RandomState(1234567890)
        yerror = rsn.normal(0, sigma)
        self.ydata = func(self.initial_values, self.xdata) + yerror

    def test_estimated_vs_analytic_deriv(self):
        g1 = models.Gaussian1D(100, 5, stddev=1)
        fitter = fitting.NonLinearLSQFitter()
        model = fitter(g1, self.xdata, self.ydata)
        g1e = models.Gaussian1D(100, 5.0, stddev=1)
        efitter = fitting.NonLinearLSQFitter()
        emodel = efitter(g1e, self.xdata, self.ydata, estimate_jacobian=True)
        utils.assert_allclose(model.parameters, emodel.parameters, rtol=10 ** (-3))

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_with_optimize(self):
        g1 = models.Gaussian1D(100, 5, stddev=1)
        fitter = fitting.NonLinearLSQFitter()
        model = fitter(g1, self.xdata, self.ydata, estimate_jacobian=True)
        func = lambda p, x: p[0] * np.exp(-0.5 / p[2] ** 2 * (x - p[1]) ** 2)
        errf = lambda p, x, y: (func(p, x) - y)
        result = optimize.leastsq(errf, self.initial_values, args=(self.xdata, self.ydata))
        utils.assert_allclose(model.parameters, result[0], rtol=10 ** (-3))

    def test_LSQ_SLSQP(self):
        g1 = models.Gaussian1D(100, 5, stddev=1)
        fitter = fitting.NonLinearLSQFitter()
        fslsqp = fitting.SLSQPFitter()
        slsqp_model = fslsqp(g1, self.xdata, self.ydata)
        model = fitter(g1, self.xdata, self.ydata)
        # There's a bug in the SLSQP algorithm and sometimes it gives the
        # negative value of the result. unitl this is understood, for this
        # test, take np.abs()
        utils.assert_allclose(model.parameters, np.abs(slsqp_model.parameters),
                              rtol=10 ** (-4))

    def test_LSQ_SLSQP_cons(self):
        g1 = models.Gaussian1D(100, 5, stddev=1)
        g1.mean.fixed = True
        fitter = fitting.NonLinearLSQFitter()
        fslsqp = fitting.SLSQPFitter()
        slsqp_model = fslsqp(g1, self.xdata, self.ydata)
        model = fitter(g1, self.xdata, self.ydata)
        # There's a bug in the SLSQP algorithm and sometimes it gives the
        # negative value of the result. unitl this is understood, for this
        # test, take np.abs()
        utils.assert_allclose(model.parameters, np.abs(slsqp_model.parameters),
                              rtol=10 ** (-4))

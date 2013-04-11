# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Module to test fitting routines
"""
from __future__ import division
import os.path
from .. import models, fitting
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

class TestPoly2D(object):
    """
    Tests for 2D polynomail fitting
    """
    def setup_class(self):
        self.model = models.Poly2DModel(2)
        self.x, self.y = np.mgrid[:5, :5]
        def poly2(x, y):
            return 1+2*x+3*x**2+4*y+5*y**2+6*x*y
        self.z = poly2(self.x, self.y)
        self.fitter = fitting.LinearLSQFitter(self.model)

    def test_poly2D_fitting(self):
        v = self.model.deriv(self.x, self.y)
        p = linalg.lstsq(v, self.z.flatten())[0]
        self.fitter(self.x, self.y, self.z)
        utils.assert_allclose(self.model.parameters, p)
                         
    def test_eval(self):
        self.fitter(self.x, self.y, self.z)
        utils.assert_allclose(self.model(self.x, self.y), self.z)
        
class TestICheb2D(object):
    """
    Tests 2D Chebyshev polynomial fitting
    
    Create a 2D polynomial (z) using Poly2DModel and default coefficients
    Fit z using a ICheb2D model
    Evaluate the ICheb2D polynomial and compare with the initial z
    """
    def setup_class(self):
        self.pmodel = models.Poly2DModel(2)
        self.x, self.y = np.mgrid[:5, :5]
        self.z = self.pmodel(self.x, self.y)
        self.ichb = models.ICheb2DModel(2, 2)
        self.fitter = fitting.LinearLSQFitter(self.ichb)
        
    def test_default_pars(self):
        self.ichb.parameters = np.arange(9)
        p = np.array([ 1344.,  1772.,   400.,  1860.,  2448.,   552.,   432.,   568.,
         128.])
        z = self.ichb(self.x, self.y)
        self.fitter(self.x, self.y, z)
        utils.assert_almost_equal(self.ichb.parameters, p)
        
    def test_poly2D_icheb2D(self):
        self.fitter(self.x, self.y, self.z)
        z1 = self.ichb(self.x, self.y)
        utils.assert_almost_equal(self.z, z1)

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
        self.g1 = models.Gauss1DModel(10, 14.9, .3)
        self.g2 = models.Gauss1DModel(10, 13, .4)
        self.jf = fitting.JointFitter([self.g1, self.g2], 
                {self.g1:['amplitude'], self.g2:['amplitude']}, [9.8])
        self.x = np.arange(10, 20, .1)
        y1 = self.g1(self.x)
        y2 = self.g2(self.x)
        n = np.random.randn(100)
        self.ny1 = y1 + 2*n
        self.ny2 = y2 + 2*n
        self.jf(self.x, self.ny1, self.x, self.ny2)
        
    def test_joint_parameter(self):
        """
        Tests that the amplitude of the two models is the same
        """
        utils.assert_allclose(self.jf.fitpars[0], self.g1.parameters[0])
        utils.assert_allclose(self.jf.fitpars[0], self.g2.parameters[0])
        
    def test_joint_fitter(self):
        """
        Tests the fitting routine with similar procedure.
        Compares the fitted parameters.
        """
        p1 = [14.9, .3]
        p2 = [13, .4]
        A = 9.8
        p = np.r_[A, p1, p2]
        compmodel = lambda A, p, x: A* np.exp((-(1/(p[1]**2)) * (x-p[0])**2))
        errf = lambda p, x1, y1, x2, y2: np.ravel(np.r_[compmodel(p[0], p[1:3], 
                                x1) -y1, compmodel(p[0], p[3:], x2) - y2])
        coeff, finfo = optimize.leastsq(errf, p, args=(self.x, self.ny1, self.x, 
                                self.ny2))
        utils.assert_allclose(coeff, self.jf.fitpars, rtol=10**(-2))
        
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
        self.model = models.ChebyshevModel(order-1)
        self.model.domain = record.get_range()
        self.lf = fitting.LinearLSQFitter(self.model)
        self.x = record.x
        self.y = record.z
        self.yy = np.array([record.z, record.z])
        
    def test_chebyshev1D(self):
        self.lf(self.x, self.y)
        utils.assert_allclose(self.model.parameters, np.array(self.icoeff),
                              rtol=10E-2)
        
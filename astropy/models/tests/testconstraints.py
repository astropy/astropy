from __future__ import division

from .. import models, fitting
import numpy as np
from scipy import optimize
from numpy.testing import utils

class TestNonLinearConstraints:
    def setUp(self):
        self.g1 = models.Gauss1DModel(10, 14.9, xsigma=.3)
        self.g2 = models.Gauss1DModel(10, 13, xsigma=.4)
        self.x = np.arange(10, 20, .1)
        self.y1 = self.g1(self.x)
        self.y2 = self.g2(self.x)
        n = np.random.randn(100)
        self.ny1 = self.y1 + 2*n
        self.ny2 = self.y2 + 2*n
        
    def testFixedPar(self):
        g1 = models.Gauss1DModel(10, 14.9, xsigma=.3, fixed={'amplitude':True})
        func = lambda p, x: 10* np.exp((-(1/(p[2]**2)) * (x-p[1])**2))
        errf = lambda p, x, y: func(p, x)-y
        p0 = [10, 14.5, 0.3]
        fitpar, s = optimize.leastsq(errf, p0, args=(self.x, self.ny1))
        fitter = fitting.NonLinearLSQFitter(g1)
        fitter(self.x, self.ny1)
        utils.assert_allclose(g1.parameters, fitpar, rtol=5*10**(-3))
        
    def testTiedPar(self):
        
        def tied(model):
            xcen = 50*model.xsigma[0]
            return xcen
        g1 = models.Gauss1DModel(10, 14.9, xsigma=.3, tied={'xcen':tied})
        fitter = fitting.NonLinearLSQFitter(g1)
        fitter(self.x, self.ny1)
        utils.assert_allclose(g1.xcen, 50*g1.xsigma[0], rtol=10**(-5))
        
    def testJointFitter(self):
        g1 = models.Gauss1DModel(10, 14.9, xsigma=.3)
        g2 = models.Gauss1DModel(10, 13, xsigma=.4)
        jf = fitting.JointFitter([g1, g2], {g1:['amplitude'], 
                                         g2:['amplitude']}, [9.8])
        x = np.arange(10, 20, .1)
        y1 = g1(x)
        y2 = g2(x)
        n = np.random.randn(100)
        ny1 = y1 + 2*n
        ny2 = y2 + 2*n
        jf(x, ny1, x, ny2)
        p1 = [14.9, .3]
        p2 = [13, .4]
        A = 9.8
        p = np.r_[A, p1, p2]
        compmodel = lambda A, p, x: A* np.exp((-(1/(p[1]**2)) * 
                                             (x-p[0])**2))
        errf = lambda p, x1, y1, x2, y2: np.ravel(
            np.r_[compmodel(p[0], p[1:3], x1) - y1, 
                  compmodel(p[0], p[3:], x2) - y2])
        fitpars, s = optimize.leastsq(errf, p, args=(x, ny1, x, ny2))
        utils.assert_allclose(jf.fitpars, fitpars, rtol=10**(-5))
        utils.assert_allclose(g1.amplitude, g2.amplitude)
        
    def testNoConstraints(self):
        g1 = models.Gauss1DModel(9.9, 14.5, xsigma=.3)
        func = lambda p, x: p[0]* np.exp((-(1/(p[2]**2)) * (x-p[1])**2))
        errf = lambda p, x, y: func(p, x) - y
        p0 = [9.9, 14.5, 0.3]
        y = g1(self.x)
        n = np.random.randn(100)
        ny = y + n
        fitpar, s = optimize.leastsq(errf, p0, args=(self.x, n+y))
        fitter = fitting.NonLinearLSQFitter(g1)
        fitter(self.x, n+y)
        utils.assert_allclose(g1.parameters, fitpar, rtol=5*10**(-3))
        
class TestLinearConstraints:
    def setUp(self):
        self.p1 = models.Poly1DModel(4)
        self.p1.c0 = 0
        self.p1.c1 = 0
        self.p1.window = [0., 9.]
        self.x = np.arange(10)
        self.y = self.p1(self.x)
        n = np.random.randn(10)
        self.ny = self.y + n
        
    def test(self):
        self.p1.c0.fixed = True
        self.p1.c1.fixed = True
        pfit=fitting.LinearLSQFitter(self.p1)
        pfit(self.x, self.y)
        utils.assert_allclose(self.y, self.p1(self.x))
        
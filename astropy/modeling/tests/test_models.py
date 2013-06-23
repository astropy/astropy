# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for model evaluation.
Compare the results of some models with other programs.
"""
from __future__ import division
from .. import models
from ..core import *
import numpy as np
from numpy.testing import utils
from ...tests.helper import pytest
from .. import fitting

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class TestSComposite(object):

    """
    Test composite models evaluation in series
    """
    def setup_class(self):
        self.x, self.y = np.mgrid[:5, :5]
        self.p1 = models.Poly1DModel(3)
        self.p11 = models.Poly1DModel(3)
        self.p2 = models.Poly2DModel(3)

    def test_single_array_input(self):
        scomptr = SCompositeModel([self.p1, self.p11])
        sresult = scomptr(self.x)
        xx = self.p11(self.p1(self.x))
        utils.assert_almost_equal(xx, sresult)

    def test_labeledinput_1(self):
        ado = LabeledInput([self.x, self.y], ['x', 'y'])
        scomptr = SCompositeModel([self.p2, self.p1],
                                  [['x', 'y'], ['z']],
                                  [['z'], ['z']])
        sresult = scomptr(ado)
        z = self.p2(self.x, self.y)
        z1 = self.p1(z)
        utils.assert_almost_equal(z1, sresult.z)

    def test_labeledinput_2(self):
        labeled_input = LabeledInput([self.x, self.y], ['x', 'y'])
        rot = models.MatrixRotation2D(angle=23.4)
        offx = models.ShiftModel(-2)
        offy = models.ShiftModel(1.2)
        scomptr = SCompositeModel([rot, offx, offy],
                                  [['x', 'y'], ['x'], ['y']],
                                  [['x', 'y'], ['x'], ['y']])
        sresult = scomptr(labeled_input)
        x, y = rot(self.x, self.y)
        x = offx(x)
        y = offy(y)
        utils.assert_almost_equal(x, sresult.x)
        utils.assert_almost_equal(y, sresult.y)

    def test_labeledinput_3(self):
        labeled_input = LabeledInput([2, 4.5], ['x', 'y'])
        rot = models.MatrixRotation2D(angle=23.4)
        offx = models.ShiftModel(-2)
        offy = models.ShiftModel(1.2)
        scomptr = SCompositeModel([rot, offx, offy],
                                  [['x', 'y'], ['x'], ['y']],
                                  [['x', 'y'], ['x'], ['y']])
        sresult = scomptr(labeled_input)
        x, y = rot(2, 4.5)
        x = offx(x)
        y = offy(y)
        utils.assert_almost_equal(x, sresult.x)
        utils.assert_almost_equal(y, sresult.y)

    def test_multiple_input(self):
        rot = models.MatrixRotation2D(angle=-60)
        scomp = SCompositeModel([rot, rot])
        xx, yy = scomp(self.x, self.y)
        iscomp = scomp.inverse()
        x1, y1 = iscomp(xx, yy)
        utils.assert_almost_equal(x1, self.x)
        utils.assert_almost_equal(y1, self.y)


class TestPComposite(object):

    """
    Test composite models evaluation in parallel
    """
    def setup_class(self):
        self.x, self.y = np.mgrid[:5, :5]
        self.p1 = models.Poly1DModel(3)
        self.p11 = models.Poly1DModel(3)
        self.p2 = models.Poly2DModel(3)

    def test_single_array_input(self):
        pcomptr = PCompositeModel([self.p1, self.p11])
        presult = pcomptr(self.x)
        delta11 = self.p11(self.x)
        delta1 = self.p1(self.x)
        xx = self.x + delta1 + delta11
        utils.assert_almost_equal(xx, presult)

    def test_labeledinput(self):
        ado = LabeledInput([self.x, self.y], ['x', 'y'])
        pcomptr = PCompositeModel([self.p1, self.p11], inmap=['x'], outmap=['x'])
        presult = pcomptr(ado)
        delta11 = self.p11(self.x)
        delta1 = self.p1(self.x)
        xx = self.x + delta1 + delta11
        utils.assert_almost_equal(xx, presult.x)

    def test_inputs_outputs_mismatch(self):
        p2 = models.Poly2DModel(1)
        ch2 = models.Chebyshev2DModel(1, 1)
        with pytest.raises(AssertionError):
            pcomp = PCompositeModel([p2, ch2])


def test_gaussian2d_eval():
    m = models.Gaussian2DModel(2., 3., 4., x_stddev=1., y_stddev=5., theta=30.)
    assert m(3., 4.) == 2.


def test_pickle():
    import copy_reg
    import types
    import cPickle

    def reduce_method(m):
        return (getattr, (m.__self__, m.__func__.__name__))

    copy_reg.pickle(types.MethodType, reduce_method)

    p1 = models.Poly1DModel(3)
    p11 = models.Poly1DModel(4)
    p2 = models.Poly2DModel(3)
    g1 = models.Gaussian1DModel(10.3, 5.4, 1.2)
    scomp_model = SCompositeModel([p1, g1])
    pcomp_model = PCompositeModel([scomp_model, p11])
    s = cPickle.dumps(pcomp_model)
    s1 = cPickle.loads(s)
    assert s1(3) == pcomp_model(3)


@pytest.mark.skipif('not HAS_SCIPY')
def test_powerlaw(scale=5., alpha=2.):
    x = np.linspace(10, 100)
    y = scale * (x) ** (-alpha)
    plm = models.PowerLawModel(1, 1)  # start with a bad guess
    fitter = fitting.NonLinearLSQFitter(plm)
    fitter(x, y)
    assert np.all((fitter.fitpars - np.array([scale, alpha])) < 0.001)

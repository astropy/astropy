# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import numpy as np
from numpy.testing import assert_allclose
from .. import models

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def test_Trapezoid1D():
    """Regression test for
    https://github.com/astropy/astropy/issues/1721
    """
    model = models.Trapezoid1D(amplitude=4.2, x_0=2.0, width=1.0, slope=3)
    xx = np.linspace(0, 4, 8)
    yy = model(xx)
    yy_ref = [0., 1.41428571, 3.12857143, 4.2, 4.2, 3.12857143, 1.41428571, 0.]
    assert_allclose(yy, yy_ref, rtol=0, atol=1e-6)

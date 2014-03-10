# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing import assert_allclose
from .. import models
from astropy.coordinates import Angle

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


def test_Gaussian2D():
    """
    Test rotated elliptical Gaussian2D model.
    https://github.com/astropy/astropy/pull/2038
    """
    model = models.Gaussian2D(100, 1.7, 3.1, x_stddev=3.3, y_stddev=5.0,
                              theta=np.pi/6.)
    y, x = np.mgrid[0:5, 0:5]
    g = model(x, y)
    g_ref = [[77.86787722, 79.8450774, 75.66323816, 66.26262987, 53.6289606],
             [86.01743889, 90.20335997, 87.419011, 78.29536082, 64.80568981],
             [90.11888627, 96.6492347, 95.79172412, 87.74139348, 74.27249818],
             [89.54601432, 98.21442091, 99.55228375, 93.25543692, 80.73169335],
             [84.38744554, 94.65710923, 98.12408063, 94.00369635, 83.22642276]]
    assert_allclose(g, g_ref, rtol=0, atol=1e-6)


def test_Gaussian2DRotation():
    amplitude = 42
    x_mean, y_mean = 0, 0
    x_stddev, y_stddev = 2, 3
    theta = Angle(10, 'deg')
    pars = dict(amplitude=amplitude, x_mean=x_mean, y_mean=y_mean,
                x_stddev=x_stddev, y_stddev=y_stddev)
    rotation_matrix = models.RotateByAngle2D(angle=theta.degree)
    point1 = (x_mean + 2 * x_stddev, y_mean + 2 * y_stddev)
    point2 = rotation_matrix(*point1)
    g1 = models.Gaussian2D(theta=0, **pars)
    g2 = models.Gaussian2D(theta=theta.radian, **pars)
    value1 = g1(*point1)
    value2 = g2(*point2)
    assert_allclose(value1, value2)

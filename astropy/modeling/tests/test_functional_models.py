# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from .. import models
from astropy.coordinates import Angle

try:
    from scipy import optimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def test_Trapezoid1D():
    """Regression test for https://github.com/astropy/astropy/issues/1721"""

    model = models.Trapezoid1D(amplitude=4.2, x_0=2.0, width=1.0, slope=3)
    xx = np.linspace(0, 4, 8)
    yy = model(xx)
    yy_ref = [0., 1.41428571, 3.12857143, 4.2, 4.2, 3.12857143, 1.41428571, 0.]
    assert_allclose(yy, yy_ref, rtol=0, atol=1e-6)


def test_GaussianAbsorption1D():
    g_em = models.Gaussian1D(0.8, 3000, 20)
    g_ab = models.GaussianAbsorption1D(0.8, 3000, 20)
    xx = np.arange(2900, 3100, 2)
    assert_allclose(g_ab(xx), 1 - g_em(xx))
    assert_allclose(g_ab.fit_deriv(xx[0], 0.8, 3000, 20),
                    -np.array(g_em.fit_deriv(xx[0], 0.8, 3000, 20)))
    assert g_ab.bounding_box_default() == g_em.bounding_box


def test_Gaussian2D():
    """
    Test rotated elliptical Gaussian2D model.
    https://github.com/astropy/astropy/pull/2038
    """

    model = models.Gaussian2D(4.2, 1.7, 3.1, x_stddev=5.1, y_stddev=3.3,
                              theta=np.pi/6.)
    y, x = np.mgrid[0:5, 0:5]
    g = model(x, y)
    g_ref = [[3.01907812, 2.99051889, 2.81271552, 2.5119566, 2.13012709],
             [3.55982239, 3.6086023, 3.4734158, 3.17454575, 2.75494838],
             [3.88059142, 4.0257528, 3.96554926, 3.70908389, 3.29410187],
             [3.91095768, 4.15212857, 4.18567526, 4.00652015, 3.64146544],
             [3.6440466, 3.95922417, 4.08454159, 4.00113878, 3.72161094]]
    assert_allclose(g, g_ref, rtol=0, atol=1e-6)


def test_Gaussian2DCovariance():
    """
    Test rotated elliptical Gaussian2D model when cov_matrix is input.
    https://github.com/astropy/astropy/pull/2199
    """

    cov_matrix = [[49., -16.], [-16., 9.]]
    model = models.Gaussian2D(17., 2.0, 2.5, cov_matrix=cov_matrix)
    y, x = np.mgrid[0:5, 0:5]
    g = model(x, y)
    g_ref = [[4.3744505, 5.8413977, 7.42988694, 9.00160175, 10.38794269],
             [8.83290201, 10.81772851, 12.61946384, 14.02225593, 14.84113227],
             [13.68528889, 15.37184621, 16.44637743, 16.76048705, 16.26953638],
             [16.26953638, 16.76048705, 16.44637743, 15.37184621, 13.68528889],
             [14.84113227, 14.02225593, 12.61946384, 10.81772851, 8.83290201]]
    assert_allclose(g, g_ref, rtol=0, atol=1e-6)


def test_Gaussian2DRotation():
    amplitude = 42
    x_mean, y_mean = 0, 0
    x_stddev, y_stddev = 2, 3
    theta = Angle(10, 'deg')
    pars = dict(amplitude=amplitude, x_mean=x_mean, y_mean=y_mean,
                x_stddev=x_stddev, y_stddev=y_stddev)
    rotation = models.Rotation2D(angle=theta.degree)
    point1 = (x_mean + 2 * x_stddev, y_mean + 2 * y_stddev)
    point2 = rotation(*point1)
    g1 = models.Gaussian2D(theta=0, **pars)
    g2 = models.Gaussian2D(theta=theta.radian, **pars)
    value1 = g1(*point1)
    value2 = g2(*point2)
    assert_allclose(value1, value2)


def test_Redshift():
    """Like ``test_ScaleModel()``."""

    # Scale by a scalar
    m = models.Redshift(0.4)
    assert m(0) == 0
    assert_array_equal(m([1, 2]), [1.4, 2.8])

    assert_allclose(m.inverse(m([1, 2])), [1, 2])

    # Scale by a list
    m = models.Redshift([-0.5, 0, 0.5], model_set_axis=0)
    assert_array_equal(m(0), 0)
    assert_array_equal(m([1, 2], model_set_axis=False),
                       [[0.5, 1], [1, 2], [1.5, 3]])

    assert_allclose(m.inverse(m([1, 2], model_set_axis=False)),
                    [[1, 2], [1, 2], [1, 2]])


def test_Ellipse2D():
    """Test Ellipse2D model."""
    amplitude = 7.5
    x0, y0 = 15, 15
    theta = Angle(45, 'deg')
    em = models.Ellipse2D(amplitude, x0, y0, 7, 3, theta.radian)
    y, x = np.mgrid[0:30, 0:30]
    e = em(x, y)
    assert np.all(e[e > 0] == amplitude)
    assert e[y0, x0] == amplitude
    assert em.bounding_box_default() == em.bounding_box

    rotation = models.Rotation2D(angle=theta.degree)
    point1 = [2, 0]      # Rotation2D center is (0, 0)
    point2 = rotation(*point1)
    point1 = np.array(point1) + [x0, y0]
    point2 = np.array(point2) + [x0, y0]
    e1 = models.Ellipse2D(amplitude, x0, y0, 7, 3, theta=0.)
    e2 = models.Ellipse2D(amplitude, x0, y0, 7, 3, theta=theta.radian)
    assert e1(*point1) == e2(*point2)


def test_Ellipse2D_circular():
    """Test that circular Ellipse2D agrees with Disk2D [3736]."""
    amplitude = 7.5
    radius = 10
    size = (radius * 2) + 1
    y, x = np.mgrid[0:size, 0:size]
    ellipse = models.Ellipse2D(amplitude, radius, radius, radius, radius,
                               theta=0)(x, y)
    disk = models.Disk2D(amplitude, radius, radius, radius)(x, y)
    assert np.all(ellipse == disk)


def test_Scale_inverse():
    m = models.Scale(1.2345)
    assert_allclose(m.inverse(m(6.789)), 6.789)


def test_Shift_inverse():
    m = models.Shift(1.2345)
    assert_allclose(m.inverse(m(6.789)), 6.789)

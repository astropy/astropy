from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import pytest

from numpy.testing import assert_allclose

from ..spatial import RipleysKEstimator
from ...utils.misc import NumpyRNGContext


a = np.array([[1, 4], [2, 5], [3, 6]])
b = np.array([[-1, 1], [-2, 2], [-3, 3]])


@pytest.mark.parametrize("points, x_min, x_max", [(a, 0, 10), (b, -5, 5)])
def test_ripley_K_implementation(points, x_min, x_max):
    """
    Test against Ripley's K function implemented in R package `spatstat`
        +-+---------+---------+----------+---------+-+
      6 +                                          * +
        |                                            |
        |                                            |
    5.5 +                                            +
        |                                            |
        |                                            |
      5 +                     *                      +
        |                                            |
    4.5 +                                            +
        |                                            |
        |                                            |
      4 + *                                          +
        +-+---------+---------+----------+---------+-+
          1        1.5        2         2.5        3

        +-+---------+---------+----------+---------+-+
      3 + *                                          +
        |                                            |
        |                                            |
    2.5 +                                            +
        |                                            |
        |                                            |
      2 +                     *                      +
        |                                            |
    1.5 +                                            +
        |                                            |
        |                                            |
      1 +                                          * +
        +-+---------+---------+----------+---------+-+
         -3       -2.5       -2        -1.5       -1
    """

    area = 100
    r = np.linspace(0, 2.5, 5)
    Kest = RipleysKEstimator(area=area, x_min=x_min, y_min=x_min, x_max=x_max,
                             y_max=x_max)

    ANS_NONE = np.array([0, 0, 0, 66.667, 66.667])
    assert_allclose(ANS_NONE, Kest(data=points, radii=r, mode='none'),
                    atol=1e-3)

    ANS_TRANS = np.array([0, 0, 0, 82.304, 82.304])
    assert_allclose(ANS_TRANS, Kest(data=points, radii=r, mode='translation'),
                    atol=1e-3)


with NumpyRNGContext(123):
    a = np.random.uniform(low=5, high=10, size=(100, 2))
    b = np.random.uniform(low=-5, high=-10, size=(100, 2))


@pytest.mark.parametrize("points", [a, b])
def test_ripley_uniform_property(points):
    # Ripley's K function without edge-correction converges to the area when
    # the number of points and the argument radii are large enough, i.e.,
    # K(x) --> area as x --> inf
        area = 50
        Kest = RipleysKEstimator(area=area)
        r = np.linspace(0, 20, 5)
        assert_allclose(area, Kest(data=points, radii=r, mode='none')[4])


with NumpyRNGContext(123):
    a = np.random.uniform(low=0, high=1, size=(500, 2))
    b = np.random.uniform(low=-1, high=0, size=(500, 2))


@pytest.mark.parametrize("points, low, high", [(a, 0, 1), (b, -1, 0)])
def test_ripley_large_density(points, low, high):
        Kest = RipleysKEstimator(area=1, x_min=low, x_max=high, y_min=low,
                                 y_max=high)
        r = np.linspace(0, 0.25, 25)
        Kpos = Kest.poisson(r)
        modes = ['ohser', 'translation', 'ripley']
        for m in modes:
            Kest_r = Kest(data=points, radii=r, mode=m)
            assert_allclose(Kpos, Kest_r, atol=1e-1)


with NumpyRNGContext(123):
    a = np.random.uniform(low=5, high=10, size=(500, 2))
    b = np.random.uniform(low=-10, high=-5, size=(500, 2))


@pytest.mark.parametrize("points, low, high", [(a, 5, 10), (b, -10, -5)])
def test_ripley_modes(points, low, high):
        Kest = RipleysKEstimator(area=25, x_max=high, y_max=high, x_min=low,
                                 y_min=low)
        r = np.linspace(0, 1.2, 25)
        Kpos_mean = np.mean(Kest.poisson(r))
        modes = ['ohser', 'translation', 'ripley']
        for m in modes:
            Kest_mean = np.mean(Kest(data=points, radii=r, mode=m))
            assert_allclose(Kpos_mean, Kest_mean, atol=1e-1, rtol=1e-1)


with NumpyRNGContext(123):
    a = np.random.uniform(low=0, high=1, size=(50, 2))
    b = np.random.uniform(low=-1, high=0, size=(50, 2))


@pytest.mark.parametrize("points, low, high", [(a, 0, 1), (b, -1, 0)])
def test_ripley_large_density_var_width(points, low, high):
        Kest = RipleysKEstimator(area=1, x_min=low, x_max=high, y_min=low,
                                 y_max=high)
        r = np.linspace(0, 0.25, 25)
        Kpos = Kest.poisson(r)
        Kest_r = Kest(data=points, radii=r, mode='var-width')
        assert_allclose(Kpos, Kest_r, atol=1e-1)


with NumpyRNGContext(123):
    a = np.random.uniform(low=5, high=10, size=(50, 2))
    b = np.random.uniform(low=-10, high=-5, size=(50, 2))


@pytest.mark.parametrize("points, low, high", [(a, 5, 10), (b, -10, -5)])
def test_ripley_var_width(points, low, high):
        Kest = RipleysKEstimator(area=25, x_max=high, y_max=high, x_min=low,
                                 y_min=low)
        r = np.linspace(0, 1.2, 25)
        Kest_ohser = np.mean(Kest(data=points, radii=r, mode='ohser'))
        Kest_var_width = np.mean(Kest(data=points, radii=r, mode='var-width'))
        assert_allclose(Kest_ohser, Kest_var_width, atol=1e-1, rtol=1e-1)

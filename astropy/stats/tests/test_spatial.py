from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from numpy.testing import assert_equal
from numpy.testing.utils import assert_allclose

from ..spatial import RipleysKEstimate

def test_ripley_isotropic_implementation():
    # Test against Ripley's K function implemented in R package `spatstat`
    a = np.array([[1, 4], [2, 5], [3, 6]])
    area = 100
    Kest = RipleysKEstimate(data=a, area=area)
    r = np.linspace(0, 2.5, 5)

    answer = np.array([0, 0, 0, 66.667, 66.667])
    assert_allclose(answer, Kest(r), atol=1e-3)

def test_ripley_isotropic_uniform_property():
    # Ripley's K function converges to the area when the number of points
    # and the argument are large enough, i.e., K(x) --> area as x --> inf
    x = np.random.uniform(low=0, high=10, size=100)
    y = np.random.uniform(low=5, high=10, size=100)
    z = np.array([x, y]).T

    area = 50
    Kest = RipleysKEstimate(data=z, area=area)
    r = np.linspace(0, 20, 5)
    assert_allclose(area, Kest(r)[4])

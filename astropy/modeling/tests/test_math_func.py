# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_allclose

from .. import math_functions

x = np.linspace(-20, 360, 100)


@pytest.mark.filterwarnings(
    r'ignore:Models in math_functions are experimental\.')
@pytest.mark.filterwarnings(r'ignore:.*:RuntimeWarning')
def test_math():
    for name in math_functions.__all__:
        model = getattr(math_functions, name)()
        func = getattr(np, model.func.__name__)
        if model.n_inputs == 1:
            assert_allclose(model(x), func(x))
        elif model.n_inputs == 2:
            assert_allclose(model(x, x), func(x, x))

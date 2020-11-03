# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.modeling import math_functions

x = np.linspace(-20, 360, 100)


@pytest.mark.filterwarnings(r'ignore:.*:RuntimeWarning')
def test_math():
    for name in math_functions.__all__:
        model_class = getattr(math_functions, name)
        assert model_class.__module__ == 'astropy.modeling.math_functions'
        model = model_class()
        func = getattr(np, model.func.__name__)
        if model.n_inputs == 1:
            assert_allclose(model(x), func(x))
        elif model.n_inputs == 2:
            assert_allclose(model(x, x), func(x, x))

    assert math_functions.ModUfunc is math_functions.RemainderUfunc
    assert math_functions.DivideUfunc is math_functions.True_divideUfunc

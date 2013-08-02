# Licensed under a 3-clause BSD style license - see LICENSE.rst
import itertools

import numpy as np

from ....tests.helper import pytest

from ..utils import discretize_model
from ....modeling.functional_models import *
from ....modeling.tests.model_lists import models_1D, models_2D
from ....modeling.tests.test_models import create_model

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


modes = ['center', 'corner', 'oversample']
test_models_1D = [Gaussian1DModel, Box1DModel]
test_models_2D = [Gaussian2DModel, Box2DModel]


@pytest.mark.parametrize(('model_class', 'mode'), list(itertools.product(test_models_1D, modes)))
def test_pixel_sum_1D(model_class, mode):
    """
    Test if the sum of all pixels corresponds nearly to the integral.
    """
    parameters = models_1D[model_class]['parameters']
    model = create_model(model_class, parameters)

    values = discretize_model(model, models_1D[model_class]['x_lim'], mode=mode)
    assert (np.abs(values.sum() - models_1D[model_class]['integral']) < 0.0001)


@pytest.mark.parametrize(('model_class', 'mode'), list(itertools.product(test_models_1D, modes)))
def test_eval_1D(model_class, mode):
    """
    Test if the sum of all pixels corresponds nearly to the integral.
    """
    parameters = models_1D[model_class]['parameters']
    model = create_model(model_class, parameters)

    x = np.arange(*models_1D[model_class]['x_lim'])
    values = model(x)
    disc_values = discretize_model(model, models_1D[model_class]['x_lim'], mode=mode)

    assert np.all(np.abs(values - disc_values) < 0.01)

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import pytest

from astropy import units as u
from astropy.modeling.core import Model, fix_inputs
from astropy.modeling.models import Polynomial1D


class _ExampleModel(Model):
    n_inputs = 1
    n_outputs = 1

    def __init__(self):
        self._input_units = {"x": u.m}
        self._return_units = {"y": u.m / u.s}
        super().__init__()

    def evaluate(self, input):
        return input / u.Quantity(1, u.s)


def _models_with_units():
    m1 = _ExampleModel() & _ExampleModel()
    m2 = _ExampleModel() + _ExampleModel()
    p = Polynomial1D(1)
    p._input_units = {"x": u.m / u.s}
    p._return_units = {"y": u.m / u.s}
    m3 = _ExampleModel() | p
    m4 = fix_inputs(m1, {"x0": 1})
    m5 = fix_inputs(m1, {0: 1})

    models = [m1, m2, m3, m4, m5]
    input_units = [
        {"x0": u.Unit("m"), "x1": u.Unit("m")},
        {"x": u.Unit("m")},
        {"x": u.Unit("m")},
        {"x1": u.Unit("m")},
        {"x1": u.Unit("m")},
    ]

    return_units = [
        {"y0": u.Unit("m / s"), "y1": u.Unit("m / s")},
        {"y": u.Unit("m / s")},
        {"y": u.Unit("m / s")},
        {"y0": u.Unit("m / s"), "y1": u.Unit("m / s")},
        {"y0": u.Unit("m / s"), "y1": u.Unit("m / s")},
    ]
    return np.array([models, input_units, return_units], dtype=object).T


@pytest.mark.parametrize(("model", "input_units", "return_units"), _models_with_units())
def test_input_units(model, input_units, return_units):
    """Test input_units on various compound models."""
    assert model.input_units == input_units
    assert model.return_units == return_units

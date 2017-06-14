import pytest

from ..polynomial import Polynomial1D
from ... import units as u
from ...tests.helper import assert_quantity_allclose

MODELS = [
    {
        'class': Polynomial1D,
        'parameters': {'degree': 2, 'c0': 3, 'c1': 2 / u.m , 'c2': 3 / u.m**2},
        'evaluation': [(3 * u.m, 36)]
     },
    {
        'class': Polynomial1D,
        'parameters': {'degree': 2, 'c0': 3 * u.kg, 'c1': 2 * u.kg / u.m , 'c2': 3 * u.kg / u.m**2},
        'evaluation': [(3 * u.m, 36 * u.kg)]
     },
    {
        'class': Polynomial1D,
        'parameters': {'degree': 2, 'c0': 3 * u.kg, 'c1': 2 * u.kg, 'c2': 3 * u.kg},
        'evaluation': [(3, 36 * u.kg)]
     }
 ]


@pytest.mark.parametrize('model', MODELS)
def test_polynomial_eval_withunits(model):
    m = model['class'](**model['parameters'])
    for x, y in model['evaluation']:
        assert_quantity_allclose(m(x), y)

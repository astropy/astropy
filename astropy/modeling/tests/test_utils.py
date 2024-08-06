# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name
import re
from inspect import Parameter

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.modeling import models
from astropy.modeling.utils import (
    _SpecialOperatorsDict,
    _validate_domain_window,
    copy_with_new_parameter_values,
    get_inputs_and_params,
    poly_map_domain,
)


def test_poly_map_domain():
    oldx = np.array([1, 2, 3, 4])

    # test shift/scale
    assert (poly_map_domain(oldx, (-4, 4), (-3, 3)) == [0.75, 1.5, 2.25, 3]).all()

    # errors
    MESSAGE = r'Expected "domain" and "window" to be a tuple of size 2'
    with pytest.raises(ValueError, match=MESSAGE):
        poly_map_domain(oldx, (-4,), (-3, 3))
    with pytest.raises(ValueError, match=MESSAGE):
        poly_map_domain(oldx, (-4, 4, -4), (-3, 3))
    with pytest.raises(ValueError, match=MESSAGE):
        poly_map_domain(oldx, (-4, 4), (-3,))
    with pytest.raises(ValueError, match=MESSAGE):
        poly_map_domain(oldx, (-4, 4), (-3, 3, -3))


def test__validate_domain_window():
    # Test if None
    assert _validate_domain_window(None) is None

    # Test normal
    assert _validate_domain_window((-2, 2)) == (-2, 2)
    assert _validate_domain_window([-2, 2]) == (-2, 2)
    assert _validate_domain_window(np.array([-2, 2])) == (-2, 2)

    # Test error
    MESSAGE = r"domain and window should be tuples of size 2"
    with pytest.raises(ValueError, match=MESSAGE):
        _validate_domain_window((-2, 2, -2))
    with pytest.raises(ValueError, match=MESSAGE):
        _validate_domain_window((-2,))
    with pytest.raises(ValueError, match=MESSAGE):
        _validate_domain_window([-2])
    with pytest.raises(ValueError, match=MESSAGE):
        _validate_domain_window(np.array([-2]))
    with pytest.raises(ValueError, match=MESSAGE):
        _validate_domain_window(-2)


def test_get_inputs_and_params():
    # test normal
    def func1(input0, input1, param0=5, param1=7):
        pass

    inputs, params = get_inputs_and_params(func1)
    for index, _input in enumerate(inputs):
        assert isinstance(_input, Parameter)
        assert _input.name == f"input{index}"
        assert _input.kind == _input.POSITIONAL_OR_KEYWORD
        assert _input.default == Parameter.empty
    default = [5, 7]
    for index, param in enumerate(params):
        assert isinstance(param, Parameter)
        assert param.name == f"param{index}"
        assert param.kind == param.POSITIONAL_OR_KEYWORD
        assert param.default == default[index]

    # Error
    MESSAGE = re.escape("Signature must not have *args or **kwargs")

    def func2(input0, input1, *args, param0=5, param1=7):
        pass

    def func3(input0, input1, param0=5, param1=7, **kwargs):
        pass

    with pytest.raises(ValueError, match=MESSAGE):
        get_inputs_and_params(func2)
    with pytest.raises(ValueError, match=MESSAGE):
        get_inputs_and_params(func3)


class Test_SpecialOperatorsDict:
    def setup_method(self):
        self.key = "test"
        self.val = "value"

    def test__set_value(self):
        special_operators = _SpecialOperatorsDict()
        assert self.key not in special_operators

        special_operators._set_value(self.key, self.val)
        assert self.key in special_operators
        assert special_operators[self.key] == self.val

        with pytest.raises(ValueError, match='Special operator "test" already exists'):
            special_operators._set_value(self.key, self.val)

    def test___setitem__(self):
        special_operators = _SpecialOperatorsDict()
        assert self.key not in special_operators

        with pytest.deprecated_call():
            special_operators[self.key] = self.val
        assert self.key in special_operators
        assert special_operators[self.key] == self.val

    def test__SpecialOperatorsDict__get_unique_id(self):
        special_operators = _SpecialOperatorsDict()
        assert special_operators._unique_id == 0

        assert special_operators._get_unique_id() == 1
        assert special_operators._unique_id == 1

        assert special_operators._get_unique_id() == 2
        assert special_operators._unique_id == 2

        assert special_operators._get_unique_id() == 3
        assert special_operators._unique_id == 3

    def test__SpecialOperatorsDict_add(self):
        special_operators = _SpecialOperatorsDict()

        operator_name = "test"
        operator = "operator"

        key0 = special_operators.add(operator_name, operator)
        assert key0 == (operator_name, special_operators._unique_id)
        assert key0 in special_operators
        assert special_operators[key0] == operator

        key1 = special_operators.add(operator_name, operator)
        assert key1 == (operator_name, special_operators._unique_id)
        assert key1 in special_operators
        assert special_operators[key1] == operator

        assert key0 != key1


from .test_models_quantities import MODELS

# For now reuse a list of models from another test file, but we should make
# sure we use a complete list of all models.
ALL_MODELS = [entry["class"](**entry["parameters"]) for entry in MODELS]

# Add a compound model
ALL_MODELS.append(models.Gaussian1D() + models.Polynomial1D(2))


@pytest.mark.parametrize("model", ALL_MODELS, ids=lambda m: m.__class__.__name__)
def test_copy_with_new_parameter_values(model):
    # Extract first parameter name
    first_parameter = model.param_names[0]

    getattr(model, first_parameter).fixed = True

    # Set up an array which will be used for the new parameter value
    array = np.ones((2, 3, 4))

    # Copy the model, replacing the first parameter with the array
    new_model = copy_with_new_parameter_values(model, **{first_parameter: array})

    # Check the result
    assert_allclose(getattr(new_model, first_parameter), array)
    assert new_model.fixed[first_parameter]

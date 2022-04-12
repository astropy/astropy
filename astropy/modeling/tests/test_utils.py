# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name
from inspect import Parameter

import numpy as np
import pytest

from astropy.modeling.utils import (
    _SpecialOperatorsDict, _validate_domain_window, get_inputs_and_params, poly_map_domain)


def test_poly_map_domain():
    oldx = np.array([1, 2, 3, 4])

    # test shift/scale
    assert (poly_map_domain(oldx, (-4, 4), (-3, 3)) == [0.75, 1.5, 2.25, 3]).all()

    # errors
    MESSAGE = 'Expected "domain" and "window" to be a tuple of size 2.'
    with pytest.raises(ValueError) as err:
        poly_map_domain(oldx, (-4,), (-3, 3))
    assert str(err.value) == MESSAGE
    with pytest.raises(ValueError) as err:
        poly_map_domain(oldx, (-4, 4, -4), (-3, 3))
    assert str(err.value) == MESSAGE
    with pytest.raises(ValueError) as err:
        poly_map_domain(oldx, (-4, 4), (-3,))
    assert str(err.value) == MESSAGE
    with pytest.raises(ValueError) as err:
        poly_map_domain(oldx, (-4, 4), (-3, 3, -3))
    assert str(err.value) == MESSAGE


def test__validate_domain_window():
    # Test if None
    assert _validate_domain_window(None) is None

    # Test normal
    assert _validate_domain_window((-2, 2)) == (-2, 2)
    assert _validate_domain_window([-2, 2]) == (-2, 2)
    assert _validate_domain_window(np.array([-2, 2])) == (-2, 2)

    # Test error
    MESSAGE = 'domain and window should be tuples of size 2.'
    with pytest.raises(ValueError) as err:
        _validate_domain_window((-2, 2, -2))
    assert str(err.value) == MESSAGE
    with pytest.raises(ValueError) as err:
        _validate_domain_window((-2,))
    assert str(err.value) == MESSAGE
    with pytest.raises(ValueError) as err:
        _validate_domain_window([-2])
    assert str(err.value) == MESSAGE
    with pytest.raises(ValueError) as err:
        _validate_domain_window(np.array([-2]))
    assert str(err.value) == MESSAGE
    with pytest.raises(ValueError) as err:
        _validate_domain_window(-2)
    assert str(err.value) == MESSAGE


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
    MESSAGE = "Signature must not have *args or **kwargs"

    def func2(input0, input1, *args, param0=5, param1=7):
        pass

    def func3(input0, input1, param0=5, param1=7, **kwargs):
        pass

    with pytest.raises(ValueError) as err:
        get_inputs_and_params(func2)
    assert str(err.value) == MESSAGE
    with pytest.raises(ValueError) as err:
        get_inputs_and_params(func3)
    assert str(err.value) == MESSAGE


class Test_SpecialOperatorsDict:
    def setup(self):
        self.key = 'test'
        self.val = 'value'

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

        operator_name = 'test'
        operator = 'operator'

        key0 = special_operators.add(operator_name, operator)
        assert key0 == (operator_name, special_operators._unique_id)
        assert key0 in special_operators
        assert special_operators[key0] == operator

        key1 = special_operators.add(operator_name, operator)
        assert key1 == (operator_name, special_operators._unique_id)
        assert key1 in special_operators
        assert special_operators[key1] == operator

        assert key0 != key1

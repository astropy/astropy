# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name
import pytest

from astropy.modeling.utils import _SpecialOperatorsDict


def test__SpecialOperatorsDict__set_value():
    key = 'test'
    val = 'value'

    special_operators = _SpecialOperatorsDict()
    assert key not in special_operators

    special_operators._set_value(key, val)
    assert key in special_operators
    assert special_operators[key] == val

    with pytest.raises(ValueError, match='Special operator "test" already exists'):
        special_operators._set_value(key, val)


def test__SpecialOperatorsDict___setitem__():
    key = 'test'
    val = 'value'

    special_operators = _SpecialOperatorsDict()
    assert key not in special_operators

    with pytest.deprecated_call():
        special_operators[key] = val
    assert key in special_operators
    assert special_operators[key] == val


def test__SpecialOperatorsDict__get_unique_id():
    special_operators = _SpecialOperatorsDict()
    assert special_operators._unique_id == 0

    assert special_operators._get_unique_id() == 1
    assert special_operators._unique_id == 1

    assert special_operators._get_unique_id() == 2
    assert special_operators._unique_id == 2

    assert special_operators._get_unique_id() == 3
    assert special_operators._unique_id == 3


def test__SpecialOperatorsDict_add():
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

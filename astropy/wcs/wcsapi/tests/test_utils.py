from pytest import raises

from ....tests.helper import assert_quantity_allclose
from .... import units as u

from ..utils import deserialize_class


def test_construct():

    result = deserialize_class(('astropy.units.Quantity', (10,), {'unit': 'deg'}))
    assert_quantity_allclose(result, 10 * u.deg)


def test_noconstruct():

    result = deserialize_class(('astropy.units.Quantity', (), {'unit': 'deg'}), construct=False)
    assert result == (u.Quantity, (), {'unit': 'deg'})


def test_invalid():

    with raises(ValueError) as exc:
        deserialize_class(('astropy.units.Quantity', (), {'unit': 'deg'}, ()))
    assert exc.value.args[0] == 'Expected a tuple of three values'

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest
import numpy as np

from astropy import units as u
from astropy.units import Quantity, UnitsError, equivalencies
from astropy.modeling.models import UnitsMapping


def test_properties():
    model = UnitsMapping(((u.dimensionless_unscaled, u.m), (u.dimensionless_unscaled, u.s)))
    assert model.n_inputs == 2
    assert model.n_outputs == 2
    assert model.inputs == ("x0", "x1")
    assert model.outputs == ("x0", "x1")
    assert model.input_units == {"x0": u.dimensionless_unscaled, "x1": u.dimensionless_unscaled}
    assert model.mapping == ((u.dimensionless_unscaled, u.m), (u.dimensionless_unscaled, u.s))


def test_add_units():
    model = UnitsMapping(((u.dimensionless_unscaled, u.m),))

    for value in [10, Quantity(10), np.arange(10), Quantity(np.arange(10))]:
        result = model(value)
        assert isinstance(result, Quantity)
        assert np.all(result.value == value)
        assert result.unit == u.m

    with pytest.raises(UnitsError):
        model(Quantity(10, u.s))


def test_remove_units():
    model = UnitsMapping(((u.m, u.dimensionless_unscaled),))

    result = model(Quantity(10, u.m))
    assert isinstance(result, Quantity)
    assert result.value == 10
    assert result.unit == u.dimensionless_unscaled

    result = model(Quantity(1000, u.cm))
    assert isinstance(result, Quantity)
    assert result.value == 10
    assert result.unit == u.dimensionless_unscaled

    with pytest.raises(UnitsError):
        model(10)

    with pytest.raises(UnitsError):
        model(Quantity(10))


def test_remove_quantity():
    model = UnitsMapping(((u.m, None),))

    result = model(Quantity(10, u.m))
    assert result == 10

    result = model(Quantity(1000, u.cm))
    assert result == 10

    with pytest.raises(UnitsError):
        model(10)

    with pytest.raises(UnitsError):
        model(Quantity(10))

    # The model shouldn't allow a mixture of None and non-None
    # output units.
    with pytest.raises(ValueError, match=r"If one return unit is None, then all must be None"):
        UnitsMapping(((u.m, None), (u.s, u.dimensionless_unscaled)))


def test_equivalencies():
    model = UnitsMapping(((u.m, u.dimensionless_unscaled),))

    with pytest.raises(UnitsError):
        model(Quantity(100, u.Hz))

    model = UnitsMapping(((u.m, u.dimensionless_unscaled),), input_units_equivalencies={"x": equivalencies.spectral()})

    result = model(Quantity(100, u.Hz))
    assert result.unit == u.dimensionless_unscaled


def test_allow_dimensionless():
    model = UnitsMapping(((u.m, u.dimensionless_unscaled),))

    with pytest.raises(UnitsError):
        model(10)

    model = UnitsMapping(((u.m, u.dimensionless_unscaled),), input_units_allow_dimensionless=True)
    result = model(10)
    assert isinstance(result, Quantity)
    assert result.value == 10
    assert result.unit == u.dimensionless_unscaled


def test_custom_inputs_and_outputs():
    model = UnitsMapping(((u.m, u.dimensionless_unscaled),))

    model.inputs = ("foo",)
    model.outputs = ("bar",)

    assert model.inputs == ("foo",)
    assert model.input_units == {"foo": u.m}
    assert model.outputs == ("bar",)


def test_repr():
    model = UnitsMapping(((u.m, None),))
    assert repr(model) == """<UnitsMapping(((Unit("m"), None),))>"""

    model = UnitsMapping(((u.m, None),), name="foo")
    assert repr(model) == """<UnitsMapping(((Unit("m"), None),), name='foo')>"""

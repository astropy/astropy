# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from ... import units as u

from .py3_test_quantity_annotations import *

# list of pairs (target unit/physical type, input unit)
x_inputs = [(u.arcsec, u.deg), ('angle', u.deg),
            (u.kpc/u.Myr, u.km/u.s), ('speed', u.km/u.s),
            ([u.arcsec, u.km], u.deg), ([u.arcsec, u.km], u.km),  # multiple allowed
            (['angle', 'length'], u.deg), (['angle', 'length'], u.km)]

y_inputs = [(u.arcsec, u.deg), ('angle', u.deg),
            (u.kpc/u.Myr, u.km/u.s), ('speed', u.km/u.s)]


@pytest.fixture(scope="module",
                params=list(range(len(x_inputs))))
def x_input(request):
    return x_inputs[request.param]


@pytest.fixture(scope="module",
                params=list(range(len(y_inputs))))
def y_input(request):
    return y_inputs[request.param]

# ---- Tests that use the fixtures defined above ----


def test_args(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, y):
        return x, y

    x, y = myfunc_args(1*x_unit, 1*y_unit)

    assert isinstance(x, u.Quantity)
    assert isinstance(y, u.Quantity)

    assert x.unit == x_unit
    assert y.unit == y_unit


def test_args_nonquantity(x_input):
    x_target, x_unit = x_input

    @u.quantity_input(x=x_target)
    def myfunc_args(x, y):
        return x, y

    x, y = myfunc_args(1*x_unit, 100)

    assert isinstance(x, u.Quantity)
    assert isinstance(y, int)

    assert x.unit == x_unit


def test_wrong_unit(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, y):
        return x, y

    with pytest.raises(u.UnitsError) as e:
        x, y = myfunc_args(1*x_unit, 100*u.Joule)  # has to be an unspecified unit

    str_to = str(y_target)
    assert str(e.value) == "Argument 'y' to function 'myfunc_args' must be in units convertible to '{0}'.".format(str_to)


def test_not_quantity(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, y):
        return x, y

    with pytest.raises(TypeError) as e:
        x, y = myfunc_args(1*x_unit, 100)
    assert str(e.value) == "Argument 'y' to function 'myfunc_args' has no 'unit' attribute. You may want to pass in an astropy Quantity instead."


def test_kwargs(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, my_arg, y=1*y_unit):
        return x, my_arg, y

    x, my_arg, y = myfunc_args(1*x_unit, 100, y=100*y_unit)

    assert isinstance(x, u.Quantity)
    assert isinstance(my_arg, int)
    assert isinstance(y, u.Quantity)

    assert y.unit == y_unit


def test_unused_kwargs(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, my_arg1, y=y_unit, my_arg2=1000):
        return x, my_arg1, y, my_arg2

    x, my_arg1, y, my_arg2 = myfunc_args(1*x_unit, 100,
                                         y=100*y_unit, my_arg2=10)

    assert isinstance(x, u.Quantity)
    assert isinstance(my_arg1, int)
    assert isinstance(y, u.Quantity)
    assert isinstance(my_arg2, int)

    assert y.unit == y_unit
    assert my_arg2 == 10


def test_kwarg_wrong_unit(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, y=10*y_unit):
        return x, y

    with pytest.raises(u.UnitsError) as e:
        x, y = myfunc_args(1*x_unit, y=100*u.Joule)

    str_to = str(y_target)
    assert str(e.value) == "Argument 'y' to function 'myfunc_args' must be in units convertible to '{0}'.".format(str_to)


def test_kwarg_not_quantity(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, y=10*y_unit):
        return x, y

    with pytest.raises(TypeError) as e:
        x, y = myfunc_args(1*x_unit, y=100)
    assert str(e.value) == "Argument 'y' to function 'myfunc_args' has no 'unit' attribute. You may want to pass in an astropy Quantity instead."


def test_kwarg_default(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, y=10*y_unit):
        return x, y

    x, y = myfunc_args(1*x_unit)

    assert isinstance(x, u.Quantity)
    assert isinstance(y, u.Quantity)

    assert x.unit == x_unit
    assert y.unit == y_unit


def test_kwargs_input(x_input, y_input):
    x_target, x_unit = x_input
    y_target, y_unit = y_input

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x=1*x_unit, y=1*y_unit):
        return x, y

    kwargs = {'x': 10*x_unit, 'y': 10*y_unit}
    x, y = myfunc_args(**kwargs)

    assert isinstance(x, u.Quantity)
    assert isinstance(y, u.Quantity)

    assert x.unit == x_unit
    assert y.unit == y_unit


def test_kwargs_extra(x_input):
    x_target, x_unit = x_input

    @u.quantity_input(x=x_target)
    def myfunc_args(x, **kwargs):
        return x

    x = myfunc_args(1*x_unit)

    assert isinstance(x, u.Quantity)

    assert x.unit == x_unit

# ---- Tests that don't used the fixtures ----


@pytest.mark.parametrize("x_unit,y_unit", [
                         (u.arcsec, u.eV),
                         ('angle', 'energy')])
def test_arg_equivalencies(x_unit, y_unit):
    @u.quantity_input(x=x_unit, y=y_unit,
                      equivalencies=u.mass_energy())
    def myfunc_args(x, y):
        return x, y+(10*u.J)  # Add an energy to check equiv is working

    x, y = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(x, u.Quantity)
    assert isinstance(y, u.Quantity)

    assert x.unit == u.arcsec
    assert y.unit == u.gram


@pytest.mark.parametrize("x_unit,energy_unit", [
                         (u.arcsec, u.eV),
                         ('angle', 'energy')])
def test_kwarg_equivalencies(x_unit, energy_unit):
    @u.quantity_input(x=x_unit, energy=energy_unit, equivalencies=u.mass_energy())
    def myfunc_args(x, energy=10*u.eV):
        return x, energy+(10*u.J)  # Add an energy to check equiv is working

    x, energy = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(x, u.Quantity)
    assert isinstance(energy, u.Quantity)

    assert x.unit == u.arcsec
    assert energy.unit == u.gram


def test_no_equivalent():
    class test_unit(object):
        pass

    class test_quantity(object):
        unit = test_unit()

    @u.quantity_input(x=u.arcsec)
    def myfunc_args(x):
        return x

    with pytest.raises(TypeError) as e:
        x, y = myfunc_args(test_quantity())

        assert str(e.value) == "Argument 'x' to function 'myfunc_args' has a 'unit' attribute without an 'is_equivalent' method. You may want to pass in an astropy Quantity instead."


def test_kwarg_invalid_physical_type():
    @u.quantity_input(x='angle', y='africanswallow')
    def myfunc_args(x, y=10*u.deg):
        return x, y

    with pytest.raises(ValueError) as e:
        x, y = myfunc_args(1*u.arcsec, y=100*u.deg)
    assert str(e.value) == "Invalid unit or physical type 'africanswallow'."


def test_default_value_check():
    x_target = u.deg
    x_unit = u.arcsec

    with pytest.raises(TypeError):
        @u.quantity_input(x=x_target)
        def myfunc_args(x=1.):
            return x

        x = myfunc_args()

    x = myfunc_args(1*x_unit)
    assert isinstance(x, u.Quantity)
    assert x.unit == x_unit


def test_args_None():
    x_target = u.deg
    x_unit = u.arcsec
    y_target = u.km
    y_unit = u.kpc

    @u.quantity_input(x=[x_target, None], y=[None, y_target])
    def myfunc_args(x, y):
        return x, y

    x, y = myfunc_args(1*x_unit, None)
    assert isinstance(x, u.Quantity)
    assert x.unit == x_unit
    assert y is None

    x, y = myfunc_args(None, 1*y_unit)
    assert isinstance(y, u.Quantity)
    assert y.unit == y_unit
    assert x is None


def test_args_None_kwarg():
    x_target = u.deg
    x_unit = u.arcsec
    y_target = u.km

    @u.quantity_input(x=x_target, y=y_target)
    def myfunc_args(x, y=None):
        return x, y

    x, y = myfunc_args(1*x_unit)
    assert isinstance(x, u.Quantity)
    assert x.unit == x_unit
    assert y is None

    x, y = myfunc_args(1*x_unit, None)
    assert isinstance(x, u.Quantity)
    assert x.unit == x_unit
    assert y is None

    with pytest.raises(TypeError):
        x, y = myfunc_args(None, None)

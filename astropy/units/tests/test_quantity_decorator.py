# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from ... import units as u


from .py3_test_quantity_annotations import *


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.arcsec),
                         ('angle', 'angle')])
def test_args(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit)
    def myfunc_args(solarx, solary):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 1*u.arcsec)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.arcsec

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.arcsec),
                         ('angle', 'angle')])
def test_args_noconvert(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit)
    def myfunc_args(solarx, solary):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.deg, 1*u.arcmin)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.deg
    assert solary.unit == u.arcmin


@pytest.mark.parametrize("solarx_unit", [
                         u.arcsec, 'angle'])
def test_args_nonquantity(solarx_unit):
    @u.quantity_input(solarx=solarx_unit)
    def myfunc_args(solarx, solary):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 100)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)

    assert solarx.unit == u.arcsec

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.eV),
                         ('angle', 'energy')])
def test_arg_equivalencies(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit, equivalencies=u.mass_energy())
    def myfunc_args(solarx, solary):
        return solarx, solary+(10*u.J)  # Add an energy to check equiv is working

    solarx, solary = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.gram

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_wrong_unit(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit)
    def myfunc_args(solarx, solary):
        return solarx, solary

    with pytest.raises(u.UnitsError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, 100*u.km)

    str_to = str(solary_unit)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' must be in units convertible to '{0}'.".format(str_to)

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_not_quantity(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit)
    def myfunc_args(solarx, solary):
        return solarx, solary

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, 100)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' has no 'unit' attribute. You may want to pass in an astropy Quantity instead."

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwargs(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, myk=solary_unit)
    def myfunc_args(solarx, solary, myk=1*u.arcsec):
        return solarx, solary, myk

    solarx, solary, myk = myfunc_args(1*u.arcsec, 100, myk=100*u.deg)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)
    assert isinstance(myk, u.Quantity)

    assert myk.unit == u.deg

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_unused_kwargs(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, myk=solary_unit)
    def myfunc_args(solarx, solary, myk=1*u.arcsec, myk2=1000):
        return solarx, solary, myk, myk2

    solarx, solary, myk, myk2 = myfunc_args(1*u.arcsec, 100, myk=100*u.deg, myk2=10)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)
    assert isinstance(myk, u.Quantity)
    assert isinstance(myk2, int)

    assert myk.unit == u.deg
    assert myk2 == 10

@pytest.mark.parametrize("solarx_unit,energy_unit", [
                         (u.arcsec, u.eV),
                         ('angle', 'energy')])
def test_kwarg_equivalencies(solarx_unit, energy_unit):
    @u.quantity_input(solarx=solarx_unit, energy=energy_unit, equivalencies=u.mass_energy())
    def myfunc_args(solarx, energy=10*u.eV):
        return solarx, energy+(10*u.J)  # Add an energy to check equiv is working

    solarx, energy = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(energy, u.Quantity)

    assert solarx.unit == u.arcsec
    assert energy.unit == u.gram

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwarg_wrong_unit(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit)
    def myfunc_args(solarx, solary=10*u.deg):
        return solarx, solary

    with pytest.raises(u.UnitsError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100*u.km)

    str_to = str(solary_unit)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' must be in units convertible to '{0}'.".format(str_to)

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwarg_not_quantity(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit)
    def myfunc_args(solarx, solary=10*u.deg):
        return solarx, solary

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' has no 'unit' attribute. You may want to pass in an astropy Quantity instead."

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwarg_default(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit)
    def myfunc_args(solarx, solary=10*u.deg):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.deg

def test_no_equivalent():
    class test_unit(object):
        pass
    class test_quantity(object):
        unit = test_unit()

    @u.quantity_input(solarx=u.arcsec)
    def myfunc_args(solarx):
        return solarx

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(test_quantity())

        assert str(e.value) == "Argument 'solarx' to function 'myfunc_args' has a 'unit' attribute without an 'is_equivalent' method. You may want to pass in an astropy Quantity instead."

@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwargs_input(solarx_unit, solary_unit):
    @u.quantity_input(solarx=solarx_unit, solary=solary_unit)
    def myfunc_args(solarx=1*u.arcsec, solary=1*u.deg):
        return solarx, solary

    kwargs = {'solarx':10*u.arcsec, 'solary':10*u.deg}
    solarx, solary = myfunc_args(**kwargs)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.deg

@pytest.mark.parametrize("solarx_unit", [
                         u.arcsec, 'angle'])
def test_kwargs_extra(solarx_unit):
    @u.quantity_input(solarx=solarx_unit)
    def myfunc_args(solarx, **kwargs):
        return solarx

    solarx = myfunc_args(1*u.deg)

    assert isinstance(solarx, u.Quantity)

    assert solarx.unit == u.deg

def test_kwarg_invalid_physical_type():
    @u.quantity_input(solarx='angle', solary='africanswallow')
    def myfunc_args(solarx, solary=10*u.deg):
        return solarx, solary

    with pytest.raises(ValueError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100*u.deg)
    assert str(e.value) == "Invalid unit of physical type 'africanswallow'."

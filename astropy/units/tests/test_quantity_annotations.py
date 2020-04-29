# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy import units as u  # pylint: disable=W0611


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.arcsec),
                         ('angle', 'angle')])
def test_args3(solarx_unit, solary_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary: solary_unit):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 1*u.arcsec)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.arcsec


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.arcsec),
                         ('angle', 'angle')])
def test_args_noconvert3(solarx_unit, solary_unit):
    @u.quantity_input()
    def myfunc_args(solarx: solarx_unit, solary: solary_unit):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.deg, 1*u.arcmin)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.deg
    assert solary.unit == u.arcmin


@pytest.mark.parametrize("solarx_unit", [
                         u.arcsec, 'angle'])
def test_args_nonquantity3(solarx_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 100)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)

    assert solarx.unit == u.arcsec


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.eV),
                         ('angle', 'energy')])
def test_arg_equivalencies3(solarx_unit, solary_unit):
    @u.quantity_input(equivalencies=u.mass_energy())
    def myfunc_args(solarx: solarx_unit, solary: solary_unit):
        return solarx, solary+(10*u.J)  # Add an energy to check equiv is working

    solarx, solary = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.gram


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_wrong_unit3(solarx_unit, solary_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary: solary_unit):
        return solarx, solary

    with pytest.raises(u.UnitsError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, 100*u.km)

    str_to = str(solary_unit)
    assert str(e.value) == f"Argument 'solary' to function 'myfunc_args' must be in units convertible to '{str_to}'."


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_not_quantity3(solarx_unit, solary_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary: solary_unit):
        return solarx, solary

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, 100)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' has no 'unit' attribute. You should pass in an astropy Quantity instead."


def test_decorator_override():
    @u.quantity_input(solarx=u.arcsec)
    def myfunc_args(solarx: u.km, solary: u.arcsec):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 1*u.arcsec)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.arcsec


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwargs3(solarx_unit, solary_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary, myk: solary_unit=1*u.arcsec):
        return solarx, solary, myk

    solarx, solary, myk = myfunc_args(1*u.arcsec, 100, myk=100*u.deg)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)
    assert isinstance(myk, u.Quantity)

    assert myk.unit == u.deg


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_unused_kwargs3(solarx_unit, solary_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary, myk: solary_unit=1*u.arcsec, myk2=1000):
        return solarx, solary, myk, myk2

    solarx, solary, myk, myk2 = myfunc_args(1*u.arcsec, 100, myk=100*u.deg, myk2=10)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)
    assert isinstance(myk, u.Quantity)
    assert isinstance(myk2, int)

    assert myk.unit == u.deg
    assert myk2 == 10


@pytest.mark.parametrize("solarx_unit,energy", [
                         (u.arcsec, u.eV),
                         ('angle', 'energy')])
def test_kwarg_equivalencies3(solarx_unit, energy):
    @u.quantity_input(equivalencies=u.mass_energy())
    def myfunc_args(solarx: solarx_unit, energy: energy=10*u.eV):
        return solarx, energy+(10*u.J)  # Add an energy to check equiv is working

    solarx, energy = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(energy, u.Quantity)

    assert solarx.unit == u.arcsec
    assert energy.unit == u.gram


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwarg_wrong_unit3(solarx_unit, solary_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary: solary_unit=10*u.deg):
        return solarx, solary

    with pytest.raises(u.UnitsError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100*u.km)

    str_to = str(solary_unit)
    assert str(e.value) == f"Argument 'solary' to function 'myfunc_args' must be in units convertible to '{str_to}'."


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwarg_not_quantity3(solarx_unit, solary_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary: solary_unit=10*u.deg):
        return solarx, solary

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' has no 'unit' attribute. You should pass in an astropy Quantity instead."


@pytest.mark.parametrize("solarx_unit,solary_unit", [
                         (u.arcsec, u.deg),
                         ('angle', 'angle')])
def test_kwarg_default3(solarx_unit, solary_unit):
    @u.quantity_input
    def myfunc_args(solarx: solarx_unit, solary: solary_unit=10*u.deg):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec)


def test_return_annotation():
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec) -> u.deg:
        return solarx

    solarx = myfunc_args(1*u.arcsec)
    assert solarx.unit is u.deg


def test_return_annotation_none():
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec) -> None:
        pass

    solarx = myfunc_args(1*u.arcsec)
    assert solarx is None


def test_enum_annotation():
    # Regression test for gh-9932
    from enum import Enum, auto

    class BasicEnum(Enum):
        AnOption = auto()

    @u.quantity_input
    def myfunc_args(a: BasicEnum, b: u.arcsec) -> None:
        pass

    myfunc_args(BasicEnum.AnOption, 1*u.arcsec)

# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.extern import six

from ... import units as u
from ...extern import six
from ...tests.helper import pytest


from .py3_test_quantity_annotations import *


def test_args():
    @u.quantity_input(solarx=u.arcsec, solary=u.arcsec)
    def myfunc_args(solarx, solary):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 1*u.arcsec)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.arcsec

def test_args_noconvert():
    @u.quantity_input(solarx=u.arcsec, solary=u.arcsec)
    def myfunc_args(solarx, solary):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.deg, 1*u.arcmin)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.deg
    assert solary.unit == u.arcmin


def test_args_nonquantity():
    @u.quantity_input(solarx=u.arcsec)
    def myfunc_args(solarx, solary):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 100)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)

    assert solarx.unit == u.arcsec

def test_arg_equivalencies():
    @u.quantity_input(solarx=u.arcsec, solary=u.eV, equivalencies=u.mass_energy())
    def myfunc_args(solarx, solary):
        return solarx, solary+(10*u.J)  # Add an energy to check equiv is working

    solarx, solary = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.gram

def test_wrong_unit():
    @u.quantity_input(solarx=u.arcsec, solary=u.deg)
    def myfunc_args(solarx, solary):
        return solarx, solary

    with pytest.raises(u.UnitsError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, 100*u.km)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' must be in units convertable to 'deg'."

def test_not_quantity():
    @u.quantity_input(solarx=u.arcsec, solary=u.deg)
    def myfunc_args(solarx, solary):
        return solarx, solary

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, 100)
    assert str(e.value) == "Argument 'solary' to function has 'myfunc_args' no 'unit' attribute. You may want to pass in an astropy Quantity instead."

def test_kwargs():
    @u.quantity_input(solarx=u.arcsec, myk=u.deg)
    def myfunc_args(solarx, solary, myk=1*u.arcsec):
        return solarx, solary, myk

    solarx, solary, myk = myfunc_args(1*u.arcsec, 100, myk=100*u.deg)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)
    assert isinstance(myk, u.Quantity)

    assert myk.unit == u.deg

def test_unused_kwargs():
    @u.quantity_input(solarx=u.arcsec, myk=u.deg)
    def myfunc_args(solarx, solary, myk=1*u.arcsec, myk2=1000):
        return solarx, solary, myk, myk2

    solarx, solary, myk, myk2 = myfunc_args(1*u.arcsec, 100, myk=100*u.deg, myk2=10)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)
    assert isinstance(myk, u.Quantity)
    assert isinstance(myk2, int)

    assert myk.unit == u.deg
    assert myk2 == 10

def test_kwarg_equivalencies():
    @u.quantity_input(solarx=u.arcsec, energy=u.eV, equivalencies=u.mass_energy())
    def myfunc_args(solarx, energy=10*u.eV):
        return solarx, energy+(10*u.J)  # Add an energy to check equiv is working

    solarx, energy = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(energy, u.Quantity)

    assert solarx.unit == u.arcsec
    assert energy.unit == u.gram

def test_kwarg_wrong_unit():
    @u.quantity_input(solarx=u.arcsec, solary=u.deg)
    def myfunc_args(solarx, solary=10*u.deg):
        return solarx, solary

    with pytest.raises(u.UnitsError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100*u.km)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' must be in units convertable to 'deg'."

def test_kwarg_not_quantity():
    @u.quantity_input(solarx=u.arcsec, solary=u.deg)
    def myfunc_args(solarx, solary=10*u.deg):
        return solarx, solary

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100)
    assert str(e.value) == "Argument 'solary' to function has 'myfunc_args' no 'unit' attribute. You may want to pass in an astropy Quantity instead."

def test_kwarg_default():
    @u.quantity_input(solarx=u.arcsec, solary=u.deg)
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

        assert str(e.value) == "Argument 'solarx' to function has 'myfunc_args' a 'unit' attribute without an 'is_equivalent' method. You may want to pass in an astropy Quantity instead."

def test_kwargs_dict():
    @u.quantity_input(solarx=u.arcsec, solary=u.deg)
    def myfunc_args(**kwargs):
        return kwargs['solarx'], kwargs['solary']

    solarx, solary = myfunc_args(solary=10*u.deg, solarx=1*u.arcsec)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.deg

def test_kwargs_input():
    @u.quantity_input(solarx=u.arcsec, solary=u.deg)
    def myfunc_args(solarx=1*u.arcsec, solary=1*u.deg):
        return solarx, solary

    kwargs = {'solarx':10*u.arcsec, 'solary':10*u.deg}
    solarx, solary = myfunc_args(**kwargs)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.deg

def test_kwargs_extra():
    @u.quantity_input(solary=u.deg)
    def myfunc_args(solarx, **kwargs):
        return solarx

    solarx, solary = myfunc_args(1*u.deg)

    assert isinstance(solarx, u.Quantity)

    assert solarx.unit == u.deg


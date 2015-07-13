# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from functools import wraps
from textwrap import dedent

from ... import units as u
from ...extern import six
from ...tests.helper import pytest


def py3only(func):
    if not six.PY3:
        return pytest.mark.skipif('not six.PY3')(func)
    else:
        @wraps(func)
        def wrapper(*args, **kwargs):
            if func.__doc__ is None:
                pytest.skip('unable to run this test due to missing '
                            'docstrings (maybe the module was compiled with '
                            'optimization flags?)')

            code = compile(dedent(func.__doc__), __file__, 'exec')
            # This uses an unqualified exec statement illegally in Python 2,
            # but perfectly allowed in Python 3 so in fact we eval the exec
            # call :)
            eval('exec(code)')

        return wrapper


@py3only
def test_args3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary: u.arcsec):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 1*u.arcsec)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.arcsec
    """


@py3only
def test_args_noconvert3():
    """
    @u.quantity_input()
    def myfunc_args(solarx: u.arcsec, solary: u.arcsec):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.deg, 1*u.arcmin)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.deg
    assert solary.unit == u.arcmin
    """


@py3only
def test_args_nonquantity3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 100)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)

    assert solarx.unit == u.arcsec
    """


@py3only
def test_arg_equivalencies3():
    """
    @u.quantity_input(equivalencies=u.mass_energy())
    def myfunc_args(solarx: u.arcsec, solary: u.eV):
        return solarx, solary+(10*u.J)  # Add an energy to check equiv is working

    solarx, solary = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.gram
    """


@py3only
def test_wrong_unit3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary: u.deg):
        return solarx, solary

    with pytest.raises(u.UnitsError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, 100*u.km)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' must be in units convertable to 'deg'."
    """


@py3only
def test_not_quantity3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary: u.deg):
        return solarx, solary

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, 100)
    assert str(e.value) == "Argument 'solary' to function has 'myfunc_args' no 'unit' attribute. You may want to pass in an astropy Quantity instead."
    """


@py3only
def test_decorator_override():
    """
    @u.quantity_input(solarx=u.arcsec)
    def myfunc_args(solarx: u.km, solary: u.arcsec):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec, 1*u.arcsec)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, u.Quantity)

    assert solarx.unit == u.arcsec
    assert solary.unit == u.arcsec
    """


@py3only
def test_kwargs3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary, myk: u.arcsec=1*u.arcsec):
        return solarx, solary, myk

    solarx, solary, myk = myfunc_args(1*u.arcsec, 100, myk=100*u.deg)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)
    assert isinstance(myk, u.Quantity)

    assert myk.unit == u.deg
    """


@py3only
def test_unused_kwargs3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary, myk: u.arcsec=1*u.arcsec, myk2=1000):
        return solarx, solary, myk, myk2

    solarx, solary, myk, myk2 = myfunc_args(1*u.arcsec, 100, myk=100*u.deg, myk2=10)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(solary, int)
    assert isinstance(myk, u.Quantity)
    assert isinstance(myk2, int)

    assert myk.unit == u.deg
    assert myk2 == 10
    """


@py3only
def test_kwarg_equivalencies3():
    """
    @u.quantity_input(equivalencies=u.mass_energy())
    def myfunc_args(solarx: u.arcsec, energy: u.eV=10*u.eV):
        return solarx, energy+(10*u.J)  # Add an energy to check equiv is working

    solarx, energy = myfunc_args(1*u.arcsec, 100*u.gram)

    assert isinstance(solarx, u.Quantity)
    assert isinstance(energy, u.Quantity)

    assert solarx.unit == u.arcsec
    assert energy.unit == u.gram
    """


@py3only
def test_kwarg_wrong_unit3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary: u.deg=10*u.deg):
        return solarx, solary

    with pytest.raises(u.UnitsError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100*u.km)
    assert str(e.value) == "Argument 'solary' to function 'myfunc_args' must be in units convertable to 'deg'."
    """


@py3only
def test_kwarg_not_quantity3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary: u.deg=10*u.deg):
        return solarx, solary

    with pytest.raises(TypeError) as e:
        solarx, solary = myfunc_args(1*u.arcsec, solary=100)
    assert str(e.value) == "Argument 'solary' to function has 'myfunc_args' no 'unit' attribute. You may want to pass in an astropy Quantity instead."
    """


@py3only
def test_kwarg_default3():
    """
    @u.quantity_input
    def myfunc_args(solarx: u.arcsec, solary: u.deg=10*u.deg):
        return solarx, solary

    solarx, solary = myfunc_args(1*u.arcsec)
    """

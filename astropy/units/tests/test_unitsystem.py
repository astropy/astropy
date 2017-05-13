# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
    Test the unit system functionality.
"""

from __future__ import absolute_import, unicode_literals, division, print_function


# Third party
from ... import units as u
from ...constants import G,c
import pytest

# This package
from ..unitsystem import UnitSystem, DimensionlessUnitSystem

def test_create():
    # simple system
    usys = UnitSystem(u.kpc, u.Myr, u.radian, u.Msun)

    # make sure length, time, mass, angle passed in
    with pytest.raises(ValueError):
        UnitSystem(u.kpc, u.Myr, u.radian)

    with pytest.raises(ValueError):
        UnitSystem(u.kpc, u.Myr, u.Msun)

    with pytest.raises(ValueError):
        UnitSystem(u.kpc, u.radian, u.Msun)

    with pytest.raises(ValueError):
        UnitSystem(u.Myr, u.radian, u.Msun)

    # as a single tuple instead of args
    usys2 = UnitSystem((u.kpc, u.Myr, u.radian, u.Msun))
    assert usys2 == usys

    # from an existing usys
    usys3 = UnitSystem(usys)
    assert usys3 == usys

def test_constants():
    usys = UnitSystem(u.kpc, u.Myr, u.radian, u.Msun)

    # try two constants
    assert usys.get_constant('G') == G.decompose([u.kpc, u.Myr, u.radian, u.Msun])
    assert usys.get_constant('c') == c.decompose([u.kpc, u.Myr, u.radian, u.Msun])

def test_dimensionless():
    usys = DimensionlessUnitSystem()
    assert usys['dimensionless'] == u.one
    assert usys['length'] == u.one

    with pytest.raises(ValueError):
        (15*u.kpc).to(usys)

    assert (1*u.one).to(usys) == 1*u.one

def test_compare_eq():
    usys1 = UnitSystem(u.kpc, u.Myr, u.radian, u.Msun, u.mas/u.yr)
    usys1_copy = UnitSystem(u.kpc, u.Myr, u.radian, u.Msun, u.mas/u.yr)

    usys2 = UnitSystem(u.kpc, u.Myr, u.radian, u.Msun, u.kiloarcsecond/u.yr)
    usys3 = UnitSystem(u.kpc, u.Myr, u.radian, u.kg, u.mas/u.yr)

    assert usys1 == usys1_copy
    assert usys1_copy == usys1

    assert usys1 != usys2
    assert usys2 != usys1

    assert usys1 != usys3
    assert usys3 != usys1

def test_convert():
    q = 15. * u.kpc/u.Myr

    usys1 = UnitSystem(u.kpc, u.Myr, u.radian, u.Msun)
    usys2 = UnitSystem(u.kpc, u.Myr, u.radian, u.Msun, u.km/u.s)

    assert q.to(usys1).unit == (u.kpc/u.Myr)
    assert q.to(usys2).unit == (u.km/u.s)

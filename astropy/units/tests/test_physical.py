# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Regression tests for the physical_type support in the units package
"""




from ... import units as u
from ...units import physical
from ...constants import hbar
from ...tests.helper import raises


def test_simple():
    assert u.m.physical_type == 'length'


def test_power():
    assert (u.cm ** 3).physical_type == 'volume'


def test_speed():
    assert (u.km / u.h).physical_type == 'speed'


def test_unknown():
    assert (u.m * u.s).physical_type == 'unknown'


def test_dimensionless():
    assert (u.m / u.m).physical_type == 'dimensionless'


def test_angular_momentum():
    assert hbar.unit.physical_type == 'angular momentum'


def test_flam():
    flam = u.erg / (u.cm**2 * u.s * u.AA)
    assert flam.physical_type == 'spectral flux density wav'


def test_photlam():
    photlam = u.photon / (u.cm ** 2 * u.s * u.AA)
    assert photlam.physical_type == 'photon flux density wav'


def test_photnu():
    photnu = u.photon / (u.cm ** 2 * u.s * u.Hz)
    assert photnu.physical_type == 'photon flux density'


@raises(ValueError)
def test_redundant_physical_type():
    physical.def_physical_type(u.m, 'utter craziness')


def test_data_quantity():
    assert u.byte.physical_type == 'data quantity'
    assert u.bit.physical_type == 'data quantity'

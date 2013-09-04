# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Regression tests for the physical_type support in the units package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)


from ... import units as u
from ...units import physical
from ...constants import hbar
from ...tests.helper import raises


def test_simple():
    assert u.mile.physical_type == 'length'


def test_power():
    assert (u.foot ** 3).physical_type == 'volume'


def test_speed():
    assert (u.mile / u.h).physical_type == 'speed'


def test_unknown():
    assert (u.m * u.s).physical_type == 'unknown'


def test_dimensionless():
    assert (u.m / u.m).physical_type == 'dimensionless'


def test_angular_momentum():
    assert hbar.unit.physical_type == 'angular momentum'


def test_flam():
    assert ((u.erg / (u.cm**2 * u.s * u.AA)).physical_type ==
            'spectral flux density wav')


@raises(ValueError)
def test_redundant_physical_type():
    physical.def_physical_type(u.m, 'utter craziness')

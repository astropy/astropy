# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Regression tests for the ptype support in the units package
"""
from __future__ import absolute_import, unicode_literals, division, print_function


from ... import units as u


def test_simple():
    assert u.mile.ptype == 'length'


def test_power():
    assert (u.foot ** 3).ptype == 'volume'


def test_speed():
    assert (u.mile / u.h).ptype == 'speed'


def test_unknown():
    assert (u.m * u.s).ptype == 'unknown'

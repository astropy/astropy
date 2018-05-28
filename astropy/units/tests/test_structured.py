# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test Structured units and quantities.
"""
from ... import units as u
from ...units import StructuredUnit


class TestStructuredUnitBasics:
    def test_initialization(self):
        p_unit = u.m
        v_unit = u.m / u.s
        t_unit = u.s
        su = StructuredUnit((p_unit, v_unit), [('p', 'O'), ('v', 'O')])
        assert su['p'] is p_unit
        assert su['v'] is v_unit
        su2 = StructuredUnit(((p_unit, v_unit), (t_unit)),
                             [('pv', [('p', 'O'), ('v', 'O')]), ('t', 'O')])
        assert isinstance(su2['pv'], StructuredUnit)
        assert su2['pv']['p'] is p_unit
        assert su2['pv']['v'] is v_unit
        assert su2['t'] is t_unit
        assert su2['pv'] == su

    def test_looks_like_unit(self):
        p_unit = u.m
        v_unit = u.m / u.s
        su = StructuredUnit((p_unit, v_unit), [('p', 'O'), ('v', 'O')])
        assert u.Unit(su) is su


class TestStructuredQuantity:
    def test_initialization(self):
        pass

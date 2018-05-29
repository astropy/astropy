# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test Structured units and quantities.
"""
import numpy as np

from ... import units as u
from ...units import StructuredUnit, StructuredQuantity, Unit, Quantity


class TestStructuredUnitBasics:
    def setup(self):
        self.pv_dtype = np.dtype([('p', 'f8'), ('v', 'f8')])
        self.pv_t_dtype = np.dtype([('pv', self.pv_dtype), ('t', 'f8')])

    def test_initialization_and_keying(self):
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
        su3 = StructuredUnit(('km', 'km/s'), [('p', 'O'), ('v', 'O')])
        assert isinstance(su3['p'], u.UnitBase)
        assert isinstance(su3['v'], u.UnitBase)

    def test_looks_like_unit(self):
        p_unit = u.m
        v_unit = u.m / u.s
        su = StructuredUnit((p_unit, v_unit), [('p', 'O'), ('v', 'O')])
        assert Unit(su) is su

    def test_initialize_with_float_dtype(self):
        su = StructuredUnit((('km', 'km/s'), 'yr'), self.pv_t_dtype)
        assert isinstance(su['pv'], StructuredUnit)
        assert isinstance(su['pv']['p'], u.UnitBase)
        assert isinstance(su['t'], u.UnitBase)
        assert su['pv']['v'] == u.km / u.s
        su = StructuredUnit((('km', 'km/s'), 'yr'), self.pv_t_dtype)
        assert isinstance(su['pv'], StructuredUnit)
        assert isinstance(su['pv']['p'], u.UnitBase)
        assert isinstance(su['t'], u.UnitBase)
        assert su['pv']['v'] == u.km / u.s


class TestStructuredQuantity:
    def setup(self):
        p_unit = u.km
        v_unit = u.km / u.s
        t_unit = u.s
        self.pv_unit = StructuredUnit((p_unit, v_unit),
                                      [('p', 'O'), ('v', 'O')])
        self.pv_t_unit = StructuredUnit(((p_unit, v_unit), (t_unit)),
                                        [('pv', [('p', 'O'), ('v', 'O')]),
                                         ('t', 'O')])
        self.pv_dtype = np.dtype([('p', 'f8'), ('v', 'f8')])
        self.pv_t_dtype = np.dtype([('pv', self.pv_dtype), ('t', 'f8')])
        self.pv = np.array([(1., 0.25), (2., 0.5), (3., 0.75)],
                           self.pv_dtype)
        self.pv_t = np.array([((4., 2.5), 0.),
                              ((5., 5.0), 1.),
                              ((6., 7.5), 2.)], self.pv_t_dtype)

    def test_initialization_and_keying(self):
        q_pv = StructuredQuantity(self.pv, self.pv_unit)
        q_p = q_pv['p']
        assert isinstance(q_p, Quantity)
        assert isinstance(q_p.unit, u.UnitBase)
        assert np.all(q_p == self.pv['p'] * self.pv_unit['p'])
        q_v = q_pv['v']
        assert isinstance(q_v, Quantity)
        assert isinstance(q_v.unit, u.UnitBase)
        assert np.all(q_v == self.pv['v'] * self.pv_unit['v'])
        q_pv_t = StructuredQuantity(self.pv_t, self.pv_t_unit)
        q_t = q_pv_t['t']
        assert np.all(q_t == self.pv_t['t'] * self.pv_t_unit['t'])
        q_pv2 = q_pv_t['pv']
        assert isinstance(q_pv2, StructuredQuantity)
        assert q_pv2.unit == self.pv_unit

    def test_initialization_with_unit_tuples(self):
        q_pv_t = StructuredQuantity(self.pv_t, (('km', 'km/s'), 's'))
        assert isinstance(q_pv_t.unit, StructuredUnit)
        assert q_pv_t.unit == self.pv_t_unit

    def test_getitem(self):
        q_pv_t = StructuredQuantity(self.pv_t, self.pv_t_unit)
        q_pv_t01 = q_pv_t[:2]
        assert isinstance(q_pv_t01, StructuredQuantity)
        assert q_pv_t01.unit == q_pv_t.unit
        assert np.all(q_pv_t01['t'] == q_pv_t['t'][:2])
        q_pv_t1 = q_pv_t[1]
        assert isinstance(q_pv_t1, StructuredQuantity)
        assert q_pv_t1.unit == q_pv_t.unit
        assert q_pv_t1.shape is ()
        assert q_pv_t1['t'] == q_pv_t['t'][1]

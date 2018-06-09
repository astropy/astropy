# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test Structured units and quantities.
"""
import numpy as np

from ... import units as u
from ...units import StructuredUnit, StructuredQuantity, Unit, Quantity


class StructuredTestBase:
    def setup(self):
        self.pv_dtype = np.dtype([('p', 'f8'), ('v', 'f8')])
        self.pv_t_dtype = np.dtype([('pv', self.pv_dtype), ('t', 'f8')])
        self.p_unit = u.km
        self.v_unit = u.km / u.s
        self.t_unit = u.s
        self.pv_dtype = np.dtype([('p', 'f8'), ('v', 'f8')])
        self.pv_t_dtype = np.dtype([('pv', self.pv_dtype), ('t', 'f8')])
        self.pv = np.array([(1., 0.25), (2., 0.5), (3., 0.75)],
                           self.pv_dtype)
        self.pv_t = np.array([((4., 2.5), 0.),
                              ((5., 5.0), 1.),
                              ((6., 7.5), 2.)], self.pv_t_dtype)


class StructuredTestBaseWithUnits(StructuredTestBase):
    def setup(self):
        super().setup()
        self.pv_unit = StructuredUnit((self.p_unit, self.v_unit),
                                      [('p', 'O'), ('v', 'O')])
        self.pv_t_unit = StructuredUnit(
            ((self.p_unit, self.v_unit), (self.t_unit)),
            [('pv', [('p', 'O'), ('v', 'O')]), ('t', 'O')])


class TestStructuredUnitBasics(StructuredTestBase):

    def test_initialization_and_keying(self):
        su = StructuredUnit((self.p_unit, self.v_unit),
                            [('p', 'O'), ('v', 'O')])
        assert su['p'] is self.p_unit
        assert su['v'] is self.v_unit
        su2 = StructuredUnit(((self.p_unit, self.v_unit), (self.t_unit)),
                             [('pv', [('p', 'O'), ('v', 'O')]), ('t', 'O')])
        assert isinstance(su2['pv'], StructuredUnit)
        assert su2['pv']['p'] is self.p_unit
        assert su2['pv']['v'] is self.v_unit
        assert su2['t'] is self.t_unit
        assert su2['pv'] == su
        su3 = StructuredUnit(('AU', 'AU/day'), [('p', 'O'), ('v', 'O')])
        assert isinstance(su3['p'], u.UnitBase)
        assert isinstance(su3['v'], u.UnitBase)

    def test_looks_like_unit(self):
        su = StructuredUnit((self.p_unit, self.v_unit),
                            [('p', 'O'), ('v', 'O')])
        assert Unit(su) is su

    def test_initialize_with_float_dtype(self):
        su = StructuredUnit(('AU', 'AU/d'), self.pv_dtype)
        assert isinstance(su['p'], u.UnitBase)
        assert isinstance(su['v'], u.UnitBase)
        assert su['p'] == u.AU
        assert su['v'] == u.AU / u.day
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

    def test_initialize_single_field(self):
        su = StructuredUnit('AU', [('p', 'O')])
        assert isinstance(su, StructuredUnit)
        assert isinstance(su['p'], u.UnitBase)
        assert su['p'] == u.AU
        su = StructuredUnit('AU')
        assert isinstance(su, StructuredUnit)
        assert isinstance(su['f0'], u.UnitBase)
        assert su['f0'] == u.AU


class TestStructuredUnitMethods(StructuredTestBaseWithUnits):
    def test_physical_type_id(self):
        pv_ptid = self.pv_unit._get_physical_type_id()
        expected = np.array((self.pv_unit['p']._get_physical_type_id(),
                             self.pv_unit['v']._get_physical_type_id()),
                            dtype=self.pv_unit.dtype)
        assert (pv_ptid == expected)
        pv_t_ptid = self.pv_t_unit._get_physical_type_id()
        expected = np.array(((self.pv_unit['p']._get_physical_type_id(),
                              self.pv_unit['v']._get_physical_type_id()),
                             self.pv_t_unit['t']._get_physical_type_id()),
                            dtype=self.pv_t_unit.dtype)
        assert (pv_t_ptid == expected)

    def test_physical_type(self):
        pv_pt = self.pv_unit.physical_type
        assert pv_pt == np.array(('length', 'speed'), self.pv_unit.dtype)
        assert pv_pt.item() == ('length', 'speed')
        pv_t_pt = self.pv_t_unit.physical_type
        assert pv_t_pt.item() == (('length', 'speed'), 'time')

    def test_si(self):
        pv_t_si = self.pv_t_unit.si
        assert pv_t_si == self.pv_t_unit
        assert pv_t_si['pv']['v'].scale == 1000

    def test_cgs(self):
        pv_t_cgs = self.pv_t_unit.cgs
        assert pv_t_cgs == self.pv_t_unit
        assert pv_t_cgs['pv']['v'].scale == 100000

    def test_decompose(self):
        pv_t_decompose = self.pv_t_unit.decompose()
        assert pv_t_decompose['pv']['v'].scale == 1000

    def test_is_equivalent(self):
        assert self.pv_unit.is_equivalent(('AU', 'AU/day'))
        assert not self.pv_unit.is_equivalent('m')
        assert not self.pv_unit.is_equivalent(('AU', 'AU'))

    def test_conversion(self):
        pv1 = self.pv_unit.to(('AU', 'AU/day'), self.pv)
        assert isinstance(pv1, np.ndarray)
        assert pv1.dtype == self.pv.dtype
        assert np.all(pv1['p'] * u.AU == self.pv['p'] * self.p_unit)
        assert np.all(pv1['v'] * u.AU / u.day == self.pv['v'] * self.v_unit)
        # Names should be from value.
        su2 = StructuredUnit((self.p_unit, self.v_unit),
                             [('position', 'O'), ('velocity', 'O')])
        pv2 = su2.to(('Mm', 'mm/s'), self.pv)
        assert pv2.dtype.names == ('p', 'v')
        assert pv2.dtype == self.pv.dtype
        # Check recursion.
        pv_t1 = self.pv_t_unit.to((('AU', 'AU/day'), 'Myr'), self.pv_t)
        assert isinstance(pv_t1, np.ndarray)
        assert pv_t1.dtype == self.pv_t.dtype
        assert np.all(pv_t1['pv']['p'] * u.AU ==
                      self.pv_t['pv']['p'] * self.p_unit)
        assert np.all(pv_t1['pv']['v'] * u.AU / u.day ==
                      self.pv_t['pv']['v'] * self.v_unit)
        assert np.all(pv_t1['t'] * u.Myr == self.pv_t['t'] * self.t_unit)


class TestStructuredQuantity(StructuredTestBaseWithUnits):
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

    def test_value(self):
        q_pv_t = StructuredQuantity(self.pv_t, self.pv_t_unit)
        value = q_pv_t.value
        assert type(value) is np.ndarray
        assert np.all(value == self.pv_t)
        value1 = q_pv_t[1].value
        assert type(value1) is np.void
        assert np.all(value1 == self.pv_t[1])

    def test_conversion(self):
        q_pv = StructuredQuantity(self.pv, self.pv_unit)
        q1 = q_pv.to(('AU', 'AU/day'))
        assert isinstance(q1, StructuredQuantity)
        assert q1['p'].unit == u.AU
        assert q1['v'].unit == u.AU / u.day
        assert np.all(q1['p'] == q_pv['p'].to(u.AU))
        assert np.all(q1['v'] == q_pv['v'].to(u.AU/u.day))
        pv1 = q_pv.to_value(('AU', 'AU/day'))
        assert type(pv1) is np.ndarray
        assert np.all(pv1['p'] == q_pv['p'].to_value(u.AU))
        assert np.all(pv1['v'] == q_pv['v'].to_value(u.AU/u.day))
        pv11 = q_pv[1].to_value(('AU', 'AU/day'))
        assert type(pv11) is np.void
        assert pv11 == pv1[1]
        q_pv_t = StructuredQuantity(self.pv_t, self.pv_t_unit)
        q2 = q_pv_t.to((('kpc', 'kpc/Myr'), 'Myr'))
        assert q2['pv']['p'].unit == u.kpc
        assert q2['pv']['v'].unit == u.kpc / u.Myr
        assert q2['t'].unit == u.Myr
        assert np.all(q2['pv']['p'] == q_pv_t['pv']['p'].to(u.kpc))
        assert np.all(q2['pv']['v'] == q_pv_t['pv']['v'].to(u.kpc/u.Myr))
        assert np.all(q2['t'] == q_pv_t['t'].to(u.Myr))

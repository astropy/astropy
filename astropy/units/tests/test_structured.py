# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test Structured units and quantities.
"""
import pytest
import numpy as np

from ... import units as u
from ...units import (StructuredUnit, StructuredQuantity,
                      Unit, UnitBase, Quantity)


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
                                      ('p', 'v'))
        self.pv_t_unit = StructuredUnit((self.pv_unit, self.t_unit),
                                        ('pv', 't'))


class TestStructuredUnitBasics(StructuredTestBase):

    def test_initialization_and_keying(self):
        su = StructuredUnit((self.p_unit, self.v_unit), ('p', 'v'))
        assert su['p'] is self.p_unit
        assert su['v'] is self.v_unit
        su2 = StructuredUnit((su, self.t_unit), ('pv', 't'))
        assert isinstance(su2['pv'], StructuredUnit)
        assert su2['pv']['p'] is self.p_unit
        assert su2['pv']['v'] is self.v_unit
        assert su2['t'] is self.t_unit
        assert su2['pv'] == su
        su3 = StructuredUnit(('AU', 'AU/day'), ('p', 'v'))
        assert isinstance(su3['p'], UnitBase)
        assert isinstance(su3['v'], UnitBase)
        su4 = StructuredUnit('AU, AU/day', ('p', 'v'))
        assert su4['p'] == u.AU
        assert su4['v'] == u.AU / u.day
        su5 = StructuredUnit(('AU', 'AU/day'))
        assert su5.dtype.names == ('f0', 'f1')
        assert su5['f0'] == u.AU
        assert su5['f1'] == u.AU / u.day

    def test_recursive_initialization(self):
        su = StructuredUnit(((self.p_unit, self.v_unit), self.t_unit),
                            (('p', 'v'), 't'))
        assert isinstance(su['pv'], StructuredUnit)
        assert su['pv']['p'] is self.p_unit
        assert su['pv']['v'] is self.v_unit
        assert su['t'] is self.t_unit
        su2 = StructuredUnit(((self.p_unit, self.v_unit), self.t_unit),
                             (['p_v', ('p', 'v')], 't'))
        assert isinstance(su2['p_v'], StructuredUnit)
        assert su2['p_v']['p'] is self.p_unit
        assert su2['p_v']['v'] is self.v_unit
        assert su2['t'] is self.t_unit
        su3 = StructuredUnit((('AU', 'AU/day'), 'yr'),
                             (['p_v', ('p', 'v')], 't'))
        assert isinstance(su3['p_v'], StructuredUnit)
        assert su3['p_v']['p'] == u.AU
        assert su3['p_v']['v'] == u.AU / u.day
        assert su3['t'] == u.yr
        su4 = StructuredUnit('(AU, AU/day), yr', (('p', 'v'), 't'))
        assert isinstance(su4['pv'], StructuredUnit)
        assert su4['pv']['p'] == u.AU
        assert su4['pv']['v'] == u.AU / u.day
        assert su4['t'] == u.yr

    def test_looks_like_unit(self):
        su = StructuredUnit((self.p_unit, self.v_unit), ('p', 'v'))
        assert Unit(su) is su

    def test_initialize_with_float_dtype(self):
        su = StructuredUnit(('AU', 'AU/d'), self.pv_dtype)
        assert isinstance(su['p'], UnitBase)
        assert isinstance(su['v'], UnitBase)
        assert su['p'] == u.AU
        assert su['v'] == u.AU / u.day
        su = StructuredUnit((('km', 'km/s'), 'yr'), self.pv_t_dtype)
        assert isinstance(su['pv'], StructuredUnit)
        assert isinstance(su['pv']['p'], UnitBase)
        assert isinstance(su['t'], UnitBase)
        assert su['pv']['v'] == u.km / u.s
        su = StructuredUnit('(km, km/s), yr', self.pv_t_dtype)
        assert isinstance(su['pv'], StructuredUnit)
        assert isinstance(su['pv']['p'], UnitBase)
        assert isinstance(su['t'], UnitBase)
        assert su['pv']['v'] == u.km / u.s

    def test_initialize_single_field(self):
        su = StructuredUnit('AU', 'p')
        assert isinstance(su, StructuredUnit)
        assert isinstance(su['p'], UnitBase)
        assert su['p'] == u.AU
        su = StructuredUnit('AU')
        assert isinstance(su, StructuredUnit)
        assert isinstance(su['f0'], UnitBase)
        assert su['f0'] == u.AU

    def test_parsing(self):
        su = Unit('AU, AU/d')
        assert isinstance(su, StructuredUnit)
        assert isinstance(su['f0'], UnitBase)
        assert isinstance(su['f1'], UnitBase)
        assert su['f0'] == u.AU
        assert su['f1'] == u.AU/u.day
        su2 = Unit('AU, AU/d, yr')
        assert isinstance(su2, StructuredUnit)
        assert su2 == StructuredUnit(('AU', 'AU/d', 'yr'))
        su2a = Unit('(AU, AU/d, yr)')
        assert isinstance(su2a, StructuredUnit)
        assert su2a == su2
        su3 = Unit('(km, km/s), yr')
        assert isinstance(su3, StructuredUnit)
        assert su3 == StructuredUnit((('km', 'km/s'), 'yr'))
        su4 = Unit('km,')
        assert isinstance(su4, StructuredUnit)
        assert su4 == StructuredUnit((u.km,))
        su5 = Unit('(m,s),')
        assert isinstance(su5, StructuredUnit)
        assert su5 == StructuredUnit(((u.m, u.s),))
        ldbody_unit = Unit('Msun, 0.5rad^2, (au, au/day)')
        assert ldbody_unit == StructuredUnit(
            (u.Msun, Unit(u.rad**2 / 2), (u.AU, u.AU / u.day)))

    def test_str(self):
        su = StructuredUnit(((u.km, u.km/u.s), u.yr))
        assert str(su) == '((km, km / s), yr)'
        assert Unit(str(su)) == su

    def test_repr(self):
        su = StructuredUnit(((u.km, u.km/u.s), u.yr))
        assert repr(su) == 'Unit("((km, km / s), yr)")'
        assert eval(repr(su)) == su


class TestStructuredUnitMethods(StructuredTestBaseWithUnits):
    def test_physical_type_id(self):
        pv_ptid = self.pv_unit._get_physical_type_id()
        expected = np.array((self.pv_unit['p']._get_physical_type_id(),
                             self.pv_unit['v']._get_physical_type_id()),
                            self.pv_unit.dtype)[()]
        assert (pv_ptid == expected)
        pv_t_ptid = self.pv_t_unit._get_physical_type_id()
        expected2 = np.array((self.pv_unit._get_physical_type_id(),
                              self.t_unit._get_physical_type_id()),
                             self.pv_t_unit.dtype)[()]
        assert (pv_t_ptid == expected2)

    def test_physical_type(self):
        pv_pt = self.pv_unit.physical_type
        assert pv_pt == np.array(('length', 'speed'), self.pv_unit.dtype)[()]
        pv_t_pt = self.pv_t_unit.physical_type
        assert pv_t_pt == np.array((pv_pt, 'time'), self.pv_t_unit.dtype)[()]

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
                             ('position', 'velocity'))
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
        # Passing in tuples should work.
        pv_t2 = self.pv_t_unit.to((('AU', 'AU/day'), 'Myr'),
                                  ((1., 0.1), 10.))
        assert pv_t2['pv']['p'] == self.p_unit.to('AU', 1.)
        assert pv_t2['pv']['v'] == self.v_unit.to('AU/day', 0.1)
        assert pv_t2['t'] == self.t_unit.to('Myr', 10.)
        pv_t3 = self.pv_t_unit.to((('AU', 'AU/day'), 'Myr'),
                                  [((1., 0.1), 10.),
                                   ((2., 0.2), 20.)])
        assert np.all(pv_t3['pv']['p'] == self.p_unit.to('AU', [1., 2.]))
        assert np.all(pv_t3['pv']['v'] == self.v_unit.to('AU/day', [0.1, 0.2]))
        assert np.all(pv_t3['t'] == self.t_unit.to('Myr', [10., 20.]))


class TestStructuredUnitArithmatic(StructuredTestBaseWithUnits):
    def test_multiplication(self):
        pv_times_au = self.pv_unit * u.au
        assert isinstance(pv_times_au, StructuredUnit)
        assert pv_times_au.dtype.names == ('p', 'v')
        assert pv_times_au['p'] == self.p_unit * u.AU
        assert pv_times_au['v'] == self.v_unit * u.AU
        au_times_pv = u.au * self.pv_unit
        assert au_times_pv == pv_times_au
        pv_times_au2 = self.pv_unit * 'au'
        assert pv_times_au2 == pv_times_au
        au_times_pv2 = 'AU' * self.pv_unit
        assert au_times_pv2 == pv_times_au
        with pytest.raises(TypeError):
            self.pv_unit * self.pv_unit
        with pytest.raises(TypeError):
            's,s' * self.pv_unit

    def test_division(self):
        pv_by_s = self.pv_unit / u.s
        assert isinstance(pv_by_s, StructuredUnit)
        assert pv_by_s.dtype.names == ('p', 'v')
        assert pv_by_s['p'] == self.p_unit / u.s
        assert pv_by_s['v'] == self.v_unit / u.s
        pv_by_s2 = self.pv_unit / 's'
        assert pv_by_s2 == pv_by_s
        with pytest.raises(TypeError):
            1. / self.pv_unit
        with pytest.raises(TypeError):
            u.s / self.pv_unit


class TestStructuredQuantity(StructuredTestBaseWithUnits):
    def test_initialization_and_keying(self):
        q_pv = StructuredQuantity(self.pv, self.pv_unit)
        q_p = q_pv['p']
        assert isinstance(q_p, Quantity)
        assert isinstance(q_p.unit, UnitBase)
        assert np.all(q_p == self.pv['p'] * self.pv_unit['p'])
        q_v = q_pv['v']
        assert isinstance(q_v, Quantity)
        assert isinstance(q_v.unit, UnitBase)
        assert np.all(q_v == self.pv['v'] * self.pv_unit['v'])
        q_pv_t = StructuredQuantity(self.pv_t, self.pv_t_unit)
        q_t = q_pv_t['t']
        assert np.all(q_t == self.pv_t['t'] * self.pv_t_unit['t'])
        q_pv2 = q_pv_t['pv']
        assert isinstance(q_pv2, StructuredQuantity)
        assert q_pv2.unit == self.pv_unit
        with pytest.raises(ValueError):
            StructuredQuantity(self.pv, self.pv_t_unit)
        with pytest.raises(ValueError):
            StructuredQuantity(self.pv_t, self.pv_unit)

    def test_initialization_with_unit_tuples(self):
        q_pv_t = StructuredQuantity(self.pv_t, (('km', 'km/s'), 's'))
        assert isinstance(q_pv_t.unit, StructuredUnit)
        assert q_pv_t.unit == self.pv_t_unit

    def test_initialization_with_string(self):
        q_pv_t = StructuredQuantity(self.pv_t, '(km, km/s), s')
        assert isinstance(q_pv_t.unit, StructuredUnit)
        assert q_pv_t.unit == self.pv_t_unit

    def test_initialization_by_multiplication_with_unit(self):
        q_pv_t = self.pv_t * self.pv_t_unit
        assert q_pv_t.unit is self.pv_t_unit
        assert np.all(q_pv_t.value == self.pv_t)
        q_pv_t2 = self.pv_t_unit * self.pv_t
        assert np.all(q_pv_t2 == q_pv_t)

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

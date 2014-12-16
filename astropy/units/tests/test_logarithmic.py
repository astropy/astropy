# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
    Test the Logarithmic Units and Quantities
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from ...extern import six

import itertools
import numpy as np
from numpy.testing.utils import assert_allclose

from ...tests.helper import pytest
from ... import units as u

lu_subclasses = [u.DexUnit, u.MagUnit, u.DecibelUnit]

lq_subclasses = [u.Dex, u.Magnitude, u.Decibel]

pu_sample = (u.dimensionless_unscaled, u.m, u.g/u.s**2, u.Jy)


class TestLogUnitCreation(object):

    def test_logarithmic_units(self):
        """Check logarithmic units are all set up correctly"""
        assert u.dB.to(u.dex) == 0.1
        assert u.dex.to(u.mag) == -2.5
        assert u.mag.to(u.dB) == -4

    @pytest.mark.parametrize('lu_cls, physical_unit',
                             itertools.product(lu_subclasses, pu_sample))
    def test_subclass_creation(self, lu_cls, physical_unit):
        """Create a LogUnit subclass object for given physcal unit,
        and do basic check that output is right"""
        lu1 = lu_cls(physical_unit)
        assert lu1.physical_unit == physical_unit
        assert lu1.function_unit == lu1._default_function_unit

        lu2 = lu_cls(physical_unit,
                     function_unit=2*lu1._default_function_unit)
        assert lu2.physical_unit == physical_unit
        assert lu2.function_unit == u.Unit(2*lu2._default_function_unit)

        with pytest.raises(ValueError):
            lu_cls(physical_unit, u.m)


class TestLogUnitConversion(object):
    @pytest.mark.parametrize('lu_cls, physical_unit',
                             itertools.product(lu_subclasses, pu_sample))
    def test_physical_unit_conversion(self, lu_cls, physical_unit):
        """Check various LogUnit subclasses are equivalent and convertible
        to their non-log counterparts"""
        lu1 = lu_cls(physical_unit)
        assert lu1.is_equivalent(physical_unit)
        assert lu1.to(physical_unit, 0.) == 1.

        assert physical_unit.is_equivalent(lu1)
        assert physical_unit.to(lu1, 1.) == 0.

        pu = u.Unit(8.*physical_unit)
        assert lu1.is_equivalent(physical_unit)
        assert lu1.to(pu, 0.) == 0.125

        assert pu.is_equivalent(lu1)
        assert_allclose(pu.to(lu1, 0.125), 0., atol=1.e-15)

        # round-trip
        value = np.linspace(0., 10., 6)
        assert_allclose(pu.to(lu1, lu1.to(pu, value)), value, atol=1.e-15)
        # and we're not just returning True all the time
        pu2 = u.g
        assert not lu1.is_equivalent(pu2)
        with pytest.raises(u.UnitsError):
            lu1.to(pu2)

        assert not pu2.is_equivalent(lu1)
        with pytest.raises(u.UnitsError):
            pu2.to(lu1)

    @pytest.mark.parametrize('lu_cls', lu_subclasses)
    def test_container_unit_conversion(self, lu_cls):
        """Check that conversion to logarithmic units (u.mag, u.dB, u.dex)
        is only possible when the physical unit is dimensionless"""
        values = np.linspace(0.,10.,6.)
        lu1 = lu_cls(u.dimensionless_unscaled)
        assert lu1.is_equivalent(lu1.function_unit)
        assert_allclose(lu1.to(lu1.function_unit, values), values)

        lu2 = lu_cls(u.Jy)
        assert not lu2.is_equivalent(lu2.function_unit)
        with pytest.raises(u.UnitsError):
            lu2.to(lu2.function_unit, values)

    @pytest.mark.parametrize(
        'flu_cls, tlu_cls, physical_unit',
        itertools.product(lu_subclasses, lu_subclasses, pu_sample))
    def test_subclass_conversion(self, flu_cls, tlu_cls, physical_unit):
        """Check various LogUnit subclasses are equivalent and convertible
        to each other if they correspond to equivalent physical units."""
        values = np.linspace(0.,10.,6.)
        flu = flu_cls(physical_unit)

        tlu = tlu_cls(physical_unit)
        assert flu.is_equivalent(tlu)
        assert_allclose(flu.to(tlu), flu.function_unit.to(tlu.function_unit))
        assert_allclose(flu.to(tlu, values),
                        values * flu.function_unit.to(tlu.function_unit))

        tlu2 = tlu_cls(u.Unit(100.*physical_unit))
        assert flu.is_equivalent(tlu2)
        # check roundtrip
        assert_allclose(flu.to(tlu2, tlu2.to(flu, values)), values, atol=1.e-15)

        tlu3 = tlu_cls(physical_unit.to_system(u.si)[0])
        assert flu.is_equivalent(tlu3)
        assert_allclose(flu.to(tlu3, tlu3.to(flu, values)), values, atol=1.e-15)

        tlu4 = tlu_cls(u.g)
        assert not flu.is_equivalent(tlu4)
        with pytest.raises(u.UnitsError):
            flu.to(tlu4, values)


class TestLogUnitArithmetic(object):
    def test_multiplication_division(self):
        """Check that multiplication/division with other units is only
        possible when the physical unit is dimensionless, and that this
        turns the unit into a normal one"""
        lu1 = u.MagUnit(u.Jy)

        with pytest.raises(u.UnitsError):
            lu1 * u.m

        with pytest.raises(u.UnitsError):
            u.m * lu1

        with pytest.raises(u.UnitsError):
            lu1 / lu1

        for unit in (u.dimensionless_unscaled, u.m, u.mag, u.dex):
            with pytest.raises(u.UnitsError):
                lu1 / unit

        lu2 = u.MagUnit(u.dimensionless_unscaled)

        with pytest.raises(u.UnitsError):
            lu2 * lu1

        with pytest.raises(u.UnitsError):
            lu2 / lu1

        # but dimensionless_unscaled can be cancelled
        assert lu2 / lu2 == u.dimensionless_unscaled

        # with dimensionless, normal units are OK, but we return a plain unit
        tf = lu2 * u.m
        tr = u.m * lu2
        for t in (tf, tr):
            assert not isinstance(t, lu2.__class__)
            assert t == lu2.function_unit * u.m
            with u.set_enabled_equivalencies(u.logarithmic()):
                with pytest.raises(u.UnitsError):
                    t.to(lu2.physical_unit)
        # now we essentially have a LogUnit with a prefactor of 100,
        # so should be equivalent again
        t = tf / u.cm
        with u.set_enabled_equivalencies(u.logarithmic()):
            assert t.is_equivalent(lu2.function_unit)
            assert_allclose(t.to(u.dimensionless_unscaled, np.arange(3.)/100.),
                            lu2.to(lu2.physical_unit, np.arange(3.)))

        # also if we divide such that the lu1 is removed, we should
        # become a normal unit
        t2 = tf / lu2
        assert not isinstance(t2, lu2.__class__)
        assert t2 == u.m
        t3 = tf / lu2.function_unit
        assert not isinstance(t3, lu2.__class__)
        assert t3 == u.m

    @pytest.mark.parametrize('power', (2, 0.5, 1, 0))
    def test_raise_to_power(self, power):
        """Check that raising LogUnits to some power is only possible when the
        physical unit is dimensionless, and that conversion is turned off when
        the resulting logarithmic unit (such as mag**2) is incompatible."""
        lu1 = u.MagUnit(u.Jy)

        if power == 0:
            assert lu1 ** power == u.dimensionless_unscaled
        elif power == 1:
            assert lu1 ** power == lu1
        else:
            with pytest.raises(u.UnitsError):
                lu1 ** power

        # with dimensionless, though, it can work, but returns a normal unit
        lu2 = u.MagUnit(u.dimensionless_unscaled)

        t = lu2**power
        if power == 0:
            assert t == u.dimensionless_unscaled
        elif power == 1:
            assert t == lu2
        else:
            assert not isinstance(t, lu2.__class__)
            assert t == lu2.function_unit**power
            # also check we roundtrip
            t2 = t**(1./power)
            assert t2 == lu2.function_unit
            with u.set_enabled_equivalencies(u.logarithmic()):
                assert_allclose(t2.to(u.dimensionless_unscaled, np.arange(3.)),
                                lu2.to(lu2.physical_unit, np.arange(3.)))

    @pytest.mark.parametrize('other', pu_sample)
    def test_addition_subtraction_to_normal_units_fails(self, other):
        lu1 = u.MagUnit(u.Jy)
        with pytest.raises(u.UnitsError):
            lu1 + other

        with pytest.raises(u.UnitsError):
            lu1 - other

        with pytest.raises(u.UnitsError):
            other - lu1

    @pytest.mark.parametrize(
        'other', (u.mag, u.MagUnit(), u.MagUnit(u.Jy), u.MagUnit(u.m),
                  u.Unit(2*u.mag), u.MagUnit('', 2.*u.mag)))
    def test_addition_subtraction(self, other):
        """Check physical units are changed appropriately"""
        lu1 = u.MagUnit(u.Jy)
        other_pu = getattr(other, 'physical_unit', u.dimensionless_unscaled)

        lu_sf = lu1 + other
        assert lu_sf.is_equivalent(lu1.physical_unit * other_pu)

        lu_sr = other + lu1
        assert lu_sr.is_equivalent(lu1.physical_unit * other_pu)

        lu_df = lu1 - other
        assert lu_df.is_equivalent(lu1.physical_unit / other_pu)

        lu_dr = other - lu1
        assert lu_dr.is_equivalent(other_pu / lu1.physical_unit)

    def test_complicated_addition_subtraction(self):
        """for fun, a more complicated example of addition and subtraction"""
        dm0 = u.Unit('DM', 1./(4.*np.pi*(10.*u.pc)**2))
        lu_dm = u.MagUnit(dm0, u.mag)
        lu_absST = u.STmag - lu_dm
        assert lu_absST.is_equivalent(u.erg/u.s/u.AA)


class TestLogQuantityCreation(object):

    @pytest.mark.parametrize('lq, lu', zip(lq_subclasses, lu_subclasses))
    def test_logarithmic_quantities(self, lq, lu):
        """Check logarithmic quantities are all set up correctly"""
        assert lq._FunctionUnit == lu
        assert type(lu().function_quantity(1.)) is lq

    @pytest.mark.parametrize('lq_cls, physical_unit',
                             itertools.product(lq_subclasses, pu_sample))
    def test_subclass_creation(self, lq_cls, physical_unit):
        """Create LogQuantity subclass objects for some physcal units,
        and basic check on transformations"""
        value = np.arange(1.,10.)
        log_q = lq_cls(value * physical_unit)
        assert log_q.unit.physical_unit == physical_unit
        assert log_q.unit.function_unit == log_q.unit._default_function_unit
        assert_allclose(log_q.quantity.value, value)
        with pytest.raises(ValueError):
            lq_cls(value, physical_unit)

    @pytest.mark.parametrize(
        'unit', (u.mag, u.MagUnit(), u.MagUnit(u.Jy), u.MagUnit(u.m),
                 u.Unit(2*u.mag), u.MagUnit('', 2.*u.mag),
                 u.MagUnit(u.Jy, -1*u.mag), u.MagUnit(u.m, -2.*u.mag)))
    def test_different_units(self, unit):
        q = u.Magnitude(1.23, unit)
        assert q.unit.function_unit == getattr(unit, 'function_unit', unit)
        assert q.unit.physical_unit is getattr(unit, 'physical_unit',
                                               u.dimensionless_unscaled)

    @pytest.mark.parametrize(
        'unit', (u.MagUnit(), u.MagUnit(u.Jy), u.MagUnit(u.m),
                 u.MagUnit('', 2.*u.mag),
                 u.MagUnit(u.Jy, -1*u.mag), u.MagUnit(u.m, -2.*u.mag)))
    def test_indirect_creation(self, unit):
        q1 = 2.5 * unit
        assert isinstance(q1, u.Magnitude)
        assert q1.value == 2.5
        assert q1.unit == unit
        pv = 100. * unit.physical_unit
        q2 = pv * unit
        assert q2.unit == unit
        assert q2.unit.physical_unit == pv.unit
        assert q2.to(unit.physical_unit).value == 100.
        assert (q2.function_value / u.mag).to(1).value == -5.


class TestLogQuantityViews(object):
    def test_value_view(self):
        lq = u.Magnitude(np.arange(10.) * u.Jy)
        lq_value = lq.value
        assert type(lq_value) is np.ndarray
        lq_value[2] = -1.
        assert np.all(lq.value == lq_value)

        lq_fv = lq.function_value
        assert type(lq_fv) is u.Quantity
        assert lq_fv.unit is lq.unit.function_unit
        lq_fv[3] = -2. * lq_fv.unit
        assert np.all(lq.value == lq_fv.value)
        assert np.all(lq.value == lq_value)


class TestLogQuantitySlicing(object):
    def test_item_get_and_set(self):
        lq1 = u.Magnitude(np.arange(1., 11.)*u.Jy)
        assert lq1[9] == u.Magnitude(10.*u.Jy)
        lq1[2] = 100.*u.Jy
        assert lq1[2] == u.Magnitude(100.*u.Jy)
        with pytest.raises(u.UnitsError):
            lq1[2] = 100.*u.m
        with pytest.raises(u.UnitsError):
            lq1[2] = 100.*u.mag
        with pytest.raises(u.UnitsError):
            lq1[2] = u.Magnitude(100.*u.m)
        assert lq1[2] == u.Magnitude(100.*u.Jy)

    def test_slice_get_and_set(self):
        lq1 = u.Magnitude(np.arange(1., 10.)*u.Jy)
        lq1[2:4] = 100.*u.Jy
        assert np.all(lq1[2:4] == u.Magnitude(100.*u.Jy))
        with pytest.raises(u.UnitsError):
            lq1[2:4] = 100.*u.m
        with pytest.raises(u.UnitsError):
            lq1[2:4] = 100.*u.mag
        with pytest.raises(u.UnitsError):
            lq1[2:4] = u.Magnitude(100.*u.m)
        assert np.all(lq1[2] == u.Magnitude(100.*u.Jy))


class TestLogQuantityArithmetic(object):
    def test_multiplication_division(self):
        """Check that multiplication/division with other quanties is only
        possible when the physical unit is dimensionless, and that this turns
        the result into a normal quantity."""
        lq = u.Magnitude(np.arange(1., 11.)*u.Jy)

        with pytest.raises(u.UnitsError):
            lq * (1.*u.m)

        with pytest.raises(u.UnitsError):
            (1.*u.m) * lq

        with pytest.raises(u.UnitsError):
            lq / lq

        for unit in (u.m, u.mag, u.dex):
            with pytest.raises(u.UnitsError):
                lq / unit

        lq2 = u.Magnitude(np.arange(1, 11.))

        with pytest.raises(u.UnitsError):
            lq2 * lq

        with pytest.raises(u.UnitsError):
            lq2 / lq

        with pytest.raises(u.UnitsError):
            lq / lq2

        # but dimensionless_unscaled can be cancelled
        r = lq2 / u.Magnitude(2.)
        assert r.unit == u.dimensionless_unscaled
        assert np.all(r.value == lq2.value/2.)

        # with dimensionless, normal units OK, but return normal quantities
        tf = lq2 * u.m
        tr = u.m * lq2
        for t in (tf, tr):
            assert not isinstance(t, lq2.__class__)
            assert t.unit == lq2.unit.function_unit * u.m
            with u.set_enabled_equivalencies(u.logarithmic()):
                with pytest.raises(u.UnitsError):
                    t.to(lq2.unit.physical_unit)

        t = tf / (50.*u.cm)
        # now we essentially have the same quantity but with a prefactor of 2
        assert t.unit.is_equivalent(lq2.unit.function_unit)
        assert_allclose(t.to(lq2.unit.function_unit), lq2.function_value*2)

    @pytest.mark.parametrize('power', (2, 0.5, 1, 0))
    def test_raise_to_power(self, power):
        """Check that raising LogQuantities to some power is only possible when
        the physical unit is dimensionless, and that conversion is turned off
        when the resulting logarithmic unit (say, mag**2) is incompatible."""
        lq = u.Magnitude(np.arange(1., 4.)*u.Jy)

        if power == 0:
            assert np.all(lq ** power == 1.)
        elif power == 1:
            assert np.all(lq ** power == lq)
        else:
            with pytest.raises(u.UnitsError):
                lq ** power

        # with dimensionless, it works, but falls back to normal quantity
        # (except for power=1)
        lq2 = u.Magnitude(np.arange(10.))

        t = lq2**power
        if power == 0:
            assert t.unit is u.dimensionless_unscaled
            assert np.all(t.value == 1.)
        elif power == 1:
            assert np.all(t == lq2)
        else:
            assert not isinstance(t, lq2.__class__)
            assert t.unit == lq2.unit.function_unit ** power
            with u.set_enabled_equivalencies(u.logarithmic()):
                with pytest.raises(u.UnitsError):
                    t.to(u.dimensionless_unscaled)

    @pytest.mark.parametrize('other', pu_sample)
    def test_addition_subtraction_to_normal_units_fails(self, other):
        lq = u.Magnitude(np.arange(1.,10.)*u.Jy)
        q = 1.23 * other
        with pytest.raises(u.UnitsError):
            lq + q

        with pytest.raises(u.UnitsError):
            lq - q

        with pytest.raises(u.UnitsError):
            q - lq

    @pytest.mark.parametrize(
        'other', (1.23 * u.mag, 2.34 * u.MagUnit(),
                  u.Magnitude(3.45 * u.Jy), u.Magnitude(4.56 * u.m),
                  5.67 * u.Unit(2*u.mag), u.Magnitude(6.78, 2.*u.mag)))
    def test_addition_subtraction(self, other):
        """Check that addition/subtraction with quantities with magnitude or
        MagUnit units works, and that it changes the physical units
        appropriately"""
        lq = u.Magnitude(np.arange(1.,10.)*u.Jy)
        other_quantity = other.to(getattr(other.unit, 'physical_unit',
                                          u.dimensionless_unscaled),
                                  equivalencies=u.logarithmic())

        lq_sf = lq + other
        assert_allclose(lq_sf.quantity, lq.quantity * other_quantity)

        lq_sr = other + lq
        assert_allclose(lq_sr.quantity, lq.quantity * other_quantity)

        lq_df = lq - other
        assert_allclose(lq_df.quantity, lq.quantity / other_quantity)

        lq_dr = other - lq
        assert_allclose(lq_dr.quantity, other_quantity / lq.quantity)

    @pytest.mark.parametrize('other', pu_sample)
    def test_inplace_addition_subtraction_unit_checks(self, other):
        lu1 = u.MagUnit(u.Jy)
        lq1 = u.Magnitude(np.arange(1.,10.), lu1)
        with pytest.raises(u.UnitsError):
            lq1 += other

        assert np.all(lq1.value == np.arange(1., 10.))
        assert lq1.unit == lu1

        with pytest.raises(u.UnitsError):
            lq1 -= other

        assert np.all(lq1.value == np.arange(1., 10.))
        assert lq1.unit == lu1

    @pytest.mark.parametrize(
        'other', (1.23 * u.mag, 2.34 * u.MagUnit(),
                  u.Magnitude(3.45 * u.Jy), u.Magnitude(4.56 * u.m),
                  5.67 * u.Unit(2*u.mag), u.Magnitude(6.78, 2.*u.mag)))
    def test_inplace_addition_subtraction(self, other):
        """Check that inplace addition/subtraction with quantities with
        magnitude or MagUnit units works, and that it changes the physical
        units appropriately"""
        lq = u.Magnitude(np.arange(1.,10.)*u.Jy)
        other_quantity = other.to(getattr(other.unit, 'physical_unit',
                                          u.dimensionless_unscaled),
                                  equivalencies=u.logarithmic())
        lq_sf = lq.copy()
        lq_sf += other
        assert_allclose(lq_sf.quantity, lq.quantity * other_quantity)

        lq_df = lq.copy()
        lq_df -= other
        assert_allclose(lq_df.quantity, lq.quantity / other_quantity)

    def test_complicated_addition_subtraction(self):
        """for fun, a more complicated example of addition and subtraction"""
        dm0 = u.Unit('DM', 1./(4.*np.pi*(10.*u.pc)**2))
        DMmag = u.MagUnit(dm0, u.mag)
        m_st = 10. * u.STmag
        dm = 5. * DMmag
        M_st = m_st - dm
        assert M_st.unit.is_equivalent(u.erg/u.s/u.AA)
        assert np.abs(M_st.quantity /
                      (m_st.quantity*4.*np.pi*(100.*u.pc)**2) - 1.) < 1.e-15


class TestLogQuantityComparisons(object):
    def test_comparison_to_non_quantities_fails(self):
        lq = u.Magnitude(np.arange(1.,10.)*u.Jy)
        # on python2, ordering operations always succeed, given essentially
        # meaningless results
        if six.PY3:
            with pytest.raises(TypeError):
                lq > 'a'

        assert not (lq == 'a')
        assert lq != 'a'

    def test_comparison(self):
        lq1 = u.Magnitude(np.arange(1.,4.)*u.Jy)
        lq2 = u.Magnitude(2.*u.Jy)
        assert np.all((lq1 > lq2) == np.array([True, False, False]))
        assert np.all((lq1 == lq2) == np.array([False, True, False]))
        lq3 = u.Dex(2.*u.Jy)
        assert np.all((lq1 > lq3) == np.array([True, False, False]))
        assert np.all((lq1 == lq3) == np.array([False, True, False]))
        lq4 = u.Magnitude(2.*u.m)
        assert not (lq1 == lq4)
        assert lq1 != lq4
        with pytest.raises(u.UnitsError):
            lq1 < lq4

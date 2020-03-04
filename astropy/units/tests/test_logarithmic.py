# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
    Test the Logarithmic Units and Quantities
"""

import pickle
import itertools

import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u, constants as c

lu_units = [u.dex, u.mag, u.decibel]

lu_subclasses = [u.DexUnit, u.MagUnit, u.DecibelUnit]

lq_subclasses = [u.Dex, u.Magnitude, u.Decibel]

pu_sample = (u.dimensionless_unscaled, u.m, u.g/u.s**2, u.Jy)


class TestLogUnitCreation:

    def test_logarithmic_units(self):
        """Check logarithmic units are set up correctly."""
        assert u.dB.to(u.dex) == 0.1
        assert u.dex.to(u.mag) == -2.5
        assert u.mag.to(u.dB) == -4

    @pytest.mark.parametrize('lu_unit, lu_cls', zip(lu_units, lu_subclasses))
    def test_callable_units(self, lu_unit, lu_cls):
        assert isinstance(lu_unit, u.UnitBase)
        assert callable(lu_unit)
        assert lu_unit._function_unit_class is lu_cls

    @pytest.mark.parametrize('lu_unit', lu_units)
    def test_equality_to_normal_unit_for_dimensionless(self, lu_unit):
        lu = lu_unit()
        assert lu == lu._default_function_unit  # eg, MagUnit() == u.mag
        assert lu._default_function_unit == lu  # and u.mag == MagUnit()

    @pytest.mark.parametrize('lu_unit, physical_unit',
                             itertools.product(lu_units, pu_sample))
    def test_call_units(self, lu_unit, physical_unit):
        """Create a LogUnit subclass using the callable unit and physical unit,
        and do basic check that output is right."""
        lu1 = lu_unit(physical_unit)
        assert lu1.physical_unit == physical_unit
        assert lu1.function_unit == lu1._default_function_unit

    def test_call_invalid_unit(self):
        with pytest.raises(TypeError):
            u.mag([])
        with pytest.raises(ValueError):
            u.mag(u.mag())

    @pytest.mark.parametrize('lu_cls, physical_unit', itertools.product(
        lu_subclasses + [u.LogUnit], pu_sample))
    def test_subclass_creation(self, lu_cls, physical_unit):
        """Create a LogUnit subclass object for given physical unit,
        and do basic check that output is right."""
        lu1 = lu_cls(physical_unit)
        assert lu1.physical_unit == physical_unit
        assert lu1.function_unit == lu1._default_function_unit

        lu2 = lu_cls(physical_unit,
                     function_unit=2*lu1._default_function_unit)
        assert lu2.physical_unit == physical_unit
        assert lu2.function_unit == u.Unit(2*lu2._default_function_unit)

        with pytest.raises(ValueError):
            lu_cls(physical_unit, u.m)

    def test_lshift_magnitude(self):
        mag = 1. << u.ABmag
        assert isinstance(mag, u.Magnitude)
        assert mag.unit == u.ABmag
        assert mag.value == 1.
        # same test for an array, which should produce a view
        a2 = np.arange(10.)
        q2 = a2 << u.ABmag
        assert isinstance(q2, u.Magnitude)
        assert q2.unit == u.ABmag
        assert np.all(q2.value == a2)
        a2[9] = 0.
        assert np.all(q2.value == a2)
        # a different magnitude unit
        mag = 10. << u.STmag
        assert isinstance(mag, u.Magnitude)
        assert mag.unit == u.STmag
        assert mag.value == 10.

    def test_ilshift_magnitude(self):
        # test in-place operation and conversion
        mag_fnu_cgs = u.mag(u.erg/u.s/u.cm**2/u.Hz)
        m = np.arange(10.0) * u.mag(u.Jy)
        jy = m.physical
        m2 = m << mag_fnu_cgs
        assert np.all(m2 == m.to(mag_fnu_cgs))
        m2 = m
        m <<= mag_fnu_cgs
        assert m is m2  # Check it was done in-place!
        assert np.all(m.value == m2.value)
        assert m.unit == mag_fnu_cgs
        # Check it works if equivalencies are in-place.
        with u.add_enabled_equivalencies(u.spectral_density(5500*u.AA)):
            st = jy.to(u.ST)
            m <<= u.STmag

        assert m is m2
        assert_quantity_allclose(m.physical, st)
        assert m.unit == u.STmag

    def test_lshift_errors(self):
        m = np.arange(10.0) * u.mag(u.Jy)
        with pytest.raises(u.UnitsError):
            m << u.STmag

        with pytest.raises(u.UnitsError):
            m << u.Jy

        with pytest.raises(u.UnitsError):
            m <<= u.STmag

        with pytest.raises(u.UnitsError):
            m <<= u.Jy


def test_predefined_magnitudes():
    assert_quantity_allclose((-21.1*u.STmag).physical,
                             1.*u.erg/u.cm**2/u.s/u.AA)
    assert_quantity_allclose((-48.6*u.ABmag).physical,
                             1.*u.erg/u.cm**2/u.s/u.Hz)

    assert_quantity_allclose((0*u.M_bol).physical, c.L_bol0)
    assert_quantity_allclose((0*u.m_bol).physical,
                             c.L_bol0/(4.*np.pi*(10.*c.pc)**2))


def test_predefined_reinitialisation():
    assert u.mag('STflux') == u.STmag
    assert u.mag('ABflux') == u.ABmag
    assert u.mag('Bol') == u.M_bol
    assert u.mag('bol') == u.m_bol

    # required for backwards-compatibility, at least unless deprecated
    assert u.mag('ST') == u.STmag
    assert u.mag('AB') == u.ABmag


def test_predefined_string_roundtrip():
    """Ensure round-tripping; see #5015"""
    assert u.Unit(u.STmag.to_string()) == u.STmag
    assert u.Unit(u.ABmag.to_string()) == u.ABmag
    assert u.Unit(u.M_bol.to_string()) == u.M_bol
    assert u.Unit(u.m_bol.to_string()) == u.m_bol


def test_inequality():
    """Check __ne__ works (regresssion for #5342)."""
    lu1 = u.mag(u.Jy)
    lu2 = u.dex(u.Jy)
    lu3 = u.mag(u.Jy**2)
    lu4 = lu3 - lu1
    assert lu1 != lu2
    assert lu1 != lu3
    assert lu1 == lu4


class TestLogUnitStrings:

    def test_str(self):
        """Do some spot checks that str, repr, etc. work as expected."""
        lu1 = u.mag(u.Jy)
        assert str(lu1) == 'mag(Jy)'
        assert repr(lu1) == 'Unit("mag(Jy)")'
        assert lu1.to_string('generic') == 'mag(Jy)'
        with pytest.raises(ValueError):
            lu1.to_string('fits')
        with pytest.raises(ValueError):
            lu1.to_string(format='cds')

        lu2 = u.dex()
        assert str(lu2) == 'dex'
        assert repr(lu2) == 'Unit("dex(1)")'
        assert lu2.to_string() == 'dex(1)'

        lu3 = u.MagUnit(u.Jy, function_unit=2*u.mag)
        assert str(lu3) == '2 mag(Jy)'
        assert repr(lu3) == 'MagUnit("Jy", unit="2 mag")'
        assert lu3.to_string() == '2 mag(Jy)'

        lu4 = u.mag(u.ct)
        assert lu4.to_string('generic') == 'mag(ct)'
        assert lu4.to_string('latex') == ('$\\mathrm{mag}$$\\mathrm{\\left( '
                                          '\\mathrm{ct} \\right)}$')
        assert lu4._repr_latex_() == lu4.to_string('latex')


class TestLogUnitConversion:
    @pytest.mark.parametrize('lu_unit, physical_unit',
                             itertools.product(lu_units, pu_sample))
    def test_physical_unit_conversion(self, lu_unit, physical_unit):
        """Check various LogUnit subclasses are equivalent and convertible
        to their non-log counterparts."""
        lu1 = lu_unit(physical_unit)
        assert lu1.is_equivalent(physical_unit)
        assert lu1.to(physical_unit, 0.) == 1.

        assert physical_unit.is_equivalent(lu1)
        assert physical_unit.to(lu1, 1.) == 0.

        pu = u.Unit(8.*physical_unit)
        assert lu1.is_equivalent(physical_unit)
        assert lu1.to(pu, 0.) == 0.125

        assert pu.is_equivalent(lu1)
        assert_allclose(pu.to(lu1, 0.125), 0., atol=1.e-15)

        # Check we round-trip.
        value = np.linspace(0., 10., 6)
        assert_allclose(pu.to(lu1, lu1.to(pu, value)), value, atol=1.e-15)
        # And that we're not just returning True all the time.
        pu2 = u.g
        assert not lu1.is_equivalent(pu2)
        with pytest.raises(u.UnitsError):
            lu1.to(pu2)

        assert not pu2.is_equivalent(lu1)
        with pytest.raises(u.UnitsError):
            pu2.to(lu1)

    @pytest.mark.parametrize('lu_unit', lu_units)
    def test_container_unit_conversion(self, lu_unit):
        """Check that conversion to logarithmic units (u.mag, u.dB, u.dex)
        is only possible when the physical unit is dimensionless."""
        values = np.linspace(0., 10., 6)
        lu1 = lu_unit(u.dimensionless_unscaled)
        assert lu1.is_equivalent(lu1.function_unit)
        assert_allclose(lu1.to(lu1.function_unit, values), values)

        lu2 = lu_unit(u.Jy)
        assert not lu2.is_equivalent(lu2.function_unit)
        with pytest.raises(u.UnitsError):
            lu2.to(lu2.function_unit, values)

    @pytest.mark.parametrize(
        'flu_unit, tlu_unit, physical_unit',
        itertools.product(lu_units, lu_units, pu_sample))
    def test_subclass_conversion(self, flu_unit, tlu_unit, physical_unit):
        """Check various LogUnit subclasses are equivalent and convertible
        to each other if they correspond to equivalent physical units."""
        values = np.linspace(0., 10., 6)
        flu = flu_unit(physical_unit)

        tlu = tlu_unit(physical_unit)
        assert flu.is_equivalent(tlu)
        assert_allclose(flu.to(tlu), flu.function_unit.to(tlu.function_unit))
        assert_allclose(flu.to(tlu, values),
                        values * flu.function_unit.to(tlu.function_unit))

        tlu2 = tlu_unit(u.Unit(100.*physical_unit))
        assert flu.is_equivalent(tlu2)
        # Check that we round-trip.
        assert_allclose(flu.to(tlu2, tlu2.to(flu, values)), values, atol=1.e-15)

        tlu3 = tlu_unit(physical_unit.to_system(u.si)[0])
        assert flu.is_equivalent(tlu3)
        assert_allclose(flu.to(tlu3, tlu3.to(flu, values)), values, atol=1.e-15)

        tlu4 = tlu_unit(u.g)
        assert not flu.is_equivalent(tlu4)
        with pytest.raises(u.UnitsError):
            flu.to(tlu4, values)

    def test_unit_decomposition(self):
        lu = u.mag(u.Jy)
        assert lu.decompose() == u.mag(u.Jy.decompose())
        assert lu.decompose().physical_unit.bases == [u.kg, u.s]
        assert lu.si == u.mag(u.Jy.si)
        assert lu.si.physical_unit.bases == [u.kg, u.s]
        assert lu.cgs == u.mag(u.Jy.cgs)
        assert lu.cgs.physical_unit.bases == [u.g, u.s]

    def test_unit_multiple_possible_equivalencies(self):
        lu = u.mag(u.Jy)
        assert lu.is_equivalent(pu_sample)

    def test_magnitude_conversion_fails_message(self):
        """Check that "dimensionless" magnitude units include a message in their
        exception text suggesting a possible cause of the problem.
        """
        with pytest.raises(u.UnitConversionError) as excinfo:
            (10*u.ABmag - 2*u.ABmag).to(u.nJy)
        assert "Did you perhaps subtract magnitudes so the unit got lost?" in str(excinfo.value)


class TestLogUnitArithmetic:
    def test_multiplication_division(self):
        """Check that multiplication/division with other units is only
        possible when the physical unit is dimensionless, and that this
        turns the unit into a normal one."""
        lu1 = u.mag(u.Jy)

        with pytest.raises(u.UnitsError):
            lu1 * u.m

        with pytest.raises(u.UnitsError):
            u.m * lu1

        with pytest.raises(u.UnitsError):
            lu1 / lu1

        for unit in (u.dimensionless_unscaled, u.m, u.mag, u.dex):
            with pytest.raises(u.UnitsError):
                lu1 / unit

        lu2 = u.mag(u.dimensionless_unscaled)

        with pytest.raises(u.UnitsError):
            lu2 * lu1

        with pytest.raises(u.UnitsError):
            lu2 / lu1

        # But dimensionless_unscaled can be cancelled.
        assert lu2 / lu2 == u.dimensionless_unscaled

        # With dimensionless, normal units are OK, but we return a plain unit.
        tf = lu2 * u.m
        tr = u.m * lu2
        for t in (tf, tr):
            assert not isinstance(t, type(lu2))
            assert t == lu2.function_unit * u.m
            with u.set_enabled_equivalencies(u.logarithmic()):
                with pytest.raises(u.UnitsError):
                    t.to(lu2.physical_unit)
        # Now we essentially have a LogUnit with a prefactor of 100,
        # so should be equivalent again.
        t = tf / u.cm
        with u.set_enabled_equivalencies(u.logarithmic()):
            assert t.is_equivalent(lu2.function_unit)
            assert_allclose(t.to(u.dimensionless_unscaled, np.arange(3.)/100.),
                            lu2.to(lu2.physical_unit, np.arange(3.)))

        # If we effectively remove lu1, a normal unit should be returned.
        t2 = tf / lu2
        assert not isinstance(t2, type(lu2))
        assert t2 == u.m
        t3 = tf / lu2.function_unit
        assert not isinstance(t3, type(lu2))
        assert t3 == u.m

        # For completeness, also ensure non-sensical operations fail
        with pytest.raises(TypeError):
            lu1 * object()
        with pytest.raises(TypeError):
            slice(None) * lu1
        with pytest.raises(TypeError):
            lu1 / []
        with pytest.raises(TypeError):
            1 / lu1

    @pytest.mark.parametrize('power', (2, 0.5, 1, 0))
    def test_raise_to_power(self, power):
        """Check that raising LogUnits to some power is only possible when the
        physical unit is dimensionless, and that conversion is turned off when
        the resulting logarithmic unit (such as mag**2) is incompatible."""
        lu1 = u.mag(u.Jy)

        if power == 0:
            assert lu1 ** power == u.dimensionless_unscaled
        elif power == 1:
            assert lu1 ** power == lu1
        else:
            with pytest.raises(u.UnitsError):
                lu1 ** power

        # With dimensionless, though, it works, but returns a normal unit.
        lu2 = u.mag(u.dimensionless_unscaled)

        t = lu2**power
        if power == 0:
            assert t == u.dimensionless_unscaled
        elif power == 1:
            assert t == lu2
        else:
            assert not isinstance(t, type(lu2))
            assert t == lu2.function_unit**power
            # also check we roundtrip
            t2 = t**(1./power)
            assert t2 == lu2.function_unit
            with u.set_enabled_equivalencies(u.logarithmic()):
                assert_allclose(t2.to(u.dimensionless_unscaled, np.arange(3.)),
                                lu2.to(lu2.physical_unit, np.arange(3.)))

    @pytest.mark.parametrize('other', pu_sample)
    def test_addition_subtraction_to_normal_units_fails(self, other):
        lu1 = u.mag(u.Jy)
        with pytest.raises(u.UnitsError):
            lu1 + other

        with pytest.raises(u.UnitsError):
            lu1 - other

        with pytest.raises(u.UnitsError):
            other - lu1

    def test_addition_subtraction_to_non_units_fails(self):
        lu1 = u.mag(u.Jy)
        with pytest.raises(TypeError):
            lu1 + 1.

        with pytest.raises(TypeError):
            lu1 - [1., 2., 3.]

    @pytest.mark.parametrize(
        'other', (u.mag, u.mag(), u.mag(u.Jy), u.mag(u.m),
                  u.Unit(2*u.mag), u.MagUnit('', 2.*u.mag)))
    def test_addition_subtraction(self, other):
        """Check physical units are changed appropriately"""
        lu1 = u.mag(u.Jy)
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
        lu_dm = u.mag(dm0)
        lu_absST = u.STmag - lu_dm
        assert lu_absST.is_equivalent(u.erg/u.s/u.AA)

    def test_neg_pos(self):
        lu1 = u.mag(u.Jy)
        neg_lu = -lu1
        assert neg_lu != lu1
        assert neg_lu.physical_unit == u.Jy**-1
        assert -neg_lu == lu1
        pos_lu = +lu1
        assert pos_lu is not lu1
        assert pos_lu == lu1


def test_pickle():
    lu1 = u.dex(u.cm/u.s**2)
    s = pickle.dumps(lu1)
    lu2 = pickle.loads(s)
    assert lu1 == lu2


def test_hashable():
    lu1 = u.dB(u.mW)
    lu2 = u.dB(u.m)
    lu3 = u.dB(u.mW)
    assert hash(lu1) != hash(lu2)
    assert hash(lu1) == hash(lu3)
    luset = {lu1, lu2, lu3}
    assert len(luset) == 2


class TestLogQuantityCreation:

    @pytest.mark.parametrize('lq, lu', zip(lq_subclasses + [u.LogQuantity],
                                           lu_subclasses + [u.LogUnit]))
    def test_logarithmic_quantities(self, lq, lu):
        """Check logarithmic quantities are all set up correctly"""
        assert lq._unit_class == lu
        assert type(lu()._quantity_class(1.)) is lq

    @pytest.mark.parametrize('lq_cls, physical_unit',
                             itertools.product(lq_subclasses, pu_sample))
    def test_subclass_creation(self, lq_cls, physical_unit):
        """Create LogQuantity subclass objects for some physical units,
        and basic check on transformations"""
        value = np.arange(1., 10.)
        log_q = lq_cls(value * physical_unit)
        assert log_q.unit.physical_unit == physical_unit
        assert log_q.unit.function_unit == log_q.unit._default_function_unit
        assert_allclose(log_q.physical.value, value)
        with pytest.raises(ValueError):
            lq_cls(value, physical_unit)

    @pytest.mark.parametrize(
        'unit', (u.mag, u.mag(), u.mag(u.Jy), u.mag(u.m),
                 u.Unit(2*u.mag), u.MagUnit('', 2.*u.mag),
                 u.MagUnit(u.Jy, -1*u.mag), u.MagUnit(u.m, -2.*u.mag)))
    def test_different_units(self, unit):
        q = u.Magnitude(1.23, unit)
        assert q.unit.function_unit == getattr(unit, 'function_unit', unit)
        assert q.unit.physical_unit is getattr(unit, 'physical_unit',
                                               u.dimensionless_unscaled)

    @pytest.mark.parametrize('value, unit', (
        (1.*u.mag(u.Jy), None),
        (1.*u.dex(u.Jy), None),
        (1.*u.mag(u.W/u.m**2/u.Hz), u.mag(u.Jy)),
        (1.*u.dex(u.W/u.m**2/u.Hz), u.mag(u.Jy))))
    def test_function_values(self, value, unit):
        lq = u.Magnitude(value, unit)
        assert lq == value
        assert lq.unit.function_unit == u.mag
        assert lq.unit.physical_unit == getattr(unit, 'physical_unit',
                                                value.unit.physical_unit)

    @pytest.mark.parametrize(
        'unit', (u.mag(), u.mag(u.Jy), u.mag(u.m), u.MagUnit('', 2.*u.mag),
                 u.MagUnit(u.Jy, -1*u.mag), u.MagUnit(u.m, -2.*u.mag)))
    def test_indirect_creation(self, unit):
        q1 = 2.5 * unit
        assert isinstance(q1, u.Magnitude)
        assert q1.value == 2.5
        assert q1.unit == unit
        pv = 100. * unit.physical_unit
        q2 = unit * pv
        assert q2.unit == unit
        assert q2.unit.physical_unit == pv.unit
        assert q2.to_value(unit.physical_unit) == 100.
        assert (q2._function_view / u.mag).to_value(1) == -5.
        q3 = unit / 0.4
        assert q3 == q1

    def test_from_view(self):
        # Cannot view a physical quantity as a function quantity, since the
        # values would change.
        q = [100., 1000.] * u.cm/u.s**2
        with pytest.raises(TypeError):
            q.view(u.Dex)
        # But fine if we have the right magnitude.
        q = [2., 3.] * u.dex
        lq = q.view(u.Dex)
        assert isinstance(lq, u.Dex)
        assert lq.unit.physical_unit == u.dimensionless_unscaled
        assert np.all(q == lq)

    def test_using_quantity_class(self):
        """Check that we can use Quantity if we have subok=True"""
        # following issue #5851
        lu = u.dex(u.AA)
        with pytest.raises(u.UnitTypeError):
            u.Quantity(1., lu)
        q = u.Quantity(1., lu, subok=True)
        assert type(q) is lu._quantity_class


def test_conversion_to_and_from_physical_quantities():
    """Ensures we can convert from regular quantities."""
    mst = [10., 12., 14.] * u.STmag
    flux_lambda = mst.physical
    mst_roundtrip = flux_lambda.to(u.STmag)
    # check we return a logquantity; see #5178.
    assert isinstance(mst_roundtrip, u.Magnitude)
    assert mst_roundtrip.unit == mst.unit
    assert_allclose(mst_roundtrip.value, mst.value)
    wave = [4956.8, 4959.55, 4962.3] * u.AA
    flux_nu = mst.to(u.Jy, equivalencies=u.spectral_density(wave))
    mst_roundtrip2 = flux_nu.to(u.STmag, u.spectral_density(wave))
    assert isinstance(mst_roundtrip2, u.Magnitude)
    assert mst_roundtrip2.unit == mst.unit
    assert_allclose(mst_roundtrip2.value, mst.value)


def test_quantity_decomposition():
    lq = 10.*u.mag(u.Jy)
    assert lq.decompose() == lq
    assert lq.decompose().unit.physical_unit.bases == [u.kg, u.s]
    assert lq.si == lq
    assert lq.si.unit.physical_unit.bases == [u.kg, u.s]
    assert lq.cgs == lq
    assert lq.cgs.unit.physical_unit.bases == [u.g, u.s]


class TestLogQuantityViews:
    def setup(self):
        self.lq = u.Magnitude(np.arange(1., 10.) * u.Jy)
        self.lq2 = u.Magnitude(np.arange(1., 5.))

    def test_value_view(self):
        lq_value = self.lq.value
        assert type(lq_value) is np.ndarray
        lq_value[2] = -1.
        assert np.all(self.lq.value == lq_value)

    def test_function_view(self):
        lq_fv = self.lq._function_view
        assert type(lq_fv) is u.Quantity
        assert lq_fv.unit is self.lq.unit.function_unit
        lq_fv[3] = -2. * lq_fv.unit
        assert np.all(self.lq.value == lq_fv.value)

    def test_quantity_view(self):
        # Cannot view as Quantity, since the unit cannot be represented.
        with pytest.raises(TypeError):
            self.lq.view(u.Quantity)
        # But a dimensionless one is fine.
        q2 = self.lq2.view(u.Quantity)
        assert q2.unit is u.mag
        assert np.all(q2.value == self.lq2.value)
        lq3 = q2.view(u.Magnitude)
        assert type(lq3.unit) is u.MagUnit
        assert lq3.unit.physical_unit == u.dimensionless_unscaled
        assert np.all(lq3 == self.lq2)


class TestLogQuantitySlicing:
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


class TestLogQuantityArithmetic:
    def test_multiplication_division(self):
        """Check that multiplication/division with other quantities is only
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
            assert not isinstance(t, type(lq2))
            assert t.unit == lq2.unit.function_unit * u.m
            with u.set_enabled_equivalencies(u.logarithmic()):
                with pytest.raises(u.UnitsError):
                    t.to(lq2.unit.physical_unit)

        t = tf / (50.*u.cm)
        # now we essentially have the same quantity but with a prefactor of 2
        assert t.unit.is_equivalent(lq2.unit.function_unit)
        assert_allclose(t.to(lq2.unit.function_unit), lq2._function_view*2)

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
            assert not isinstance(t, type(lq2))
            assert t.unit == lq2.unit.function_unit ** power
            with u.set_enabled_equivalencies(u.logarithmic()):
                with pytest.raises(u.UnitsError):
                    t.to(u.dimensionless_unscaled)

    def test_error_on_lq_as_power(self):
        lq = u.Magnitude(np.arange(1., 4.)*u.Jy)
        with pytest.raises(TypeError):
            lq ** lq

    @pytest.mark.parametrize('other', pu_sample)
    def test_addition_subtraction_to_normal_units_fails(self, other):
        lq = u.Magnitude(np.arange(1., 10.)*u.Jy)
        q = 1.23 * other
        with pytest.raises(u.UnitsError):
            lq + q

        with pytest.raises(u.UnitsError):
            lq - q

        with pytest.raises(u.UnitsError):
            q - lq

    @pytest.mark.parametrize(
        'other', (1.23 * u.mag, 2.34 * u.mag(),
                  u.Magnitude(3.45 * u.Jy), u.Magnitude(4.56 * u.m),
                  5.67 * u.Unit(2*u.mag), u.Magnitude(6.78, 2.*u.mag)))
    def test_addition_subtraction(self, other):
        """Check that addition/subtraction with quantities with magnitude or
        MagUnit units works, and that it changes the physical units
        appropriately."""
        lq = u.Magnitude(np.arange(1., 10.)*u.Jy)
        other_physical = other.to(getattr(other.unit, 'physical_unit',
                                          u.dimensionless_unscaled),
                                  equivalencies=u.logarithmic())

        lq_sf = lq + other
        assert_allclose(lq_sf.physical, lq.physical * other_physical)

        lq_sr = other + lq
        assert_allclose(lq_sr.physical, lq.physical * other_physical)

        lq_df = lq - other
        assert_allclose(lq_df.physical, lq.physical / other_physical)

        lq_dr = other - lq
        assert_allclose(lq_dr.physical, other_physical / lq.physical)

    @pytest.mark.parametrize('other', pu_sample)
    def test_inplace_addition_subtraction_unit_checks(self, other):
        lu1 = u.mag(u.Jy)
        lq1 = u.Magnitude(np.arange(1., 10.), lu1)
        with pytest.raises(u.UnitsError):
            lq1 += other

        assert np.all(lq1.value == np.arange(1., 10.))
        assert lq1.unit == lu1

        with pytest.raises(u.UnitsError):
            lq1 -= other

        assert np.all(lq1.value == np.arange(1., 10.))
        assert lq1.unit == lu1

    @pytest.mark.parametrize(
        'other', (1.23 * u.mag, 2.34 * u.mag(),
                  u.Magnitude(3.45 * u.Jy), u.Magnitude(4.56 * u.m),
                  5.67 * u.Unit(2*u.mag), u.Magnitude(6.78, 2.*u.mag)))
    def test_inplace_addition_subtraction(self, other):
        """Check that inplace addition/subtraction with quantities with
        magnitude or MagUnit units works, and that it changes the physical
        units appropriately."""
        lq = u.Magnitude(np.arange(1., 10.)*u.Jy)
        other_physical = other.to(getattr(other.unit, 'physical_unit',
                                          u.dimensionless_unscaled),
                                  equivalencies=u.logarithmic())
        lq_sf = lq.copy()
        lq_sf += other
        assert_allclose(lq_sf.physical, lq.physical * other_physical)

        lq_df = lq.copy()
        lq_df -= other
        assert_allclose(lq_df.physical, lq.physical / other_physical)

    def test_complicated_addition_subtraction(self):
        """For fun, a more complicated example of addition and subtraction."""
        dm0 = u.Unit('DM', 1./(4.*np.pi*(10.*u.pc)**2))
        DMmag = u.mag(dm0)
        m_st = 10. * u.STmag
        dm = 5. * DMmag
        M_st = m_st - dm
        assert M_st.unit.is_equivalent(u.erg/u.s/u.AA)
        assert np.abs(M_st.physical /
                      (m_st.physical*4.*np.pi*(100.*u.pc)**2) - 1.) < 1.e-15


class TestLogQuantityComparisons:
    def test_comparison_to_non_quantities_fails(self):
        lq = u.Magnitude(np.arange(1., 10.)*u.Jy)
        with pytest.raises(TypeError):
            lq > 'a'

        assert not (lq == 'a')
        assert lq != 'a'

    def test_comparison(self):
        lq1 = u.Magnitude(np.arange(1., 4.)*u.Jy)
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
        q5 = 1.5 * u.Jy
        assert np.all((lq1 > q5) == np.array([True, False, False]))
        assert np.all((q5 < lq1) == np.array([True, False, False]))
        with pytest.raises(u.UnitsError):
            lq1 >= 2.*u.m
        with pytest.raises(u.UnitsError):
            lq1 <= lq1.value * u.mag
        # For physically dimensionless, we can compare with the function unit.
        lq6 = u.Magnitude(np.arange(1., 4.))
        fv6 = lq6.value * u.mag
        assert np.all(lq6 == fv6)
        # but not some arbitrary unit, of course.
        with pytest.raises(u.UnitsError):
            lq6 < 2.*u.m


class TestLogQuantityMethods:
    def setup(self):
        self.mJy = np.arange(1., 5.).reshape(2, 2) * u.mag(u.Jy)
        self.m1 = np.arange(1., 5.5, 0.5).reshape(3, 3) * u.mag()
        self.mags = (self.mJy, self.m1)

    @pytest.mark.parametrize('method', ('mean', 'min', 'max', 'round', 'trace',
                                        'std', 'var', 'ptp', 'diff', 'ediff1d'))
    def test_always_ok(self, method):
        for mag in self.mags:
            res = getattr(mag, method)()
            assert np.all(res.value ==
                          getattr(mag._function_view, method)().value)
            if method in ('std', 'ptp', 'diff', 'ediff1d'):
                assert res.unit == u.mag()
            elif method == 'var':
                assert res.unit == u.mag**2
            else:
                assert res.unit == mag.unit

    def test_clip(self):
        for mag in self.mags:
            assert np.all(mag.clip(2. * mag.unit, 4. * mag.unit).value ==
                          mag.value.clip(2., 4.))

    @pytest.mark.parametrize('method', ('sum', 'cumsum', 'nansum'))
    def test_only_ok_if_dimensionless(self, method):
        res = getattr(self.m1, method)()
        assert np.all(res.value ==
                      getattr(self.m1._function_view, method)().value)
        assert res.unit == self.m1.unit
        with pytest.raises(TypeError):
            getattr(self.mJy, method)()

    def test_dot(self):
        assert np.all(self.m1.dot(self.m1).value ==
                      self.m1.value.dot(self.m1.value))

    @pytest.mark.parametrize('method', ('prod', 'cumprod'))
    def test_never_ok(self, method):
        with pytest.raises(TypeError):
            getattr(self.mJy, method)()
        with pytest.raises(TypeError):
            getattr(self.m1, method)()

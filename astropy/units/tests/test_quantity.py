# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
    Test the Quantity class and related.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing.utils import assert_allclose

from ...tests.helper import pytest
from ...tests.helper import raises, pytest
from ...utils import isiterable
from ... import units as u
from ...extern.six.moves import xrange
from ...extern import six

""" The Quantity class will represent a number + unit + uncertainty """


class TestQuantityCreation(object):

    def test_1(self):
        # create objects through operations with Unit objects:

        quantity = 11.42 * u.meter  # returns a Quantity object
        assert isinstance(quantity, u.Quantity)
        quantity = u.meter * 11.42  # returns a Quantity object
        assert isinstance(quantity, u.Quantity)

        quantity = 11.42 / u.meter
        assert isinstance(quantity, u.Quantity)
        quantity = u.meter / 11.42
        assert isinstance(quantity, u.Quantity)

        quantity = 11.42 * u.meter / u.second
        assert isinstance(quantity, u.Quantity)

        with pytest.raises(TypeError):
            quantity = 182.234 + u.meter

        with pytest.raises(TypeError):
            quantity = 182.234 - u.meter

        with pytest.raises(TypeError):
            quantity = 182.234 % u.meter

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_2(self, QuantityClass):

        # create objects using the Quantity constructor:
        q1 = QuantityClass(11.412, unit=u.meter)
        q2 = QuantityClass(21.52, "cm")
        q3 = QuantityClass(11.412)

        assert q1.unit == u.meter
        assert q2.unit == u.cm
        # By default quantities that don't specify a unit are unscaled
        # dimensionless
        assert q3.unit == u.Unit(1)

        with pytest.raises(TypeError):
            QuantityClass(object(), unit=u.m)

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_3(self, QuantityClass):
        # with pytest.raises(u.UnitsError):
        with pytest.raises(ValueError):  # Until @mdboom fixes the errors in units
            QuantityClass(11.412, unit="testingggg")

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_unit_property(self, QuantityClass):
        # test getting and setting 'unit' attribute
        q1 = QuantityClass(11.4, unit=u.meter)

        with pytest.raises(AttributeError):
            q1.unit = u.cm

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_preserve_dtype(self, QuantityClass):

        # If unit is not sepcified, preserve dtype (at least to the extent
        # that Numpy does when copying, i.e. int32 -> int64, not float64)

        q1 = QuantityClass(12, unit=u.m / u.s, dtype=int)
        q2 = QuantityClass(q1)

        assert q1.value == q2.value
        assert q1.unit == q2.unit
        assert q1.dtype == q2.dtype


class TestQuantityOperations(object):
    q1 = u.Quantity(11.42, u.meter)
    q2 = u.Quantity(8.0, u.centimeter)

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_addition_and_subtraction(self, QuantityClass):
        q1 = QuantityClass(11.42, u.meter)
        q2 = QuantityClass(8.0, u.centimeter)
        # Take units from left object, q1
        new_quantity = q1 + q2
        assert new_quantity.value == 11.5
        assert new_quantity.unit == u.meter

        # Take units from left object, q2
        new_quantity = q2 + q1
        assert new_quantity.value == 1150.0
        assert new_quantity.unit == u.centimeter

        new_q = QuantityClass(1500.1, u.m) + QuantityClass(13.5, u.km)
        assert new_q.unit == u.m
        assert new_q.value == 15000.1

        # Take units from left object, q1
        new_quantity = q1 - q2
        assert new_quantity.value == 11.34
        assert new_quantity.unit == u.meter

        # Take units from left object, q2
        new_quantity = q2 - q1
        assert new_quantity.value == -1134.0
        assert new_quantity.unit == u.centimeter

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_multiplication_and_division(self, QuantityClass):
        q1 = QuantityClass(11.42, u.meter)
        q2 = QuantityClass(8.0, u.centimeter)
        # Take units from left object, q1
        new_quantity = q1 * q2
        assert new_quantity.value == 91.36
        assert new_quantity.unit == (u.meter * u.centimeter)

        # Take units from left object, q2
        new_quantity = q2 * q1
        assert new_quantity.value == 91.36
        assert new_quantity.unit == (u.centimeter * u.meter)

        # Multiply with a number
        new_quantity = 15. * q1
        assert new_quantity.value == 171.3
        assert new_quantity.unit == u.meter

        # Multiply with a number
        new_quantity = q1 * 15.
        assert new_quantity.value == 171.3
        assert new_quantity.unit == u.meter

        # Take units from left object, q1
        new_quantity = q1 / q2
        np.testing.assert_array_almost_equal(
            new_quantity.value, 1.4275, decimal=5)
        assert new_quantity.unit == (u.meter / u.centimeter)

        # Take units from left object, q2
        new_quantity = q2 / q1
        np.testing.assert_array_almost_equal(
            new_quantity.value, 0.70052539404553416, decimal=16)
        assert new_quantity.unit == (u.centimeter / u.meter)

        q3 = QuantityClass(11.4, unit=u.meter)
        q4 = QuantityClass(10.0, unit=u.second)
        new_quantity = q3 / q4
        np.testing.assert_array_almost_equal(
            new_quantity.value, 1.14, decimal=10)
        assert new_quantity.unit == (u.meter / u.second)

        # divide with a number
        new_quantity = q1 / 10.
        np.testing.assert_array_almost_equal(
            new_quantity.value, 1.142, decimal=10)
        assert new_quantity.unit == u.meter

        # divide with a number
        new_quantity = 11.42 / q1
        assert new_quantity.value == 1.
        assert new_quantity.unit == u.Unit("1/m")

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_commutativity(self, QuantityClass):
        """Regression test for issue #587."""

        q1 = QuantityClass(11.42, u.meter)
        new_q = QuantityClass(11.42, 'm*s')

        assert q1 * u.s == u.s * q1 == new_q
        assert q1 / u.s == QuantityClass(11.42, 'm/s')
        if QuantityClass is u.MaskedQuantity:
            pytest.xfail()
        assert u.s / q1 == QuantityClass(1 / 11.42, 's/m')

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_power(self, QuantityClass):
        q1 = QuantityClass(11.42, u.meter)
        # raise quantity to a power
        new_quantity = q1 ** 2
        np.testing.assert_array_almost_equal(
            new_quantity.value, 130.4164, decimal=5)
        assert new_quantity.unit == u.Unit("m^2")

        new_quantity = q1 ** 3
        np.testing.assert_array_almost_equal(
            new_quantity.value, 1489.355288, decimal=7)
        assert new_quantity.unit == u.Unit("m^3")

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_unary_and_abs(self, QuantityClass):
        q1 = QuantityClass(11.42, u.meter)

        # Test the minus unary operator

        new_quantity = -q1
        assert new_quantity.value == -q1.value
        assert new_quantity.unit == q1.unit

        new_quantity = -(-q1)
        assert new_quantity.value == q1.value
        assert new_quantity.unit == q1.unit

        # Test the plus unary operator

        new_quantity = +q1
        assert new_quantity.value == q1.value
        assert new_quantity.unit == q1.unit

        new_quantity = abs(q1)
        assert new_quantity.value == q1.value
        assert new_quantity.unit == q1.unit

        new_quantity = abs(-q1)
        assert new_quantity.value == q1.value
        assert new_quantity.unit == q1.unit

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_incompatible_units(self, QuantityClass):
        """ When trying to add or subtract units that aren't compatible,
        throw an error """

        q1 = QuantityClass(11.412, unit=u.meter)
        q2 = QuantityClass(21.52, unit=u.second)

        with pytest.raises(u.UnitsError):
            q1 + q2

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_dimensionless_operations(self, QuantityClass):
        # test conversion to dimensionless
        dq = QuantityClass(3., u.m / u.km)
        dq1 = dq + QuantityClass(1., u.mm / u.km)
        assert dq1.value == 3.001
        assert dq1.unit == dq.unit

        dq2 = dq + 1.
        assert dq2.value == 1.003
        assert dq2.unit == u.dimensionless_unscaled

        # this test will check that operations with dimensionless Quantities
        # don't work
        with pytest.raises(u.UnitsError):
            QuantityClass(1., "m") + QuantityClass(0.1, "")

        with pytest.raises(u.UnitsError):
            QuantityClass(1., "m") - QuantityClass(0.1, "")

        # and test that scaling of integers works
        q = QuantityClass(np.array([1, 2, 3]), u.m / u.km)
        q2 = q + np.array([4, 5, 6])
        assert q2.unit == u.dimensionless_unscaled
        assert_allclose(q2.value, np.array([4.001, 5.002, 6.003]))
        # but not if doing it inplace
        with pytest.raises(TypeError):
            q += np.array([1, 2, 3])
        # except if it is actually possible
        q = QuantityClass(np.array([1, 2, 3]), u.km / u.m)
        q += np.array([4, 5, 6])
        assert q.unit == u.dimensionless_unscaled
        assert np.all(q.value == np.array([1004, 2005, 3006]))

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_complicated_operation(self, QuantityClass):
        """ Perform a more complicated test """
        from .. import imperial

        # Multiple units
        distance = QuantityClass(15., u.meter)
        time = QuantityClass(11., u.second)

        velocity = (distance / time).to(imperial.mile / u.hour)
        np.testing.assert_array_almost_equal(
            velocity.value, 3.05037, decimal=5)

        G = QuantityClass(6.673E-11, u.m ** 3 / u.kg / u.s ** 2)
        new_q = ((1. / (4. * np.pi * G)).to(u.pc ** -3 / u.s ** -2 * u.kg))
        np.testing.assert_array_almost_equal(
            (1./(4.*np.pi*new_q)).to(G.unit).value, G.value, decimal=9)

        # Area
        side1 = QuantityClass(11., u.centimeter)
        side2 = QuantityClass(7., u.centimeter)
        area = side1 * side2
        np.testing.assert_array_almost_equal(area.value, 77., decimal=15)
        assert area.unit == u.cm * u.cm

    def test_comparison(self):
        # equality/ non-equality is straightforward for quantity objects
        assert (1 / (u.cm * u.cm)) == 1 * u.cm ** -2
        assert 1 * u.m == 100 * u.cm
        assert 1 * u.m != 1 * u.cm

        # here one is a unit, which is an invalid comparison
        assert 1. * u.cm * u.cm * u.cm != u.cm ** 3

        # mismatched types should never work
        assert not 1. * u.cm == 1.
        assert 1. * u.cm != 1.

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_converters(self, QuantityClass):

        q1 = u.Quantity(1, u.m)

        converter_err_msg = ("Only dimensionless scalar quantities "
                             "can be converted to Python scalars")
        index_err_msg = ("Only integer dimensionless scalar quantities "
                         "can be converted to a Python index")
        with pytest.raises(TypeError) as exc:
            float(q1)
        assert exc.value.args[0] == converter_err_msg

        with pytest.raises(TypeError) as exc:
            int(q1)
        assert exc.value.args[0] == converter_err_msg

        if six.PY2:
            with pytest.raises(TypeError) as exc:
                long(q1)
            assert exc.value.args[0] == converter_err_msg

        with pytest.raises(TypeError) as exc:
            q1 * ['a', 'b', 'c']
        assert exc.value.args[0] == index_err_msg

        # dimensionless but scaled is also not OK
        q2 = QuantityClass(1.23, u.m / u.km)

        with pytest.raises(TypeError) as exc:
            float(q2)
        assert exc.value.args[0] == converter_err_msg

        with pytest.raises(TypeError) as exc:
            int(q2)
        assert exc.value.args[0] == converter_err_msg

        if six.PY2:
            with pytest.raises(TypeError) as exc:
                long(q2)
            assert exc.value.args[0] == converter_err_msg

        with pytest.raises(TypeError) as exc:
            q2 * ['a', 'b', 'c']
        assert exc.value.args[0] == index_err_msg

        # dimensionless unscaled is OK, though for index needs to be int
        q3 = QuantityClass(1.23, u.dimensionless_unscaled)

        assert float(q3) == 1.23
        assert int(q3) == 1
        if six.PY2:
            assert long(q3) == 1

        with pytest.raises(TypeError) as exc:
            q1 * ['a', 'b', 'c']
        assert exc.value.args[0] == index_err_msg

        # integer dimensionless unscaled is good for all
        q4 = QuantityClass(2, u.dimensionless_unscaled)

        assert float(q4) == 2.
        assert int(q4) == 2
        if six.PY2:
            assert long(q4) == 2

        assert q4 * ['a', 'b', 'c'] == ['a', 'b', 'c', 'a', 'b', 'c']

        # but arrays are not OK
        q5 = QuantityClass([1, 2], u.m)
        with pytest.raises(TypeError) as exc:
            float(q5)
        assert exc.value.args[0] == converter_err_msg

        with pytest.raises(TypeError) as exc:
            int(q5)
        assert exc.value.args[0] == converter_err_msg

        if six.PY2:
            with pytest.raises(TypeError) as exc:
                long(q5)
            assert exc.value.args[0] == converter_err_msg

        with pytest.raises(TypeError) as exc:
            q5 * ['a', 'b', 'c']
        assert exc.value.args[0] == index_err_msg

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_array_converters(self, QuantityClass):
        # Scalar quantity
        q = QuantityClass(1.23, u.m)
        assert np.all(np.array(q) == np.array([1.23]))

        # Array quantity
        q = QuantityClass([1., 2., 3.], u.m)
        assert np.all(np.array(q) == np.array([1., 2., 3.]))


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_conversion(QuantityClass):
    q1 = QuantityClass(0.1, unit=u.meter)

    new_quantity = q1.to(u.kilometer)
    assert new_quantity.value == 0.0001

    with pytest.raises(u.UnitsError):
        q1.to(u.zettastokes)


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_conversion_with_equiv(QuantityClass):
    q1 = QuantityClass(0.1, unit=u.meter)
    q2 = q1.to(u.Hz, equivalencies=u.spectral())
    assert_allclose(q2.value, 2997924580.0)

    q1 = QuantityClass(0.4, unit=u.arcsecond)
    q2 = q1.to(u.au, equivalencies=u.parallax())
    q3 = q2.to(u.arcminute, equivalencies=u.parallax())

    assert_allclose(q2.value, 515662.015)
    assert q2.unit == u.au
    assert_allclose(q3.value, 0.0066666667)
    assert q3.unit == u.arcminute


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_conversion_equivalency_passed_on(QuantityClass):
    q1 = QuantityClass([1000,2000], unit=u.Hz, equivalencies=u.spectral())
    q2 = q1.to(u.nm).to(u.Hz)
    assert q2.unit == u.Hz
    assert_allclose(q2.value, q1.value)
    q3 = QuantityClass([1000, 2000], unit=u.nm)
    q4 = q3.to(u.Hz, equivalencies=u.spectral()).to(u.nm)
    assert q4.unit == u.nm
    assert_allclose(q4.value, q3.value)


def test_si():
    q1 = 10. * u.m * u.s ** 2 / (200. * u.ms) ** 2  # 250 meters
    assert q1.si.value == 250
    assert q1.si.unit == u.m

    q = 10. * u.m  # 10 meters
    assert q.si.value == 10
    assert q.si.unit == u.m

    q = 10. / u.m  # 10 1 / meters
    assert q.si.value == 10
    assert q.si.unit == (1 / u.m)


def test_cgs():
    q1 = 10. * u.cm * u.s ** 2 / (200. * u.ms) ** 2  # 250 centimeters
    assert q1.cgs.value == 250
    assert q1.cgs.unit == u.cm

    q = 10. * u.m  # 10 centimeters
    assert q.cgs.value == 1000
    assert q.cgs.unit == u.cm

    q = 10. / u.cm  # 10 1 / centimeters
    assert q.cgs.value == 10
    assert q.cgs.unit == (1 / u.cm)

    q = 10. * u.Pa  # 10 pascals
    assert q.cgs.value == 100
    assert q.cgs.unit == u.barye


class TestQuantityComparison(object):

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_quantity_equality(self, QuantityClass):
        assert QuantityClass(1000, unit='m') == QuantityClass(1, unit='km')
        assert not (QuantityClass(1, unit='m') == QuantityClass(1, unit='km'))
        with pytest.raises(u.UnitsError):
            QuantityClass(1, unit='m') == QuantityClass(1, unit='s')

    @pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
    def test_quantity_comparison(self, QuantityClass):
        assert QuantityClass(1100, unit=u.meter) > QuantityClass(1, unit=u.kilometer)
        assert QuantityClass(900, unit=u.meter) < QuantityClass(1, unit=u.kilometer)

        with pytest.raises(u.UnitsError):
            assert QuantityClass(1100, unit=u.meter) > QuantityClass(1, unit=u.second)

        with pytest.raises(u.UnitsError):
            assert QuantityClass(1100, unit=u.meter) < QuantityClass(1, unit=u.second)

        assert QuantityClass(1100, unit=u.meter) >= QuantityClass(1, unit=u.kilometer)
        assert QuantityClass(1000, unit=u.meter) >= QuantityClass(1, unit=u.kilometer)

        assert QuantityClass(900, unit=u.meter) <= QuantityClass(1, unit=u.kilometer)
        assert QuantityClass(1000, unit=u.meter) <= QuantityClass(1, unit=u.kilometer)

        with pytest.raises(u.UnitsError):
            assert QuantityClass(1100, unit=u.meter) >= QuantityClass(1, unit=u.second)

        with pytest.raises(u.UnitsError):
            assert QuantityClass(1100, unit=u.meter) <= QuantityClass(1, unit=u.second)

        assert QuantityClass(1200, unit=u.meter) != QuantityClass(1, unit=u.kilometer)

        with pytest.raises(u.UnitsError):
            assert QuantityClass(1100, unit=u.meter) != QuantityClass(1, unit=u.second)


class TestQuantityDisplay(object):
    scalarintq = u.Quantity(1, unit='m')
    scalarfloatq = u.Quantity(1.3, unit='m')
    arrq = u.Quantity([1, 2.3, 8.9], unit='m')

    def test_scalar_quantity_str(self):
        assert str(self.scalarintq) == "1 m"
        assert str(self.scalarfloatq) == "1.3 m"

    def test_scalar_quantity_repr(self):
        assert repr(self.scalarintq) == "<Quantity 1 m>"
        assert repr(self.scalarfloatq) == "<Quantity 1.3 m>"

    def test_array_quantity_str(self):
        assert str(self.arrq) == "[ 1.   2.3  8.9] m"

    def test_array_quantity_repr(self):
        assert repr(self.arrq) == "<Quantity [ 1. , 2.3, 8.9] m>"
        
    def test_scalar_quantity_format(self):
        assert format(self.scalarintq, '02d') == "01 m"
        assert format(self.scalarfloatq, '.1f') == "1.3 m"
        assert format(self.scalarfloatq, '.0f') == "1 m"


class TestMaskedQuantityDisplay(object):
    scalarintq = u.MaskedQuantity(1, unit='m')
    scalarfloatq = u.MaskedQuantity(1.3, unit='m')
    arrq = u.MaskedQuantity([1, 2.3, 8.9], unit='m', mask=[0, 1, 0])

    def test_scalar_quantity_str(self):
        assert str(self.scalarintq) == "1 m"
        assert str(self.scalarfloatq) == "1.3 m"

    def test_scalar_quantity_repr(self):
        assert (repr(self.scalarintq) ==
                "masked_Quantity(data = 1 m,\n"
                "                mask = False,\n"
                "          fill_value = 999999)\n")
        assert (repr(self.scalarfloatq) ==
                "masked_Quantity(data = 1.3 m,\n"
                "                mask = False,\n"
                "          fill_value = 1e+20)\n")

    def test_array_quantity_str(self):
        assert str(self.arrq) == "[1.0 -- 8.9] m"

    def test_array_quantity_repr(self):
        assert (repr(self.arrq) ==
                "masked_Quantity(data = [1.0 -- 8.9] m,\n"
                "                mask = [False  True False],\n"
                "          fill_value = 1e+20)\n")


def test_decompose():
    q1 = 5 * u.N
    assert q1.decompose() == (5 * u.kg * u.m * u.s ** -2)


def test_decompose_regression():
    """
    Regression test for bug #1163

    If decompose was called multiple times on a Quantity with an array and a
    scale != 1, the result changed every time. This is because the value was
    being referenced not copied, then modified, which changed the original
    value.
    """
    q = np.array([1, 2, 3]) * u.m / (2. * u.km)
    assert np.all(q.decompose().value == np.array([0.0005, 0.001, 0.0015]))
    assert np.all(q == np.array([1, 2, 3]) * u.m / (2. * u.km))
    assert np.all(q.decompose().value == np.array([0.0005, 0.001, 0.0015]))


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_arrays(QuantityClass):
    """
    Test using quantites with array values
    """
    from numpy.testing import assert_array_equal

    qsec = QuantityClass(np.arange(10), u.second)
    assert isinstance(qsec.value, np.ndarray)
    assert not qsec.isscalar

    # len and indexing should work for arrays
    assert len(qsec) == len(qsec.value)
    qsecsub25 = qsec[2:5]
    assert qsecsub25.unit == qsec.unit
    assert isinstance(qsecsub25, QuantityClass)
    assert len(qsecsub25) == 3

    # make sure isscalar, len, and indexing behave correcly for non-arrays.
    qsecnotarray = QuantityClass(10., u.second)
    assert qsecnotarray.isscalar
    with pytest.raises(TypeError):
        len(qsecnotarray)
    with pytest.raises(TypeError):
        qsecnotarray[0]

    qseclen0array = QuantityClass(np.array(10), u.second)
    # 0d numpy array should act basically like a scalar
    assert qseclen0array.isscalar
    with pytest.raises(TypeError):
        len(qseclen0array)
    with pytest.raises(TypeError):
        qseclen0array[0]
    assert isinstance(qseclen0array.value, int)

    # can also create from lists, will auto-convert to arrays
    qsec = QuantityClass(list(xrange(10)), u.second)
    assert isinstance(qsec.value, np.ndarray)

    # quantity math should work with arrays
    assert_array_equal((qsec * 2).value, (np.arange(10) * 2))
    assert_array_equal((qsec / 2).value, (np.arange(10) / 2))
    # quantity addition/subtraction should *not* work with arrays b/c unit
    # ambiguous
    with pytest.raises(u.UnitsError):
        assert_array_equal((qsec + 2).value, (np.arange(10) + 2))
    with pytest.raises(u.UnitsError):
        assert_array_equal((qsec - 2).value, (np.arange(10) + 2))

    value = qsec.value
    assert not isinstance(value, QuantityClass)
    # should create by unit multiplication, too
    if QuantityClass is u.MaskedQuantity:
        pytest.xfail()
    qsec2 = value * u.second
    assert np.all(qsec == qsec2)

    qsec3 = u.second * value
    assert np.all(qsec == qsec3)

    # make sure numerical-converters fail when arrays are present
    with pytest.raises(TypeError):
        float(qsec)
    with pytest.raises(TypeError):
        int(qsec)
    if six.PY2:
        with pytest.raises(TypeError):
            long(qsec)


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_array_indexing_slicing(QuantityClass):
    q = QuantityClass(np.array([1., 2., 3.]), u.m)
    assert q[0] == QuantityClass(1., u.m)
    assert np.all(q[0:2] == QuantityClass([1., 2.], u.m))

    q[1:2] = np.array([400.]) * u.cm
    assert np.all(q == QuantityClass(np.array([1., 4., 3.]), u.m))


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_inverse_quantity(QuantityClass):
    """
    Regression test from issue #679
    """
    q = QuantityClass(4., u.meter / u.second)
    qot = q / 2
    toq = 2 / q
    npqot = q / np.array(2)

    assert npqot.value == 2.0
    assert npqot.unit == (u.meter / u.second)

    assert qot.value == 2.0
    assert qot.unit == (u.meter / u.second)

    assert toq.value == 0.5
    assert toq.unit == (u.second / u.meter)


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_mutability(QuantityClass):
    q = QuantityClass(9.8, u.meter / u.second / u.second)

    with pytest.raises(AttributeError):
        q.value = 3

    with pytest.raises(AttributeError):
        q.unit = u.kg


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_initialized_with_quantity(QuantityClass):
    q1 = QuantityClass(60, u.second)

    q2 = QuantityClass(q1, u.minute)
    assert q2.value == 1

    q3 = QuantityClass([q1, q2], u.second)
    assert q3[0].value == 60
    assert q3[1].value == 60

    q4 = QuantityClass([q2, q1])
    assert q4.unit == q2.unit
    assert q4[0].value == 1
    assert q4[1].value == 1


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_string_unit(QuantityClass):
    q1 = QuantityClass(1., u.m) / 's'
    assert q1.value == 1
    assert q1.unit == (u.m / u.s)

    q2 = q1 * "m"
    assert q2.unit == ((u.m * u.m) / u.s)


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_implicit_conversion(QuantityClass):
    q = QuantityClass(1.0, u.meter)
    # Manually turn this on to simulate what might happen in a subclass
    q._include_easy_conversion_members = True
    assert_allclose(q.centimeter, 100)
    assert_allclose(q.cm, 100)
    assert_allclose(q.parsec, 3.240779289469756e-17)


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_implicit_conversion_autocomplete(QuantityClass):
    q = QuantityClass(1.0, u.meter)
    # Manually turn this on to simulate what might happen in a subclass
    q._include_easy_conversion_members = True
    q.foo = 42

    attrs = dir(q)
    assert 'centimeter' in attrs
    assert 'cm' in attrs
    assert 'parsec' in attrs
    assert 'foo' in attrs
    assert 'to' in attrs
    assert 'value' in attrs
    # Something from the base class, object
    assert '__setattr__' in attrs

    with pytest.raises(AttributeError):
        q.l


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_iterability(QuantityClass):
    """Regressiont est for issue #878.

    Scalar quantities should not be iterable and should raise a type error on
    iteration.
    """

    q1 = QuantityClass([15.0, 17.0], u.m)
    assert isiterable(q1)

    q2 = six.next(iter(q1))
    assert q2 == 15.0 * u.m
    assert not isiterable(q2)
    pytest.raises(TypeError, iter, q2)


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_copy(QuantityClass):

    q1 = QuantityClass(np.array([1., 2., 3.]), unit=u.m)
    q2 = q1.copy()
    assert q2.__class__ is QuantityClass
    assert np.all(q1.value == q2.value)
    assert q1.unit == q2.unit
    assert q1.dtype == q2.dtype

    assert q1.value is not q2.value


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_deepcopy(QuantityClass):
    from copy import deepcopy

    q1 = QuantityClass(np.array([1., 2., 3.]), unit=u.m)
    q2 = deepcopy(q1)

    assert q2.__class__ is QuantityClass
    assert np.all(q1.value == q2.value)
    assert q1.unit == q2.unit
    assert q1.dtype == q2.dtype

    assert q1.value is not q2.value


def test_equality_numpy_scalar():
    """
    A regression test to ensure that numpy scalars are correctly compared
    (which originally failed due to the lack of ``__array_priority__``).
    """
    assert 10 != 10. * u.m
    assert np.int64(10) != 10 * u.m
    assert 10 * u.m != np.int64(10)


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_pickelability(QuantityClass):
    """
    Testing pickleability of quantity
    """
    from ...extern.six.moves import cPickle as pickle

    q1 = QuantityClass(np.arange(10), u.m)

    q2 = pickle.loads(pickle.dumps(q1))

    assert np.all(q1.value == q2.value)
    assert q1.unit.is_equivalent(q2.unit)
    assert q1.unit == q2.unit


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_quantity_from_string(QuantityClass):
    with pytest.raises(TypeError):
        QuantityClass(u.m * "5")
        # the reverse should also fail once #1408 is in

    with pytest.raises(TypeError):
        QuantityClass('5', u.m)

    with pytest.raises(TypeError):
        QuantityClass(['5'], u.m)

    with pytest.raises(TypeError):
        QuantityClass(np.array(['5']), u.m)


@pytest.mark.parametrize('QuantityClass', (u.Quantity, u.MaskedQuantity))
def test_unsupported(QuantityClass):
    q1 = QuantityClass(np.arange(10), u.m)

    with pytest.raises(TypeError):
        np.bitwise_and(q1, q1)

# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
    Test the Quantity class and related.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import copy

import numpy as np
from numpy.testing import (assert_allclose, assert_array_equal,
                           assert_array_almost_equal)

from ...tests.helper import raises, pytest
from ...utils import isiterable
from ... import units as u
from ...extern.six.moves import xrange
from ...extern.six.moves import cPickle as pickle
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

    def test_2(self):

        # create objects using the Quantity constructor:
        q1 = u.Quantity(11.412, unit=u.meter)
        q2 = u.Quantity(21.52, "cm")
        q3 = u.Quantity(11.412)

        # By default quantities that don't specify a unit are unscaled
        # dimensionless
        assert q3.unit == u.Unit(1)

        with pytest.raises(TypeError):
            q4 = u.Quantity(object(), unit=u.m)

    def test_3(self):
        # with pytest.raises(u.UnitsError):
        with pytest.raises(ValueError):  # Until @mdboom fixes the errors in units
            q1 = u.Quantity(11.412, unit="testingggg")

    def test_unit_property(self):
        # test getting and setting 'unit' attribute
        q1 = u.Quantity(11.4, unit=u.meter)

        with pytest.raises(AttributeError):
            q1.unit = u.cm

    def test_preserve_dtype(self):

        # If unit is not sepcified, preserve dtype (at least to the extent
        # that Numpy does when copying, i.e. int32 -> int64, not float64)

        q1 = u.Quantity(12, unit=u.m / u.s, dtype=int)
        q2 = u.Quantity(q1)

        assert q1.value == q2.value
        assert q1.unit == q2.unit
        assert q1.dtype == q2.dtype

    def test_copy(self):

        # By default, a new quantity is constructed, but not if copy=False

        a = np.arange(10.)

        q0 = u.Quantity(a, unit=u.m / u.s)
        assert q0.base is not a

        q1 = u.Quantity(a, unit=u.m / u.s, copy=False)
        assert q1.base is a

        q2 = u.Quantity(q0)
        assert q2 is not q0
        assert q2.base is not q0.base

        q2 = u.Quantity(q0, copy=False)
        assert q2 is q0
        assert q2.base is q0.base

        q3 = u.Quantity(q0, q0.unit, copy=False)
        assert q3 is q0
        assert q3.base is q0.base

        q4 = u.Quantity(q0, u.cm / u.s, copy=False)
        assert q4 is not q0
        assert q4.base is not q0.base


class TestQuantityOperations(object):
    q1 = u.Quantity(11.42, u.meter)
    q2 = u.Quantity(8.0, u.centimeter)

    def test_addition(self):
        # Take units from left object, q1
        new_quantity = self.q1 + self.q2
        assert new_quantity.value == 11.5
        assert new_quantity.unit == u.meter

        # Take units from left object, q2
        new_quantity = self.q2 + self.q1
        assert new_quantity.value == 1150.0
        assert new_quantity.unit == u.centimeter

        new_q = u.Quantity(1500.1, u.m) + u.Quantity(13.5, u.km)
        assert new_q.unit == u.m
        assert new_q.value == 15000.1

    def test_subtraction(self):
        # Take units from left object, q1
        new_quantity = self.q1 - self.q2
        assert new_quantity.value == 11.34
        assert new_quantity.unit == u.meter

        # Take units from left object, q2
        new_quantity = self.q2 - self.q1
        assert new_quantity.value == -1134.0
        assert new_quantity.unit == u.centimeter

    def test_multiplication(self):
        # Take units from left object, q1
        new_quantity = self.q1 * self.q2
        assert new_quantity.value == 91.36
        assert new_quantity.unit == (u.meter * u.centimeter)

        # Take units from left object, q2
        new_quantity = self.q2 * self.q1
        assert new_quantity.value == 91.36
        assert new_quantity.unit == (u.centimeter * u.meter)

        # Multiply with a number
        new_quantity = 15. * self.q1
        assert new_quantity.value == 171.3
        assert new_quantity.unit == u.meter

        # Multiply with a number
        new_quantity = self.q1 * 15.
        assert new_quantity.value == 171.3
        assert new_quantity.unit == u.meter

    def test_division(self):
        # Take units from left object, q1
        new_quantity = self.q1 / self.q2
        assert_array_almost_equal(new_quantity.value, 1.4275, decimal=5)
        assert new_quantity.unit == (u.meter / u.centimeter)

        # Take units from left object, q2
        new_quantity = self.q2 / self.q1
        assert_array_almost_equal(new_quantity.value, 0.70052539404553416,
                                  decimal=16)
        assert new_quantity.unit == (u.centimeter / u.meter)

        q1 = u.Quantity(11.4, unit=u.meter)
        q2 = u.Quantity(10.0, unit=u.second)
        new_quantity = q1 / q2
        assert_array_almost_equal(new_quantity.value, 1.14, decimal=10)
        assert new_quantity.unit == (u.meter / u.second)

        # divide with a number
        new_quantity = self.q1 / 10.
        assert new_quantity.value == 1.142
        assert new_quantity.unit == u.meter

        # divide with a number
        new_quantity = 11.42 / self.q1
        assert new_quantity.value == 1.
        assert new_quantity.unit == u.Unit("1/m")

    def test_commutativity(self):
        """Regression test for issue #587."""

        new_q = u.Quantity(11.42, 'm*s')

        assert self.q1 * u.s == u.s * self.q1 == new_q
        assert self.q1 / u.s == u.Quantity(11.42, 'm/s')
        assert u.s / self.q1 == u.Quantity(1 / 11.42, 's/m')

    def test_power(self):
        # raise quantity to a power
        new_quantity = self.q1 ** 2
        assert_array_almost_equal(new_quantity.value, 130.4164, decimal=5)
        assert new_quantity.unit == u.Unit("m^2")

        new_quantity = self.q1 ** 3
        assert_array_almost_equal(new_quantity.value, 1489.355288, decimal=7)
        assert new_quantity.unit == u.Unit("m^3")

    def test_unary(self):

        # Test the minus unary operator

        new_quantity = -self.q1
        assert new_quantity.value == -self.q1.value
        assert new_quantity.unit == self.q1.unit

        new_quantity = -(-self.q1)
        assert new_quantity.value == self.q1.value
        assert new_quantity.unit == self.q1.unit

        # Test the plus unary operator

        new_quantity = +self.q1
        assert new_quantity.value == self.q1.value
        assert new_quantity.unit == self.q1.unit

    def test_abs(self):

        q = 1. * u.m / u.s
        new_quantity = abs(q)
        assert new_quantity.value == q.value
        assert new_quantity.unit == q.unit

        q = -1. * u.m / u.s
        new_quantity = abs(q)
        assert new_quantity.value == -q.value
        assert new_quantity.unit == q.unit

    def test_incompatible_units(self):
        """ When trying to add or subtract units that aren't compatible, throw an error """

        q1 = u.Quantity(11.412, unit=u.meter)
        q2 = u.Quantity(21.52, unit=u.second)

        with pytest.raises(u.UnitsError):
            new_q = q1 + q2

    def test_dimensionless_operations(self):
        # test conversion to dimensionless
        dq = 3. * u.m / u.km
        dq1 = dq + 1. * u.mm / u.km
        assert dq1.value == 3.001
        assert dq1.unit == dq.unit

        dq2 = dq + 1.
        assert dq2.value == 1.003
        assert dq2.unit == u.dimensionless_unscaled

        # this test will check that operations with dimensionless Quantities
        # don't work
        with pytest.raises(u.UnitsError):
            self.q1 + u.Quantity(0.1, unit=u.Unit(""))

        with pytest.raises(u.UnitsError):
            self.q1 - u.Quantity(0.1, unit=u.Unit(""))

        # and test that scaling of integers works
        q = np.array([1, 2, 3]) * u.m / u.km
        q2 = q + np.array([4, 5, 6])
        assert q2.unit == u.dimensionless_unscaled
        assert_allclose(q2.value, np.array([4.001, 5.002, 6.003]))
        # but not if doing it inplace
        with pytest.raises(TypeError):
            q += np.array([1, 2, 3])
        # except if it is actually possible
        q = np.array([1, 2, 3]) * u.km / u.m
        q += np.array([4, 5, 6])
        assert q.unit == u.dimensionless_unscaled
        assert np.all(q.value == np.array([1004, 2005, 3006]))

    def test_complicated_operation(self):
        """ Perform a more complicated test """
        from .. import imperial

        # Multiple units
        distance = u.Quantity(15., u.meter)
        time = u.Quantity(11., u.second)

        velocity = (distance / time).to(imperial.mile / u.hour)
        assert_array_almost_equal(
            velocity.value, 3.05037, decimal=5)

        G = u.Quantity(6.673E-11, u.m ** 3 / u.kg / u.s ** 2)
        new_q = ((1. / (4. * np.pi * G)).to(u.pc ** -3 / u.s ** -2 * u.kg))

        # Area
        side1 = u.Quantity(11., u.centimeter)
        side2 = u.Quantity(7., u.centimeter)
        area = side1 * side2
        assert_array_almost_equal(area.value, 77., decimal=15)
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

    def test_numeric_converters(self):
        # float, int, long, and __index__ should only work for single
        # quantities, of appropriate type, and only if they are dimensionless.
        # (Check on __index__ is also a regression test for #1557)

        # quantities with units should never convert, or be usable as an index
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
        q2 = u.Quantity(1.23, u.m / u.km)

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
        q3 = u.Quantity(1.23, u.dimensionless_unscaled)

        assert float(q3) == 1.23
        assert int(q3) == 1
        if six.PY2:
            assert long(q3) == 1

        with pytest.raises(TypeError) as exc:
            q1 * ['a', 'b', 'c']
        assert exc.value.args[0] == index_err_msg

        # integer dimensionless unscaled is good for all
        q4 = u.Quantity(2, u.dimensionless_unscaled)

        assert float(q4) == 2.
        assert int(q4) == 2
        if six.PY2:
            assert long(q4) == 2

        assert q4 * ['a', 'b', 'c'] == ['a', 'b', 'c', 'a', 'b', 'c']

        # but arrays are not OK
        q5 = u.Quantity([1, 2], u.m)
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

    def test_array_converters(self):

        # Scalar quantity
        q = u.Quantity(1.23, u.m)
        assert np.all(np.array(q) == np.array([1.23]))

        # Array quantity
        q = u.Quantity([1., 2., 3.], u.m)
        assert np.all(np.array(q) == np.array([1., 2., 3.]))


def test_quantity_conversion():
    q1 = u.Quantity(0.1, unit=u.meter)

    new_quantity = q1.to(u.kilometer)
    assert new_quantity.value == 0.0001

    with pytest.raises(u.UnitsError):
        q1.to(u.zettastokes)


def test_quantity_conversion_with_equiv():
    q1 = u.Quantity(0.1, unit=u.meter)
    q2 = q1.to(u.Hz, equivalencies=u.spectral())
    assert_allclose(q2.value, 2997924580.0)

    q1 = u.Quantity(0.4, unit=u.arcsecond)
    q2 = q1.to(u.au, equivalencies=u.parallax())
    q3 = q2.to(u.arcminute, equivalencies=u.parallax())

    assert_allclose(q2.value, 515662.015)
    assert q2.unit == u.au
    assert_allclose(q3.value, 0.0066666667)
    assert q3.unit == u.arcminute


def test_quantity_conversion_equivalency_passed_on():
    class MySpectral(u.Quantity):
        _equivalencies = u.spectral()

        def __quantity_view__(self, obj, unit):
            return obj.view(MySpectral)

        def __quantity_instance__(self, *args, **kwargs):
            return MySpectral(*args, **kwargs)

    q1 = MySpectral([1000, 2000], unit=u.Hz)
    q2 = q1.to(u.nm)
    assert q2.unit == u.nm
    q3 = q2.to(u.Hz)
    assert q3.unit == u.Hz
    assert_allclose(q3.value, q1.value)
    q4 = MySpectral([1000, 2000], unit=u.nm)
    q5 = q4.to(u.Hz).to(u.nm)
    assert q5.unit == u.nm
    assert_allclose(q4.value, q5.value)


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

    def test_quantity_equality(self):
        assert u.Quantity(1000, unit='m') == u.Quantity(1, unit='km')
        assert not (u.Quantity(1, unit='m') == u.Quantity(1, unit='km'))
        with pytest.raises(u.UnitsError):
            u.Quantity(1, unit='m') == u.Quantity(1, unit='s')

    def test_quantity_comparison(self):
        assert u.Quantity(1100, unit=u.meter) > u.Quantity(1, unit=u.kilometer)
        assert u.Quantity(900, unit=u.meter) < u.Quantity(1, unit=u.kilometer)

        with pytest.raises(u.UnitsError):
            assert u.Quantity(1100, unit=u.meter) > u.Quantity(1, unit=u.second)

        with pytest.raises(u.UnitsError):
            assert u.Quantity(1100, unit=u.meter) < u.Quantity(1, unit=u.second)

        assert u.Quantity(1100, unit=u.meter) >= u.Quantity(1, unit=u.kilometer)
        assert u.Quantity(1000, unit=u.meter) >= u.Quantity(1, unit=u.kilometer)

        assert u.Quantity(900, unit=u.meter) <= u.Quantity(1, unit=u.kilometer)
        assert u.Quantity(1000, unit=u.meter) <= u.Quantity(1, unit=u.kilometer)

        with pytest.raises(u.UnitsError):
            assert u.Quantity(
                1100, unit=u.meter) >= u.Quantity(1, unit=u.second)

        with pytest.raises(u.UnitsError):
            assert u.Quantity(1100, unit=u.meter) <= u.Quantity(1, unit=u.second)

        assert u.Quantity(1200, unit=u.meter) != u.Quantity(1, unit=u.kilometer)

        with pytest.raises(u.UnitsError):
            assert u.Quantity(1100, unit=u.meter) != u.Quantity(1, unit=u.second)


class TestQuantityDisplay(object):
    scalarintq = u.Quantity(1, unit='m')
    scalarfloatq = u.Quantity(1.3, unit='m')
    arrq = u.Quantity([1, 2.3, 8.9], unit='m')

    def test_dimensionless_quantity_repr(self):
        q2 = u.Quantity(1, unit='m-1')
        assert repr(self.scalarintq * q2) == "<Quantity 1>"
        assert repr(self.arrq * q2) == "<Quantity [ 1. , 2.3, 8.9]>"

    def test_dimensionless_quantity_str(self):
        q2 = u.Quantity(1, unit='m-1')
        assert str(self.scalarintq * q2) == "1"
        assert str(self.arrq * q2) == "[ 1.   2.3  8.9]"

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


def test_arrays():
    """
    Test using quantites with array values
    """

    qsec = u.Quantity(np.arange(10), u.second)
    assert isinstance(qsec.value, np.ndarray)
    assert not qsec.isscalar

    # len and indexing should work for arrays
    assert len(qsec) == len(qsec.value)
    qsecsub25 = qsec[2:5]
    assert qsecsub25.unit == qsec.unit
    assert isinstance(qsecsub25, u.Quantity)
    assert len(qsecsub25) == 3

    # make sure isscalar, len, and indexing behave correcly for non-arrays.
    qsecnotarray = u.Quantity(10., u.second)
    assert qsecnotarray.isscalar
    with pytest.raises(TypeError):
        len(qsecnotarray)
    with pytest.raises(TypeError):
        qsecnotarray[0]

    qseclen0array = u.Quantity(np.array(10), u.second)
    # 0d numpy array should act basically like a scalar
    assert qseclen0array.isscalar
    with pytest.raises(TypeError):
        len(qseclen0array)
    with pytest.raises(TypeError):
        qseclen0array[0]
    assert isinstance(qseclen0array.value, int)

    # can also create from lists, will auto-convert to arrays
    qsec = u.Quantity(list(xrange(10)), u.second)
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

    # should create by unit multiplication, too
    qsec2 = np.arange(10) * u.second
    qsec3 = u.second * np.arange(10)

    assert np.all(qsec == qsec2)
    assert np.all(qsec2 == qsec3)

    # make sure numerical-converters fail when arrays are present
    with pytest.raises(TypeError):
        float(qsec)
    with pytest.raises(TypeError):
        int(qsec)
    if six.PY2:
        with pytest.raises(TypeError):
            long(qsec)


def test_array_indexing_slicing():
    q = np.array([1., 2., 3.]) * u.m
    assert q[0] == 1. * u.m
    assert np.all(q[0:2] == u.Quantity([1., 2.], u.m))


def test_array_setslice():
    q = np.array([1., 2., 3. ]) * u.m
    q[1:2] = np.array([400.]) * u.cm
    assert np.all(q == np.array([1., 4., 3.]) * u.m)


def test_inverse_quantity():
    """
    Regression test from issue #679
    """
    q = u.Quantity(4., u.meter / u.second)
    qot = q / 2
    toq = 2 / q
    npqot = q / np.array(2)

    assert npqot.value == 2.0
    assert npqot.unit == (u.meter / u.second)

    assert qot.value == 2.0
    assert qot.unit == (u.meter / u.second)

    assert toq.value == 0.5
    assert toq.unit == (u.second / u.meter)


def test_quantity_mutability():
    q = u.Quantity(9.8, u.meter / u.second / u.second)

    with pytest.raises(AttributeError):
        q.value = 3

    with pytest.raises(AttributeError):
        q.unit = u.kg


def test_quantity_initialized_with_quantity():
    q1 = u.Quantity(60, u.second)

    q2 = u.Quantity(q1, u.minute)
    assert q2.value == 1

    q3 = u.Quantity([q1, q2], u.second)
    assert q3[0].value == 60
    assert q3[1].value == 60

    q4 = u.Quantity([q2, q1])
    assert q4.unit == q2.unit
    assert q4[0].value == 1
    assert q4[1].value == 1


def test_quantity_string_unit():
    q1 = 1. * u.m / 's'
    assert q1.value == 1
    assert q1.unit == (u.m / u.s)

    q2 = q1 * "m"
    assert q2.unit == ((u.m * u.m) / u.s)


@raises(ValueError)
def test_quantity_invalid_unit_string():
    "foo" * u.m


def test_implicit_conversion():
    q = u.Quantity(1.0, u.meter)
    # Manually turn this on to simulate what might happen in a subclass
    q._include_easy_conversion_members = True
    assert_allclose(q.centimeter, 100)
    assert_allclose(q.cm, 100)
    assert_allclose(q.parsec, 3.240779289469756e-17)


def test_implicit_conversion_autocomplete():
    q = u.Quantity(1.0, u.meter)
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


def test_quantity_iterability():
    """Regressiont est for issue #878.

    Scalar quantities should not be iterable and should raise a type error on
    iteration.
    """

    q1 = [15.0, 17.0] * u.m
    assert isiterable(q1)

    q2 = six.next(iter(q1))
    assert q2 == 15.0 * u.m
    assert not isiterable(q2)
    pytest.raises(TypeError, iter, q2)


def test_copy():

    q1 = u.Quantity(np.array([1., 2., 3.]), unit=u.m)
    q2 = q1.copy()

    assert np.all(q1.value == q2.value)
    assert q1.unit == q2.unit
    assert q1.dtype == q2.dtype

    assert q1.value is not q2.value


def test_deepcopy():
    q1 = u.Quantity(np.array([1., 2., 3.]), unit=u.m)
    q2 = copy.deepcopy(q1)

    assert isinstance(q2, u.Quantity)
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


def test_quantity_pickelability():
    """
    Testing pickleability of quantity
    """

    q1 = np.arange(10) * u.m

    q2 = pickle.loads(pickle.dumps(q1))

    assert np.all(q1.value == q2.value)
    assert q1.unit.is_equivalent(q2.unit)
    assert q1.unit == q2.unit


def test_quantity_from_string():
    with pytest.raises(TypeError):
        q = u.Quantity(u.m * "5")
        # the reverse should also fail once #1408 is in

    with pytest.raises(TypeError):
        q = u.Quantity('5', u.m)

    with pytest.raises(TypeError):
        q = u.Quantity(['5'], u.m)

    with pytest.raises(TypeError):
        q = u.Quantity(np.array(['5']), u.m)


def test_unsupported():
    q1 = np.arange(10) * u.m

    with pytest.raises(TypeError):
        q2 = np.bitwise_and(q1, q1)


def test_unit_identity():
    q = 1.0 * u.hour
    assert q.unit is u.hour


def test_quantity_to_view():
    q1 = np.array([1000, 2000]) * u.m
    q2 = q1.to(u.km)
    assert q1.value[0] == 1000
    assert q2.value[0] == 1


@raises(ValueError)
def test_quantity_tuple_power():
    (5.0 * u.m) ** (1, 2)

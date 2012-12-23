# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
    Test the Quantity class and related.
"""

from __future__ import absolute_import, unicode_literals, division, print_function

import pytest
import numpy as np

from ...tests.helper import raises
from ...tests.compat import assert_allclose
from ... import units as u

""" The Quantity class will represent a number + unit + uncertainty """

class TestQuantityCreation():

    def test_1(self):
        # create objects through operations with Unit objects:

        quantity = 11.42 * u.meter # returns a Quantity object
        assert isinstance(quantity,u.Quantity)
        quantity = u.meter * 11.42 # returns a Quantity object
        assert isinstance(quantity,u.Quantity)

        quantity = 11.42 / u.meter
        assert isinstance(quantity,u.Quantity)
        quantity = u.meter / 11.42
        assert isinstance(quantity,u.Quantity)

        quantity = 11.42 * u.meter / u.second
        assert isinstance(quantity,u.Quantity)

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

        with pytest.raises(TypeError):
            q3 = u.Quantity(11.412)

    def test_3(self):
        #with pytest.raises(u.UnitsException):
        with pytest.raises(ValueError): # Until @mdboom fixes the errors in units
            q1 = u.Quantity(11.412, unit="testingggg")

    def test_unit_property(self):
        # test getting and setting 'unit' attribute
        q1 = u.Quantity(11.4, unit=u.meter)

        with pytest.raises(AttributeError):
            q1.unit = u.centimeter

class TestQuantityOperations():
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

        new_q = u.Quantity(1500.1,u.m) + u.Quantity(13.5,u.km)
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
        assert new_quantity.unit == (u.meter*u.centimeter)

        # Take units from left object, q2
        new_quantity = self.q2 * self.q1
        assert new_quantity.value == 91.36
        assert new_quantity.unit == (u.centimeter*u.meter)

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
        np.testing.assert_array_almost_equal(new_quantity.value, 1.4275, decimal=5)
        assert new_quantity.unit == (u.meter / u.centimeter)

        # Take units from left object, q2
        new_quantity = self.q2 / self.q1
        np.testing.assert_array_almost_equal(new_quantity.value, 0.70052539404553416, decimal=16)
        assert new_quantity.unit == (u.centimeter / u.meter)

        q1 = u.Quantity(11.4, unit=u.meter)
        q2 = u.Quantity(10.0, unit=u.second)
        new_quantity = q1 / q2
        np.testing.assert_array_almost_equal(new_quantity.value, 1.14, decimal=10)
        assert new_quantity.unit == (u.meter / u.second)

        # divide with a number
        new_quantity = self.q1 / 10.
        assert new_quantity.value == 1.142
        assert new_quantity.unit == u.meter

        # divide with a number
        new_quantity = 11.42 / self.q1
        assert new_quantity.value == 1.
        assert new_quantity.unit == u.Unit("1/m")

    def test_power(self):
        # raise quantity to a power
        new_quantity = self.q1**2
        np.testing.assert_array_almost_equal(new_quantity.value, 130.4164, decimal=5)
        assert new_quantity.unit == u.Unit("m^2")

        new_quantity = self.q1**3
        np.testing.assert_array_almost_equal(new_quantity.value, 1489.355288, decimal=7)
        assert new_quantity.unit == u.Unit("m^3")

    def test_incompatible_units(self):
        """ When trying to add or subtract units that aren't compatible, throw an error """

        q1 = u.Quantity(11.412, unit=u.meter)
        q2 = u.Quantity(21.52, unit=u.second)

        with pytest.raises(u.UnitsException):
            new_q = q1 + q2

    def test_dimensionless_operations(self):
        # this test will check that operations with dimensionless Quantities don't work

        with pytest.raises(u.UnitsException):
            self.q1 + u.Quantity(0.1, unit=u.Unit(""))

        with pytest.raises(u.UnitsException):
            self.q1 - u.Quantity(0.1, unit=u.Unit(""))

    def test_complicated_operation(self):
        """ Perform a more complicated test """

        # Multiple units
        distance = u.Quantity(15., u.meter)
        time = u.Quantity(11., u.second)

        velocity = (distance / time).to(u.mile/u.hour)
        np.testing.assert_array_almost_equal(velocity.value, 3.05037, decimal=5)

        G = u.Quantity(6.673E-11, u.m**3/u.kg/u.s**2)
        new_q = ((1. / (4.*np.pi*G)).to(u.pc**-3/u.s**-2*u.kg))

        # Area
        side1 = u.Quantity(11., u.centimeter)
        side2 = u.Quantity(7., u.centimeter)
        area = side1 * side2
        np.testing.assert_array_almost_equal(area.value, 77., decimal=15)
        assert area.unit == u.cm*u.cm

    def test_comparison(self):
        # equality/ non-equality is straightforward for quantity objects
        assert (1 / (u.cm * u.cm)) == 1 * u.cm ** -2
        assert 1 * u.m == 100 * u.cm
        assert 1 * u.m != 1 * u.cm

        #here one is a unit, which is an invalid comparison
        assert 1. * u.cm * u.cm * u.cm != u.cm ** 3

        #mismatched types should never work
        assert not 1. * u.cm == 1.
        assert 1. * u.cm != 1.

    def test_numeric_converters(self):
        q = u.Quantity(1.23, u.m)

        assert float(q) == 1.23
        assert int(q) == 1
        assert long(q) == 1L


def test_quantity_conversion():
    q1 = u.Quantity(0.1, unit=u.meter)

    new_quantity = q1.to(u.kilometer)
    assert new_quantity.value == 0.0001

    with pytest.raises(u.UnitsException):
        q1.to(u.zettastokes)

def test_si():
    q1 = 10. * u.m * u.s**2 / (200. * u.ms)**2 # 250 meters
    assert q1.si.value == 250
    assert q1.si.unit == u.m

    q = 10.*u.m # 10 meters
    assert q.si.value == 10
    assert q.si.unit == u.m

    q = 10./u.m # 10 1 / meters
    assert q.si.value == 10
    assert q.si.unit == (1/u.m)

def test_cgs():
    q1 = 10. * u.cm * u.s**2 / (200. * u.ms)**2 # 250 centimeters
    assert q1.cgs.value == 250
    assert q1.cgs.unit == u.cm

    q = 10.*u.m # 10 centimeters
    assert q.cgs.value == 1000
    assert q.cgs.unit == u.cm

    q = 10./u.cm # 10 1 / centimeters
    assert q.cgs.value == 10
    assert q.cgs.unit == (1/u.cm)

    q = 10.*u.Pa # 10 pascals
    assert q.cgs.value == 100
    assert q.cgs.unit == u.barye

class TestQuantityComparison():
    def test_quantity_equality(self):
        assert u.Quantity(1000, unit='m') == u.Quantity(1, unit='km')
        assert not (u.Quantity(1, unit='m') == u.Quantity(1, unit='km'))
        with pytest.raises(u.UnitsException):
            u.Quantity(1, unit='m') == u.Quantity(1, unit='s')

    def test_quantity_comparison(self):
        assert u.Quantity(1100, unit=u.meter) > u.Quantity(1, unit=u.kilometer)
        assert u.Quantity(900, unit=u.meter) < u.Quantity(1, unit=u.kilometer)

        with pytest.raises(u.UnitsException):
            assert u.Quantity(1100, unit=u.meter) > u.Quantity(1, unit=u.second)

        with pytest.raises(u.UnitsException):
            assert u.Quantity(1100, unit=u.meter) < u.Quantity(1, unit=u.second)

        assert u.Quantity(1100, unit=u.meter) >= u.Quantity(1, unit=u.kilometer)
        assert u.Quantity(1000, unit=u.meter) >= u.Quantity(1, unit=u.kilometer)

        assert u.Quantity(900, unit=u.meter) <= u.Quantity(1, unit=u.kilometer)
        assert u.Quantity(1000, unit=u.meter) <= u.Quantity(1, unit=u.kilometer)

        with pytest.raises(u.UnitsException):
            assert u.Quantity(1100, unit=u.meter) >= u.Quantity(1, unit=u.second)

        with pytest.raises(u.UnitsException):
            assert u.Quantity(1100, unit=u.meter) <= u.Quantity(1, unit=u.second)

        assert u.Quantity(1200, unit=u.meter) != u.Quantity(1, unit=u.kilometer)

        with pytest.raises(u.UnitsException):
            assert u.Quantity(1100, unit=u.meter) != u.Quantity(1, unit=u.second)

class TestQuantityDisplay():

    def test_quantity_str(self):
        q1 = u.Quantity(1, unit='m')
        assert str(q1) == "1 m"

    def test_quantity_repr(self):
        q1 = u.Quantity(1, unit='m')
        assert repr(q1) == "<Quantity 1 m>"


def test_decompose():
    q1 = 5 * u.N
    assert q1.decomposed_unit == (5 * u.kg * u.m * u.s ** -2)

def test_arrays():
    """
    Test using quantites with array values
    """
    from numpy.testing import assert_array_equal

    qsec = u.Quantity(np.arange(10), u.second)
    assert isinstance(qsec.value, np.ndarray)
    assert not qsec.isscalar

    #len and indexing should work for arrays
    assert len(qsec) == len(qsec.value)
    qsecsub25 = qsec[2:5]
    assert qsecsub25.unit == qsec.unit
    assert isinstance(qsecsub25, u.Quantity)
    assert len(qsecsub25) == 3

    #make sure isscalar, len, and indexing behave correcly for non-arrays.
    qsecnotarray = u.Quantity(10., u.second)
    assert qsecnotarray.isscalar
    with pytest.raises(TypeError):
        len(qsecnotarray)
    with pytest.raises(TypeError):
        qsecnotarray[0]

    qseclen0array = u.Quantity(np.array(10), u.second)
    # 0d numpy array should act basically like a scalar, but still keep its identity as a numpy array
    assert qseclen0array.isscalar
    with pytest.raises(TypeError):
        len(qseclen0array)
    with pytest.raises(TypeError):
        qseclen0array[0]
    assert isinstance(qseclen0array.value, np.ndarray)


    #can also create from lists, will auto-convert to arrays
    qsec = u.Quantity(range(10), u.second)
    assert isinstance(qsec.value, np.ndarray)

    #quantity math should work with arrays
    assert_array_equal((qsec * 2).value, (np.arange(10) * 2))
    assert_array_equal((qsec / 2).value, (np.arange(10) / 2))
    #quantity addition/subtraction should *not* work with arrays b/c unit ambiguous
    with pytest.raises(TypeError):
        assert_array_equal((qsec + 2).value, (np.arange(10) + 2))
    with pytest.raises(TypeError):
        assert_array_equal((qsec - 2).value, (np.arange(10) + 2))

    #should create by unit multiplication, too
    qsec2 = np.arange(10) * u.second
    qsec3 = u.second * np.arange(10)

    assert np.all(qsec == qsec2)
    assert np.all(qsec2 == qsec3)

    # make sure numerical-converters fail when arrays are present
    with pytest.raises(TypeError):
        float(qsec)
    with pytest.raises(TypeError):
        int(qsec)
    with pytest.raises(TypeError):
        long(qsec)

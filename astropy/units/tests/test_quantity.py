# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
    Test the Quantity class and related.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import sys
import copy
import decimal
from fractions import Fraction

import pytest
import numpy as np
from numpy.testing import (assert_allclose, assert_array_equal,
                           assert_array_almost_equal)

from ...tests.helper import raises
from ...utils import isiterable, minversion
from ...utils.compat import NUMPY_LT_1_10, NUMPY_LT_1_14
from ...utils.compat.numpy import matmul
from ... import units as u
from ...units.quantity import _UNIT_NOT_INITIALISED
from ...extern.six.moves import range
from ...extern.six.moves import cPickle as pickle
from ...extern import six

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from distutils.version import LooseVersion
    MATPLOTLIB_LT_15 = LooseVersion(matplotlib.__version__) < LooseVersion("1.5")
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


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

    def test_nan_inf(self):
        # Not-a-number
        q = u.Quantity('nan', unit='cm')
        assert np.isnan(q.value)

        q = u.Quantity('NaN', unit='cm')
        assert np.isnan(q.value)

        q = u.Quantity('-nan', unit='cm')  # float() allows this
        assert np.isnan(q.value)

        q = u.Quantity('nan cm')
        assert np.isnan(q.value)
        assert q.unit == u.cm

        # Infinity
        q = u.Quantity('inf', unit='cm')
        assert np.isinf(q.value)

        q = u.Quantity('-inf', unit='cm')
        assert np.isinf(q.value)

        q = u.Quantity('inf cm')
        assert np.isinf(q.value)
        assert q.unit == u.cm

        q = u.Quantity('Infinity', unit='cm')  # float() allows this
        assert np.isinf(q.value)

        # make sure these strings don't parse...
        with pytest.raises(TypeError):
            q = u.Quantity('', unit='cm')

        with pytest.raises(TypeError):
            q = u.Quantity('spam', unit='cm')

    def test_unit_property(self):
        # test getting and setting 'unit' attribute
        q1 = u.Quantity(11.4, unit=u.meter)

        with pytest.raises(AttributeError):
            q1.unit = u.cm

    def test_preserve_dtype(self):
        """Test that if an explicit dtype is given, it is used, while if not,
        numbers are converted to float (including decimal.Decimal, which
        numpy converts to an object; closes #1419)
        """
        # If dtype is specified, use it, but if not, convert int, bool to float
        q1 = u.Quantity(12, unit=u.m / u.s, dtype=int)
        assert q1.dtype == int

        q2 = u.Quantity(q1)
        assert q2.dtype == float
        assert q2.value == float(q1.value)
        assert q2.unit == q1.unit

        # but we should preserve float32
        a3 = np.array([1., 2.], dtype=np.float32)
        q3 = u.Quantity(a3, u.yr)
        assert q3.dtype == a3.dtype
        # items stored as objects by numpy should be converted to float
        # by default
        q4 = u.Quantity(decimal.Decimal('10.25'), u.m)
        assert q4.dtype == float

        q5 = u.Quantity(decimal.Decimal('10.25'), u.m, dtype=object)
        assert q5.dtype == object

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

    def test_subok(self):
        """Test subok can be used to keep class, or to insist on Quantity"""
        class MyQuantitySubclass(u.Quantity):
            pass

        myq = MyQuantitySubclass(np.arange(10.), u.m)
        # try both with and without changing the unit
        assert type(u.Quantity(myq)) is u.Quantity
        assert type(u.Quantity(myq, subok=True)) is MyQuantitySubclass
        assert type(u.Quantity(myq, u.km)) is u.Quantity
        assert type(u.Quantity(myq, u.km, subok=True)) is MyQuantitySubclass

    def test_order(self):
        """Test that order is correctly propagated to np.array"""
        ac = np.array(np.arange(10.), order='C')
        qcc = u.Quantity(ac, u.m, order='C')
        assert qcc.flags['C_CONTIGUOUS']
        qcf = u.Quantity(ac, u.m, order='F')
        assert qcf.flags['F_CONTIGUOUS']
        qca = u.Quantity(ac, u.m, order='A')
        assert qca.flags['C_CONTIGUOUS']
        # check it works also when passing in a quantity
        assert u.Quantity(qcc, order='C').flags['C_CONTIGUOUS']
        assert u.Quantity(qcc, order='A').flags['C_CONTIGUOUS']
        assert u.Quantity(qcc, order='F').flags['F_CONTIGUOUS']

        af = np.array(np.arange(10.), order='F')
        qfc = u.Quantity(af, u.m, order='C')
        assert qfc.flags['C_CONTIGUOUS']
        qff = u.Quantity(ac, u.m, order='F')
        assert qff.flags['F_CONTIGUOUS']
        qfa = u.Quantity(af, u.m, order='A')
        assert qfa.flags['F_CONTIGUOUS']
        assert u.Quantity(qff, order='C').flags['C_CONTIGUOUS']
        assert u.Quantity(qff, order='A').flags['F_CONTIGUOUS']
        assert u.Quantity(qff, order='F').flags['F_CONTIGUOUS']

    def test_ndmin(self):
        """Test that ndmin is correctly propagated to np.array"""
        a = np.arange(10.)
        q1 = u.Quantity(a, u.m, ndmin=1)
        assert q1.ndim == 1 and q1.shape == (10,)
        q2 = u.Quantity(a, u.m, ndmin=2)
        assert q2.ndim == 2 and q2.shape == (1, 10)
        # check it works also when passing in a quantity
        q3 = u.Quantity(q1, u.m, ndmin=3)
        assert q3.ndim == 3 and q3.shape == (1, 1, 10)

    def test_non_quantity_with_unit(self):
        """Test that unit attributes in objects get recognized."""
        class MyQuantityLookalike(np.ndarray):
            pass

        a = np.arange(3.)
        mylookalike = a.copy().view(MyQuantityLookalike)
        mylookalike.unit = 'm'
        q1 = u.Quantity(mylookalike)
        assert isinstance(q1, u.Quantity)
        assert q1.unit is u.m
        assert np.all(q1.value == a)

        q2 = u.Quantity(mylookalike, u.mm)
        assert q2.unit is u.mm
        assert np.all(q2.value == 1000.*a)

        q3 = u.Quantity(mylookalike, copy=False)
        assert np.all(q3.value == mylookalike)
        q3[2] = 0
        assert q3[2] == 0.
        assert mylookalike[2] == 0.

        mylookalike = a.copy().view(MyQuantityLookalike)
        mylookalike.unit = u.m
        q4 = u.Quantity(mylookalike, u.mm, copy=False)
        q4[2] = 0
        assert q4[2] == 0.
        assert mylookalike[2] == 2.

        mylookalike.unit = 'nonsense'
        with pytest.raises(TypeError):
            u.Quantity(mylookalike)


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

    @pytest.mark.skipif(sys.version_info[:2] < (3, 5),
                        reason="__matmul__ only introduced in Python 3.5")
    def test_matrix_multiplication(self):
        # We cannot write @ as an operator since class gets parsed also on
        # Python <3.5, so we use eval instead.
        a = np.eye(3)
        q = a * u.m
        result1 = eval("q @ a")
        assert np.all(result1 == q)
        result2 = eval("a @ q")
        assert np.all(result2 == q)
        result3 = eval("q @ q")
        assert np.all(result3 == a * u.m ** 2)
        # less trivial case.
        q2 = np.array([[[1., 0., 0.],
                        [0., 1., 0.],
                        [0., 0., 1.]],
                       [[0., 1., 0.],
                        [0., 0., 1.],
                        [1., 0., 0.]],
                       [[0., 0., 1.],
                        [1., 0., 0.],
                        [0., 1., 0.]]]) / u.s
        result4 = eval("q @ q2")
        assert np.all(result4 == matmul(a, q2.value) * q.unit * q2.unit)

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

    def test_non_number_type(self):
        q1 = u.Quantity(11.412, unit=u.meter)
        type_err_msg = ("Unsupported operand type(s) for ufunc add: "
                        "'Quantity' and 'dict'")
        with pytest.raises(TypeError) as exc:
            q1 + {'a': 1}
        assert exc.value.args[0] == type_err_msg

        with pytest.raises(TypeError):
            q1 + u.meter

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
        q = u.Quantity(np.array([1, 2, 3]), u.m / u.km, dtype=int)
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

        # when one is a unit, Quantity does not know what to do,
        # but unit is fine with it, so it still works
        unit = u.cm**3
        q = 1. * unit
        assert q.__eq__(unit) is NotImplemented
        assert unit.__eq__(q) is True
        assert q == unit
        q = 1000. * u.mm**3
        assert q == unit

        # mismatched types should never work
        assert not 1. * u.cm == 1.
        assert 1. * u.cm != 1.

    def test_numeric_converters(self):
        # float, int, long, and __index__ should only work for single
        # quantities, of appropriate type, and only if they are dimensionless.
        # for index, this should be unscaled as well
        # (Check on __index__ is also a regression test for #1557)

        # quantities with units should never convert, or be usable as an index
        q1 = u.Quantity(1, u.m)

        converter_err_msg = ("only dimensionless scalar quantities "
                             "can be converted to Python scalars")
        index_err_msg = ("only integer dimensionless scalar quantities "
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

        # We used to test `q1 * ['a', 'b', 'c'] here, but that that worked
        # at all was a really odd confluence of bugs.  Since it doesn't work
        # in numpy >=1.10 any more, just go directly for `__index__` (which
        # makes the test more similar to the `int`, `long`, etc., tests).
        with pytest.raises(TypeError) as exc:
            q1.__index__()
        assert exc.value.args[0] == index_err_msg

        # dimensionless but scaled is OK, however
        q2 = u.Quantity(1.23, u.m / u.km)

        assert float(q2) == float(q2.to_value(u.dimensionless_unscaled))
        assert int(q2) == int(q2.to_value(u.dimensionless_unscaled))

        if six.PY2:
            assert long(q2) == long(q2.to_value(u.dimensionless_unscaled))

        with pytest.raises(TypeError) as exc:
            q2.__index__()
        assert exc.value.args[0] == index_err_msg

        # dimensionless unscaled is OK, though for index needs to be int
        q3 = u.Quantity(1.23, u.dimensionless_unscaled)

        assert float(q3) == 1.23
        assert int(q3) == 1
        if six.PY2:
            assert long(q3) == 1

        with pytest.raises(TypeError) as exc:
            q3.__index__()
        assert exc.value.args[0] == index_err_msg

        # integer dimensionless unscaled is good for all
        q4 = u.Quantity(2, u.dimensionless_unscaled, dtype=int)

        assert float(q4) == 2.
        assert int(q4) == 2
        if six.PY2:
            assert long(q4) == 2

        assert q4.__index__() == 2

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
            q5.__index__()
        assert exc.value.args[0] == index_err_msg

    # Fails for numpy >=1.10; see https://github.com/numpy/numpy/issues/5074
    # It seems unlikely this will be resolved, so xfail'ing it.
    @pytest.mark.xfail(not NUMPY_LT_1_10,
                       reason="list multiplication only works for numpy <=1.10")
    def test_numeric_converter_to_index_in_practice(self):
        """Test that use of __index__ actually works."""
        q4 = u.Quantity(2, u.dimensionless_unscaled, dtype=int)
        assert q4 * ['a', 'b', 'c'] == ['a', 'b', 'c', 'a', 'b', 'c']

    def test_array_converters(self):

        # Scalar quantity
        q = u.Quantity(1.23, u.m)
        assert np.all(np.array(q) == np.array([1.23]))

        # Array quantity
        q = u.Quantity([1., 2., 3.], u.m)
        assert np.all(np.array(q) == np.array([1., 2., 3.]))


def test_quantity_conversion():
    q1 = u.Quantity(0.1, unit=u.meter)
    value = q1.value
    assert value == 0.1
    value_in_km = q1.to_value(u.kilometer)
    assert value_in_km == 0.0001
    new_quantity = q1.to(u.kilometer)
    assert new_quantity.value == 0.0001

    with pytest.raises(u.UnitsError):
        q1.to(u.zettastokes)
    with pytest.raises(u.UnitsError):
        q1.to_value(u.zettastokes)


def test_quantity_value_views():
    q1 = u.Quantity([1., 2.], unit=u.meter)
    # views if the unit is the same.
    v1 = q1.value
    v1[0] = 0.
    assert np.all(q1 == [0., 2.] * u.meter)
    v2 = q1.to_value()
    v2[1] = 3.
    assert np.all(q1 == [0., 3.] * u.meter)
    v3 = q1.to_value('m')
    v3[0] = 1.
    assert np.all(q1 == [1., 3.] * u.meter)
    v4 = q1.to_value('cm')
    v4[0] = 0.
    # copy if different unit.
    assert np.all(q1 == [1., 3.] * u.meter)


def test_quantity_conversion_with_equiv():
    q1 = u.Quantity(0.1, unit=u.meter)
    v2 = q1.to_value(u.Hz, equivalencies=u.spectral())
    assert_allclose(v2, 2997924580.0)
    q2 = q1.to(u.Hz, equivalencies=u.spectral())
    assert_allclose(q2.value, v2)

    q1 = u.Quantity(0.4, unit=u.arcsecond)
    v2 = q1.to_value(u.au, equivalencies=u.parallax())
    q2 = q1.to(u.au, equivalencies=u.parallax())
    v3 = q2.to_value(u.arcminute, equivalencies=u.parallax())
    q3 = q2.to(u.arcminute, equivalencies=u.parallax())

    assert_allclose(v2, 515662.015)
    assert_allclose(q2.value, v2)
    assert q2.unit == u.au
    assert_allclose(v3, 0.0066666667)
    assert_allclose(q3.value, v3)
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

# Regression test for issue #2315, divide-by-zero error when examining 0*unit


def test_self_equivalency():
    assert u.deg.is_equivalent(0*u.radian)
    assert u.deg.is_equivalent(1*u.radian)


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
        # for ==, !=, return False, True if units do not match
        assert (u.Quantity(1100, unit=u.m) != u.Quantity(1, unit=u.s)) is True
        assert (u.Quantity(1100, unit=u.m) == u.Quantity(1, unit=u.s)) is False

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


class TestQuantityDisplay(object):
    scalarintq = u.Quantity(1, unit='m', dtype=int)
    scalarfloatq = u.Quantity(1.3, unit='m')
    arrq = u.Quantity([1, 2.3, 8.9], unit='m')

    def test_dimensionless_quantity_repr(self):
        q2 = u.Quantity(1., unit='m-1')
        q3 = u.Quantity(1, unit='m-1', dtype=int)
        if NUMPY_LT_1_14:
            assert repr(self.scalarintq * q2) == "<Quantity 1.0>"
            assert repr(self.arrq * q2) == "<Quantity [ 1. , 2.3, 8.9]>"
        else:
            assert repr(self.scalarintq * q2) == "<Quantity 1.>"
            assert repr(self.arrq * q2) == "<Quantity [1. , 2.3, 8.9]>"
        assert repr(self.scalarintq * q3) == "<Quantity 1>"

    def test_dimensionless_quantity_str(self):
        q2 = u.Quantity(1., unit='m-1')
        q3 = u.Quantity(1, unit='m-1', dtype=int)
        assert str(self.scalarintq * q2) == "1.0"
        assert str(self.scalarintq * q3) == "1"
        if NUMPY_LT_1_14:
            assert str(self.arrq * q2) == "[ 1.   2.3  8.9]"
        else:
            assert str(self.arrq * q2) == "[1.  2.3 8.9]"

    def test_dimensionless_quantity_format(self):
        q1 = u.Quantity(3.14)
        assert format(q1, '.2f') == '3.14'

    def test_scalar_quantity_str(self):
        assert str(self.scalarintq) == "1 m"
        assert str(self.scalarfloatq) == "1.3 m"

    def test_scalar_quantity_repr(self):
        assert repr(self.scalarintq) == "<Quantity 1 m>"
        assert repr(self.scalarfloatq) == "<Quantity 1.3 m>"

    def test_array_quantity_str(self):
        if NUMPY_LT_1_14:
            assert str(self.arrq) == "[ 1.   2.3  8.9] m"
        else:
            assert str(self.arrq) == "[1.  2.3 8.9] m"

    def test_array_quantity_repr(self):
        if NUMPY_LT_1_14:
            assert repr(self.arrq) == "<Quantity [ 1. , 2.3, 8.9] m>"
        else:
            assert repr(self.arrq) == "<Quantity [1. , 2.3, 8.9] m>"

    def test_scalar_quantity_format(self):
        assert format(self.scalarintq, '02d') == "01 m"
        assert format(self.scalarfloatq, '.1f') == "1.3 m"
        assert format(self.scalarfloatq, '.0f') == "1 m"

    def test_uninitialized_unit_format(self):
        bad_quantity = np.arange(10.).view(u.Quantity)
        assert str(bad_quantity).endswith(_UNIT_NOT_INITIALISED)
        assert repr(bad_quantity).endswith(_UNIT_NOT_INITIALISED + '>')

    def test_repr_latex(self):
        from ...units.quantity import conf

        q2scalar = u.Quantity(1.5e14, 'm/s')
        assert self.scalarintq._repr_latex_() == r'$1 \; \mathrm{m}$'
        assert self.scalarfloatq._repr_latex_() == r'$1.3 \; \mathrm{m}$'
        assert (q2scalar._repr_latex_() ==
                r'$1.5 \times 10^{14} \; \mathrm{\frac{m}{s}}$')
        assert self.arrq._repr_latex_() == r'$[1,~2.3,~8.9] \; \mathrm{m}$'

        qmed = np.arange(100)*u.m
        qbig = np.arange(1000)*u.m
        qvbig = np.arange(10000)*1e9*u.m

        pops = np.get_printoptions()
        oldlat = conf.latex_array_threshold
        try:
            # check precision behavior
            q = u.Quantity(987654321.123456789, 'm/s')
            qa = np.array([7.89123, 123456789.987654321, 0]) * u.cm
            np.set_printoptions(precision=8)
            assert q._repr_latex_() == r'$9.8765432 \times 10^{8} \; \mathrm{\frac{m}{s}}$'
            assert qa._repr_latex_() == r'$[7.89123,~1.2345679 \times 10^{8},~0] \; \mathrm{cm}$'
            np.set_printoptions(precision=2)
            assert q._repr_latex_() == r'$9.9 \times 10^{8} \; \mathrm{\frac{m}{s}}$'
            assert qa._repr_latex_() == r'$[7.9,~1.2 \times 10^{8},~0] \; \mathrm{cm}$'

            # check thresholding behavior
            conf.latex_array_threshold = 100  # should be default
            lsmed = qmed._repr_latex_()
            assert r'\dots' not in lsmed
            lsbig = qbig._repr_latex_()
            assert r'\dots' in lsbig
            lsvbig = qvbig._repr_latex_()
            assert r'\dots' in lsvbig

            conf.latex_array_threshold = 1001
            lsmed = qmed._repr_latex_()
            assert r'\dots' not in lsmed
            lsbig = qbig._repr_latex_()
            assert r'\dots' not in lsbig
            lsvbig = qvbig._repr_latex_()
            assert r'\dots' in lsvbig

            conf.latex_array_threshold = -1  # means use the numpy threshold
            np.set_printoptions(threshold=99)
            lsmed = qmed._repr_latex_()
            assert r'\dots' in lsmed
            lsbig = qbig._repr_latex_()
            assert r'\dots' in lsbig
            lsvbig = qvbig._repr_latex_()
            assert r'\dots' in lsvbig
        finally:
            # prevent side-effects from influencing other tests
            np.set_printoptions(**pops)
            conf.latex_array_threshold = oldlat

        qinfnan = [np.inf, -np.inf, np.nan] * u.m
        assert qinfnan._repr_latex_() == r'$[\infty,~-\infty,~{\rm NaN}] \; \mathrm{m}$'


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

    qseclen0array = u.Quantity(np.array(10), u.second, dtype=int)
    # 0d numpy array should act basically like a scalar
    assert qseclen0array.isscalar
    with pytest.raises(TypeError):
        len(qseclen0array)
    with pytest.raises(TypeError):
        qseclen0array[0]
    assert isinstance(qseclen0array.value, int)

    # but with multiple dtypes, single elements are OK; need to use str()
    # since numpy under python2 cannot handle unicode literals
    a = np.array([(1., 2., 3.), (4., 5., 6.), (7., 8., 9.)],
                 dtype=[(str('x'), np.float),
                        (str('y'), np.float),
                        (str('z'), np.float)])
    qkpc = u.Quantity(a, u.kpc)
    assert not qkpc.isscalar
    qkpc0 = qkpc[0]
    assert qkpc0.value == a[0].item()
    assert qkpc0.unit == qkpc.unit
    assert isinstance(qkpc0, u.Quantity)
    assert qkpc0.isscalar
    qkpcx = qkpc['x']
    assert np.all(qkpcx.value == a['x'])
    assert qkpcx.unit == qkpc.unit
    assert isinstance(qkpcx, u.Quantity)
    assert not qkpcx.isscalar
    qkpcx1 = qkpc['x'][1]
    assert qkpcx1.unit == qkpc.unit
    assert isinstance(qkpcx1, u.Quantity)
    assert qkpcx1.isscalar
    qkpc1x = qkpc[1]['x']
    assert qkpc1x.isscalar
    assert qkpc1x == qkpcx1

    # can also create from lists, will auto-convert to arrays
    qsec = u.Quantity(list(range(10)), u.second)
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
    q = np.array([1., 2., 3.]) * u.m
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

    q1 = u.Quantity(np.array([[1., 2., 3.], [4., 5., 6.]]), unit=u.m)
    q2 = q1.copy()

    assert np.all(q1.value == q2.value)
    assert q1.unit == q2.unit
    assert q1.dtype == q2.dtype
    assert q1.value is not q2.value

    q3 = q1.copy(order='F')
    assert q3.flags['F_CONTIGUOUS']
    assert np.all(q1.value == q3.value)
    assert q1.unit == q3.unit
    assert q1.dtype == q3.dtype
    assert q1.value is not q3.value

    q4 = q1.copy(order='C')
    assert q4.flags['C_CONTIGUOUS']
    assert np.all(q1.value == q4.value)
    assert q1.unit == q4.unit
    assert q1.dtype == q4.dtype
    assert q1.value is not q4.value


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


def test_quantity_initialisation_from_string():
    q = u.Quantity('1')
    assert q.unit == u.dimensionless_unscaled
    assert q.value == 1.
    q = u.Quantity('1.5 m/s')
    assert q.unit == u.m/u.s
    assert q.value == 1.5
    assert u.Unit(q) == u.Unit('1.5 m/s')
    q = u.Quantity('.5 m')
    assert q == u.Quantity(0.5, u.m)
    q = u.Quantity('-1e1km')
    assert q == u.Quantity(-10, u.km)
    q = u.Quantity('-1e+1km')
    assert q == u.Quantity(-10, u.km)
    q = u.Quantity('+.5km')
    assert q == u.Quantity(.5, u.km)
    q = u.Quantity('+5e-1km')
    assert q == u.Quantity(.5, u.km)
    q = u.Quantity('5', u.m)
    assert q == u.Quantity(5., u.m)
    q = u.Quantity('5 km', u.m)
    assert q.value == 5000.
    assert q.unit == u.m
    q = u.Quantity('5Em')
    assert q == u.Quantity(5., u.Em)

    with pytest.raises(TypeError):
        u.Quantity('')
    with pytest.raises(TypeError):
        u.Quantity('m')
    with pytest.raises(TypeError):
        u.Quantity('1.2.3 deg')
    with pytest.raises(TypeError):
        u.Quantity('1+deg')
    with pytest.raises(TypeError):
        u.Quantity('1-2deg')
    with pytest.raises(TypeError):
        u.Quantity('1.2e-13.3m')
    with pytest.raises(TypeError):
        u.Quantity(['5'])
    with pytest.raises(TypeError):
        u.Quantity(np.array(['5']))
    with pytest.raises(ValueError):
        u.Quantity('5E')
    with pytest.raises(ValueError):
        u.Quantity('5 foo')


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


def test_quantity_fraction_power():
    q = (25.0 * u.m**2) ** Fraction(1, 2)
    assert q.value == 5.
    assert q.unit == u.m
    # Regression check to ensure we didn't create an object type by raising
    # the value of the quantity to a Fraction. [#3922]
    assert q.dtype.kind == 'f'


def test_inherit_docstrings():
    assert u.Quantity.argmax.__doc__ == np.ndarray.argmax.__doc__


def test_quantity_from_table():
    """
    Checks that units from tables are respected when converted to a Quantity.
    This also generically checks the use of *anything* with a `unit` attribute
    passed into Quantity
    """
    from... table import Table

    t = Table(data=[np.arange(5), np.arange(5)], names=['a', 'b'])
    t['a'].unit = u.kpc

    qa = u.Quantity(t['a'])
    assert qa.unit == u.kpc
    assert_array_equal(qa.value, t['a'])

    qb = u.Quantity(t['b'])
    assert qb.unit == u.dimensionless_unscaled
    assert_array_equal(qb.value, t['b'])

    # This does *not* auto-convert, because it's not necessarily obvious that's
    # desired.  Instead we revert to standard `Quantity` behavior
    qap = u.Quantity(t['a'], u.pc)
    assert qap.unit == u.pc
    assert_array_equal(qap.value, t['a'] * 1000)

    qbp = u.Quantity(t['b'], u.pc)
    assert qbp.unit == u.pc
    assert_array_equal(qbp.value, t['b'])


def test_assign_slice_with_quantity_like():
    # Regression tests for gh-5961
    from ... table import Table, Column
    # first check directly that we can use a Column to assign to a slice.
    c = Column(np.arange(10.), unit=u.mm)
    q = u.Quantity(c)
    q[:2] = c[:2]
    # next check that we do not fail the original problem.
    t = Table()
    t['x'] = np.arange(10) * u.mm
    t['y'] = np.ones(10) * u.mm
    assert type(t['x']) is Column

    xy = np.vstack([t['x'], t['y']]).T * u.mm
    ii = [0, 2, 4]

    assert xy[ii, 0].unit == t['x'][ii].unit
    # should not raise anything
    xy[ii, 0] = t['x'][ii]


def test_insert():
    """
    Test Quantity.insert method.  This does not test the full capabilities
    of the underlying np.insert, but hits the key functionality for
    Quantity.
    """
    q = [1, 2] * u.m

    # Insert a compatible float with different units
    q2 = q.insert(0, 1 * u.km)
    assert np.all(q2.value == [1000, 1, 2])
    assert q2.unit is u.m
    assert q2.dtype.kind == 'f'

    if minversion(np, '1.8.0'):
        q2 = q.insert(1, [1, 2] * u.km)
        assert np.all(q2.value == [1, 1000, 2000, 2])
        assert q2.unit is u.m

    # Cannot convert 1.5 * u.s to m
    with pytest.raises(u.UnitsError):
        q.insert(1, 1.5 * u.s)

    # Tests with multi-dim quantity
    q = [[1, 2], [3, 4]] * u.m
    q2 = q.insert(1, [10, 20] * u.m, axis=0)
    assert np.all(q2.value == [[1, 2],
                               [10, 20],
                               [3, 4]])

    q2 = q.insert(1, [10, 20] * u.m, axis=1)
    assert np.all(q2.value == [[1, 10, 2],
                               [3, 20, 4]])

    q2 = q.insert(1, 10 * u.m, axis=1)
    assert np.all(q2.value == [[1, 10, 2],
                               [3, 10, 4]])


def test_repr_array_of_quantity():
    """
    Test print/repr of object arrays of Quantity objects with different
    units.

    Regression test for the issue first reported in
    https://github.com/astropy/astropy/issues/3777
    """

    a = np.array([1 * u.m, 2 * u.s], dtype=object)
    if NUMPY_LT_1_14:
        assert repr(a) == 'array([<Quantity 1.0 m>, <Quantity 2.0 s>], dtype=object)'
        assert str(a) == '[<Quantity 1.0 m> <Quantity 2.0 s>]'
    else:
        assert repr(a) == 'array([<Quantity 1. m>, <Quantity 2. s>], dtype=object)'
        assert str(a) == '[<Quantity 1. m> <Quantity 2. s>]'


class TestSpecificTypeQuantity(object):
    def setup(self):
        class Length(u.SpecificTypeQuantity):
            _equivalent_unit = u.m

        class Length2(Length):
            _default_unit = u.m

        class Length3(Length):
            _unit = u.m

        self.Length = Length
        self.Length2 = Length2
        self.Length3 = Length3

    def test_creation(self):
        l = self.Length(np.arange(10.)*u.km)
        assert type(l) is self.Length
        with pytest.raises(u.UnitTypeError):
            self.Length(np.arange(10.) * u.hour)

        with pytest.raises(u.UnitTypeError):
            self.Length(np.arange(10.))

        l2 = self.Length2(np.arange(5.))
        assert type(l2) is self.Length2
        assert l2._default_unit is self.Length2._default_unit

        with pytest.raises(u.UnitTypeError):
            self.Length3(np.arange(10.))

    def test_view(self):
        l = (np.arange(5.) * u.km).view(self.Length)
        assert type(l) is self.Length
        with pytest.raises(u.UnitTypeError):
            (np.arange(5.) * u.s).view(self.Length)

        v = np.arange(5.).view(self.Length)
        assert type(v) is self.Length
        assert v._unit is None

        l3 = np.ones((2, 2)).view(self.Length3)
        assert type(l3) is self.Length3
        assert l3.unit is self.Length3._unit

    def test_operation_precedence_and_fallback(self):
        l = self.Length(np.arange(5.)*u.cm)
        sum1 = l + 1.*u.m
        assert type(sum1) is self.Length
        sum2 = 1.*u.km + l
        assert type(sum2) is self.Length
        sum3 = l + l
        assert type(sum3) is self.Length
        res1 = l * (1.*u.m)
        assert type(res1) is u.Quantity
        res2 = l * l
        assert type(res2) is u.Quantity


@pytest.mark.skipif('not HAS_MATPLOTLIB')
@pytest.mark.xfail('MATPLOTLIB_LT_15')
class TestQuantityMatplotlib(object):
    """Test if passing matplotlib quantities works.

    TODO: create PNG output and check against reference image
          once `astropy.wcsaxes` is merged, which provides
          the machinery for this.

    See https://github.com/astropy/astropy/issues/1881
    See https://github.com/astropy/astropy/pull/2139
    """

    def test_plot(self):
        data = u.Quantity([4, 5, 6], 's')
        plt.plot(data)

    def test_scatter(self):
        x = u.Quantity([4, 5, 6], 'second')
        y = [1, 3, 4] * u.m
        plt.scatter(x, y)


def test_unit_class_override():
    class MyQuantity(u.Quantity):
        pass

    my_unit = u.Unit("my_deg", u.deg)
    my_unit._quantity_class = MyQuantity
    q1 = u.Quantity(1., my_unit)
    assert type(q1) is u.Quantity
    q2 = u.Quantity(1., my_unit, subok=True)
    assert type(q2) is MyQuantity

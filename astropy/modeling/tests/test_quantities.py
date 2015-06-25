# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests specifically for models that use units and quantities."""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing import assert_allclose


from ..core import Fittable1DModel, ModelDefinitionError, InputParameterError
from ..parameters import Parameter, ParameterDefinitionError
from ..models import Gaussian1D
from ... import units as u
from ...units import UnitsError, Quantity
from ...tests.helper import pytest


def test_quantities_as_parameters():
    """
    Basic tests for initializing general models (that do not require units)
    with parameters that have units attached.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)
    assert g.amplitude.value == 1.0
    assert g.amplitude.unit is u.J
    assert g.mean.value == 1.0
    assert g.mean.unit is u.m
    assert g.stddev.value == 0.1
    assert g.stddev.unit is u.m


def test_quantity_parameter_descriptors():
    """
    Test Model classes that specify units in their Parameter descriptors,
    either via a Quantity default, or explicit use of the unit argument.
    """

    def tests(TestModel):
        assert TestModel.a.unit == u.m
        assert TestModel.a.default == 1.0

        m = TestModel()
        assert m.a.unit == u.m
        assert m.a.default == m.a.value == 1.0

        m = TestModel(2.0 * u.m)
        assert m.a.unit == u.m
        assert m.a.value == 2.0
        assert m.a.default == 1.0

        # Instantiate with a different, but compatible unit
        m = TestModel(2.0 * u.pc)
        assert m.a.unit == u.pc
        assert m.a.value == 2.0
        # The default is still in the original units
        assert m.a.default == 1.0

        # Instantiating with incompatible units is in error
        with pytest.raises(InputParameterError):
            TestModel(1.0 * u.Jy)

    class TestA(Fittable1DModel):
        a = Parameter(default=1.0, unit=u.m)
        @staticmethod
        def evaluate(x, a):
            return x

    tests(TestA)

    class TestB(Fittable1DModel):
        a = Parameter(default=1.0 * u.m)
        @staticmethod
        def evaluate(x, a):
            return x

    tests(TestB)

    # Conflicting default and units arguments
    with pytest.raises(ParameterDefinitionError):
        class TestC(Fittable1DModel):
            a = Parameter(default=1.0 * u.m, unit=u.Jy)
            @staticmethod
            def evaluate(x, a):
                return x


def test_quantity_parameter_arithmetic():
    """
    Test that arithmetic operations with properties that have units return the
    appropriate Quantities.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    assert g.mean + (1 * u.m) == 2 * u.m
    with pytest.raises(UnitsError):
        g.mean + 1
    assert (1 * u.m) + g.mean == 2 * u.m
    with pytest.raises(UnitsError):
        1 + g.mean
    assert g.mean * 2 == (2 * u.m)
    assert 2 * g.mean == (2 * u.m)
    assert g.mean * (2 * u.m) == (2 * (u.m ** 2))
    assert (2 * u.m) * g.mean == (2 * (u.m ** 2))

    assert -g.mean == (-1 * u.m)
    assert abs(-g.mean) == g.mean


def test_quantity_parameter_comparison():
    """
    Basic test of comparison operations on properties with units.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    assert g.mean == 1 * u.m
    assert 1 * u.m == g.mean
    assert g.mean != 1
    assert 1 != g.mean

    assert g.mean < 2 * u.m
    assert 2 * u.m > g.mean
    with pytest.raises(UnitsError):
        g.mean < 2
    with pytest.raises(UnitsError):
        2 > g.mean

    g = Gaussian1D([1, 2] * u.J, [1, 2] * u.m, [0.1, 0.2] * u.m)

    assert g.mean == [1, 2] * u.m
    assert np.all([1, 2] * u.m == g.mean)
    assert g.mean != [1, 2]
    assert np.all([1, 2] != g.mean)
    with pytest.raises(UnitsError):
        g.mean < [3, 4]
    with pytest.raises(UnitsError):
        [3, 4] > g.mean


def test_basic_evaluate_with_quantities():
    """
    Test evaluation of a single model with Quantity parameters, that does
    not explicitly require units.
    """

    g = Gaussian1D(1, 1, 0.1)
    gq = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    assert isinstance(gq(0), Quantity)
    assert gq(0).unit is u.J
    assert g(0) == gq(0).value

    # zero is allowed without explicit units, but other unitless quantities
    # should be an exception
    with pytest.raises(UnitsError):
        gq(1)

    assert gq(1 * u.m).value == g(1)

    # Should get the same numeric result as if we multiplied by 1000
    assert_allclose(gq(0.0005 * u.km).value, g(0.5))


def test_output_units():
    """
    Test multiple allowed settings for the output_units attribute
    on a single-output model.
    """

    class TestModelA(Fittable1DModel):
        a = Parameter()
        output_units = 'a'

        @staticmethod
        def evaluate(x, a):
            # For the sake of this test, the dimensions of x
            # shouldn't matter, and the return value should be
            # in the dimensions of the 'a' param
            return (x / x) * a

    m = TestModelA(a=1 * u.m)
    assert m(0).unit is m.a.unit
    assert m(1 * u.m).unit is m.a.unit
    assert m(27 * u.s).unit is m.a.unit

    class TestModelB(Fittable1DModel):
        a = Parameter()
        output_units = 'x'

        @staticmethod
        def evaluate(x, a):
            # In this cfase only the input units should matter
            return (a / a) * x

    m = TestModelB(a=1 / u.s)
    assert m(0).unit == u.dimensionless_unscaled
    assert m(1 * u.m).unit is u.m

    class TestModelC(Fittable1DModel):
        a = Parameter()
        output_units = lambda a, x: a.unit * x.unit

        @staticmethod
        def evaluate(x, a):
            # In this case the output's units are some compound
            # involving both the input and parameter's units
            return a * x

    m = TestModelC(a=1 / u.s)
    assert m(0).unit == (1 / u.s)
    assert m(2).unit == (1 / u.s)
    assert m(0 * u.dimensionless_unscaled).unit == (1 / u.s)
    assert m(2 * u.dimensionless_unscaled).unit == (1 / u.s)
    assert m(2 * u.m).unit == (u.m / u.s)

    class TestModelD(Fittable1DModel):
        output_units = u.m

        @staticmethod
        def evaluate(x):
            # This is a no-op model that just always forces the output to be in
            # meters (if the input is a length)
            return 2 * x

    m = TestModelD()
    assert m(0).unit is u.m
    assert m(1 * u.m) == 2 * u.m
    assert m(1 * u.km) == 2000 * u.m

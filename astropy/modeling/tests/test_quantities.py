# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests specifically for models that use units and quantities."""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing import assert_allclose


from ..core import Model, Fittable1DModel, ModelDefinitionError, InputParameterError
from ..parameters import Parameter, ParameterDefinitionError
from ..models import Gaussian1D
from ... import units as u
from ...units import UnitsError, Quantity
from ...tests.helper import pytest, assert_quantity_allclose


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


def test_parameter_unit_immutable():
    """
    Check that units can't be changed once set
    """

    g = Gaussian1D(1 * u.J, 3, 0.1)

    with pytest.raises(UnitsError) as exc:
        g.amplitude = 2
    assert exc.value.args[0] == "The 'amplitude' parameter should be given as a Quantity with units equivalent to J"

    with pytest.raises(UnitsError) as exc:
        g.amplitude = 2 * u.Jy
    assert exc.value.args[0] == "The 'amplitude' parameter should be given as a Quantity with units equivalent to J"

    with pytest.raises(UnitsError) as exc:
        g.mean = 2 * u.Jy
    assert exc.value.args[0] == "The 'mean' parameter should be given as a unitless value"

def test_parameter_unit_conversion():
    """
    Check that units are converted on-the-fly if compatible
    """

    g = Gaussian1D(1 * u.Jy, 3, 0.1)
    g.amplitude = 3000 * u.mJy
    assert g.amplitude.value == 3
    assert g.amplitude.unit is u.Jy


def test_parameter_unit_conversion_equivalencies():
    """
    Check that units are converted on-the-fly if compatible, including
    equivalencies.
    """

    class TestModel(Model):
        a = Parameter(default=1.0, unit=u.m)
        b = Parameter(default=1.0, unit=u.m, equivalencies=u.spectral())
        c = Parameter(default=1.0, unit=u.K, equivalencies=u.temperature_energy())
        @staticmethod
        def evaluate(x, a):
            return x

    model = TestModel()

    with pytest.raises(UnitsError) as exc:
        model.a = 3 * u.Hz
    assert exc.value.args[0] == "The 'a' parameter should be given as a Quantity with units equivalent to m"

    model.b = 3 * u.Hz

    with pytest.raises(UnitsError) as exc:
        model.c = 3 * u.Hz
    assert exc.value.args[0] == "The 'c' parameter should be given as a Quantity with units equivalent to K"

    model.c = 3 * u.eV


def test_init_equivalencies():
    """
    Check that equivalencies are taken into account when initializing a model
    with quantities, if using a model where the parameters have specific
    equivalencies defined.
    """

    class TestModel(Model):
        a = Parameter(default=1.0, unit=u.m)
        b = Parameter(default=1.0, unit=u.m, equivalencies=u.spectral())
        @staticmethod
        def evaluate(x, a):
            return x

    with pytest.raises(InputParameterError) as exc:
        model = TestModel(a=2 * u.Hz)
    assert exc.value.args[0] == ('TestModel.__init__() requires parameter \'a\' to '
                                 'be in units equivalent to Unit("m") (got Unit("Hz"))')

    model = TestModel(b=2 * u.Hz)

    with pytest.raises(InputParameterError) as exc:
        model = TestModel(b=2 * u.s)
    assert exc.value.args[0] == ('TestModel.__init__() requires parameter \'b\' to '
                                 'be in units equivalent to Unit("m") (got Unit("s"))')


def test_parameter_unit_equivalency():
    """
    Test that we can set equivalencies on an existing model
    """
    g = Gaussian1D(1 * u.Jy, 3 * u.m, 0.1 * u.nm)

    with pytest.raises(UnitsError) as exc:
        g.mean = 3 * u.Hz
    assert exc.value.args[0] == "The 'mean' parameter should be given as a Quantity with units equivalent to m"

    g.mean.equivalencies = u.spectral()

    g.mean = 3 * u.Hz


def test_parameter_change_unit():
    """
    Test that changing the unit on a parameter works as expected.
    """

    g = Gaussian1D(1, 1 * u.m, 0.1 * u.m)

    # Setting a unit on a unitless parameter should not work
    with pytest.raises(ValueError) as exc:
        g.amplitude.unit = u.Jy
    assert exc.value.args[0] == "Cannot attach units to parameters that were not initially specified with units"

    # Changing to an equivalent unit should work
    g.mean.unit = u.cm

    # But changing to another unit should not
    with pytest.raises(UnitsError) as exc:
        g.mean.unit = u.Hz
    assert exc.value.args[0] == ("Cannot set parameter units to Hz since it is "
                                 "not equivalent with the original units of cm")

    # Now we test a more complex example with equivalencies
    class TestModel(Model):
        a = Parameter(default=1.0, unit=u.m, equivalencies=u.spectral())
        @staticmethod
        def evaluate(x, a):
            return x

    model = TestModel()

    # This should work as before
    model.a.unit = u.cm

    # This should now work because of the equivalency
    model.a.unit = u.Hz

    # But this should still not work
    # But changing to another unit should not
    with pytest.raises(UnitsError) as exc:
        model.a.unit = u.Jy
    assert exc.value.args[0] == ("Cannot set parameter units to Jy since it is "
                                 "not equivalent with the original units of Hz")


def test_parameter_change_value():
    """
    Test that changing the value on a parameter works as expected.
    """

    g = Gaussian1D(1 * u.Jy, 1 * u.m, 0.1 * u.m)

    # To set a parameter to a quantity, we simply do
    g.amplitude = 2 * u.Jy

    # If we try setting the value, we need to pass a non-quantity value:
    g.amplitude.value = 4
    assert_quantity_allclose(g.amplitude, 4 * u.Jy)
    assert g.amplitude.value == 4
    assert g.amplitude.unit is u.Jy

    # If we try setting it to a Quantity, we raise an error
    with pytest.raises(TypeError) as exc:
        g.amplitude.value = 3 * u.Jy
    assert exc.value.args[0] == ("The .value property on parameters should be set to "
                                 "unitless values, not Quantity objects. To set a "
                                 "parameter to a quantity simply set the parameter "
                                 "directly without using .value")


def test_parameter_quantity_property():
    """
    Test that the quantity property of Parameters behaves as expected
    """

    g = Gaussian1D(1 * u.Jy, 1 * u.m, 0.1 * u.m)

    assert_quantity_allclose(g.amplitude.quantity, 1 * u.Jy)

    # Setting a parameter to a quantity changes the value and the default unit
    g.amplitude.quantity = 5 * u.mJy
    assert g.amplitude.value == 5
    assert g.amplitude.unit is u.mJy

    # But we shouldn't be able to set .quantity to a Quantity with
    # non-equivalent units
    with pytest.raises(UnitsError) as exc:
        g.amplitude.quantity = 3 * u.m
    assert exc.value.args[0] == ("Cannot set parameter units to m since it is not "
                                 "equivalent with the original units of mJy")


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

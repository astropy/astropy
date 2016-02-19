# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests that relate to using quantities/units on parameters of models.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ..core import Model, Fittable1DModel, InputParameterError
from ..parameters import Parameter, ParameterDefinitionError
from ..models import Gaussian1D
from ... import units as u
from ...units import UnitsError
from ...tests.helper import pytest, assert_quantity_allclose

class BaseTestModel(Fittable1DModel):
    @staticmethod
    def evaluate(x, a):
        return x


def test_parameter_quantity():
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


@pytest.mark.xfail
def test_parameter_set_compatible_units():
    """
    Make sure that parameters can be set to values with compatible units
    """
    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)
    g.amplitude = 4 * u.kJ
    assert_quantity_allclose(g.amplitude, 4 * u.kJ)
    g.mean = 3 * u.km
    assert_quantity_allclose(g.mean, 3 * u.km)
    g.stdev = 2 * u.mm
    assert_quantity_allclose(g.stddev, 2 * u.mm)


def test_parameter_unit_conversion():
    """
    Check that units are converted on-the-fly if compatible. Note that this is
    a slightly different test to above, because this is checking that the
    quantity is converted, whether the above test makes no such assumption.
    """
    g = Gaussian1D(1 * u.Jy, 3, 0.1)
    g.amplitude = 3000 * u.mJy
    assert g.amplitude.value == 3
    assert g.amplitude.unit is u.Jy


def test_parameter_set_incompatible_units():
    """
    Check that parameters can only be set to quantities with equivalent units.
    """

    g = Gaussian1D(1 * u.Jy, 3, 0.1)

    # The amplitude should be equivalent to flux density units, but in the
    # examples below, this is not the case.

    with pytest.raises(UnitsError) as exc:
        g.amplitude = 2
    assert exc.value.args[0] == ("The 'amplitude' parameter should be given as "
                                 "a Quantity with units equivalent to Jy")

    with pytest.raises(UnitsError) as exc:
        g.amplitude = 2 * u.s
    assert exc.value.args[0] == ("The 'amplitude' parameter should be given as "
                                 "a Quantity with units equivalent to Jy")

    with pytest.raises(UnitsError) as exc:
        g.mean = 2 * u.Jy
    assert exc.value.args[0] == ("The 'mean' parameter should be given as a "
                                 "unitless value")


def test_parameter_set_incompatible_units_with_equivalencies():
    """
    Check that parameters can only be set to quantities with equivalent units,
    taking into account equivalencies defined in the parameter initializers.
    """

    class TestModel(BaseTestModel):
        a = Parameter(default=1.0, unit=u.m)
        b = Parameter(default=1.0, unit=u.m, equivalencies=u.spectral())
        c = Parameter(default=1.0, unit=u.K, equivalencies=u.temperature_energy())

    model = TestModel()

    with pytest.raises(UnitsError) as exc:
        model.a = 3 * u.Hz
    assert exc.value.args[0] == "The 'a' parameter should be given as a Quantity with units equivalent to m"

    model.b = 3 * u.Hz

    with pytest.raises(UnitsError) as exc:
        model.c = 3 * u.Hz
    assert exc.value.args[0] == "The 'c' parameter should be given as a Quantity with units equivalent to K"

    model.c = 3 * u.eV


def test_parameter_equivalencies_used_in_init():
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
        TestModel(b=2 * u.s)
    assert exc.value.args[0] == ('TestModel.__init__() requires parameter \'b\' to '
                                 'be in units equivalent to Unit("m") (got Unit("s"))')


def test_parameter_set_equivalency():
    """
    Test that we can set equivalencies on an existing model.
    """
    g = Gaussian1D(1 * u.Jy, 3 * u.m, 0.1 * u.nm)

    with pytest.raises(UnitsError) as exc:
        g.mean = 3 * u.Hz
    assert exc.value.args[0] == ("The 'mean' parameter should be given as a "
                                 "Quantity with units equivalent to m")

    g.mean.equivalencies = u.spectral()

    g.mean = 3 * u.Hz

    # Note that we can't use this for equivalencies that rely on the values at
    # which the model is evaluated (such as brightness_temperature which
    # depends on spectral coordinate) because we won't know what values to use
    # until evaluation. See the test_output_units_equivalencies_with_parameters
    # test for how we can deal with this during evaluation.


def test_parameter_change_unit():
    """
    Test that changing the unit on a parameter works as expected (units can
    only be changed to a compatible unit).
    """

    g = Gaussian1D(1, 1 * u.m, 0.1 * u.m)

    # Setting a unit on a unitless parameter should not work
    with pytest.raises(ValueError) as exc:
        g.amplitude.unit = u.Jy
    assert exc.value.args[0] == ("Cannot attach units to parameters that were "
                                 "not initially specified with units")

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
    with pytest.raises(UnitsError) as exc:
        model.a.unit = u.Jy
    assert exc.value.args[0] == ("Cannot set parameter units to Jy since it is "
                                 "not equivalent with the original units of Hz")


def test_parameter_set_value():
    """
    Test that changing the value on a parameter works as expected.
    """

    g = Gaussian1D(1 * u.Jy, 1 * u.m, 0.1 * u.m)

    # To set a parameter to a quantity, we simply do
    g.amplitude = 2 * u.Jy

    # If we try setting the value, we need to pass a non-quantity value
    # TODO: determine whether this is the desired behavior?
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

    # Since parameters have a .value and .unit parameter that return just the
    # value and unit respectively, we also have a .quantity parameter that
    # returns a Quantity instance.

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


def test_parameter_default_units_match():

    # If the unit and default quantity units are different, raise an error
    with pytest.raises(ParameterDefinitionError) as exc:
        class TestC(Fittable1DModel):
            a = Parameter(default=1.0 * u.m, unit=u.Jy)
    assert exc.value.args[0] == ("parameter default 1.0 m does not have units "
                                 "equivalent to the required unit Jy")



@pytest.mark.parametrize(('unit', 'default'), ((u.m, 1.0), (None, 1 * u.m)))
def test_parameter_defaults(unit, default):
    """
    Test that default quantities are correctly taken into account
    """

    class TestModel(BaseTestModel):
        a = Parameter(default=default, unit=unit)

    # TODO: decide whether the default property should return a value or
    # a quantity?

    # The default unit and value should be set on the class
    assert TestModel.a.unit == u.m
    assert TestModel.a.default == 1.0

    # Check that the default unit and value are also set on a class instance
    m = TestModel()
    assert m.a.unit == u.m
    assert m.a.default == m.a.value == 1.0

    # If the parameter is set to a different value, the default is still the
    # internal default
    m = TestModel(2.0 * u.m)
    assert m.a.unit == u.m
    assert m.a.value == 2.0
    assert m.a.default == 1.0

    # Instantiate with a different, but compatible unit
    m = TestModel(2.0 * u.pc)
    assert m.a.unit == u.pc
    assert m.a.value == 2.0
    # The default is still in the original units
    # TODO: but how do we know what those units are if we don't return a
    # quantity?
    assert m.a.default == 1.0

    # Instantiating with incompatible units raises an error
    with pytest.raises(InputParameterError) as exc:
        TestModel(1.0 * u.Jy)
    assert exc.value.args[0] == ('TestModel.__init__() requires parameter \'a\' '
                                 'to be in units equivalent to Unit("m") '
                                 '(got Unit("Jy"))')


def test_parameter_quantity_arithmetic():
    """
    Test that arithmetic operations with properties that have units return the
    appropriate Quantities.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    # Addition should work if units are compatible
    assert g.mean + (1 * u.m) == 2 * u.m
    assert (1 * u.m) + g.mean == 2 * u.m


    # Multiplication by a scalar should also preserve the quantity-ness
    assert g.mean * 2 == (2 * u.m)
    assert 2 * g.mean == (2 * u.m)

    # Multiplication by a quantity should result in units being multiplied
    assert g.mean * (2 * u.m) == (2 * (u.m ** 2))
    assert (2 * u.m) * g.mean == (2 * (u.m ** 2))

    # Negation should work properly too
    assert -g.mean == (-1 * u.m)
    assert abs(-g.mean) == g.mean

    # However, addition of a quantity + scalar should not work
    with pytest.raises(UnitsError) as exc:
        g.mean + 1
    assert exc.value.args[0] == ("Can only apply 'add' function to "
                                 "dimensionless quantities when other argument "
                                 "is not a quantity (unless the latter is all "
                                 "zero/infinity/nan)")
    with pytest.raises(UnitsError) as exc:
        1 + g.mean
    assert exc.value.args[0] == ("Can only apply 'add' function to "
                                 "dimensionless quantities when other argument "
                                 "is not a quantity (unless the latter is all "
                                 "zero/infinity/nan)")


def test_parameter_quantity_comparison():
    """
    Basic test of comparison operations on properties with units.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    # Essentially here we are checking that parameters behave like Quantity

    assert g.mean == 1 * u.m
    assert 1 * u.m == g.mean
    assert g.mean != 1
    assert 1 != g.mean

    assert g.mean < 2 * u.m
    assert 2 * u.m > g.mean

    with pytest.raises(UnitsError) as exc:
        g.mean < 2
    assert exc.value.args[0] == ("Can only apply 'less' function to "
                                 "dimensionless quantities when other argument "
                                 "is not a quantity (unless the latter is all "
                                 "zero/infinity/nan)")

    with pytest.raises(UnitsError) as exc:
        2 > g.mean
    assert exc.value.args[0] == ("Can only apply 'less' function to "
                                 "dimensionless quantities when other argument "
                                 "is not a quantity (unless the latter is all "
                                 "zero/infinity/nan)")

    g = Gaussian1D([1, 2] * u.J, [1, 2] * u.m, [0.1, 0.2] * u.m)

    assert g.mean == [1, 2] * u.m
    assert np.all([1, 2] * u.m == g.mean)
    assert g.mean != [1, 2]
    assert np.all([1, 2] != g.mean)

    with pytest.raises(UnitsError) as exc:
        g.mean < [3, 4]
    assert exc.value.args[0] == ("Can only apply 'less' function to "
                                 "dimensionless quantities when other argument "
                                 "is not a quantity (unless the latter is all "
                                 "zero/infinity/nan)")

    with pytest.raises(UnitsError) as exc:
        [3, 4] > g.mean
    assert exc.value.args[0] == ("Can only apply 'less' function to "
                                 "dimensionless quantities when other argument "
                                 "is not a quantity (unless the latter is all "
                                 "zero/infinity/nan)")


def test_parameter_unit_linking_gaussian():
    """
    For some models, units of parameters should be related/consistent. We check
    the case of the Gaussian here, where the mean and stddev should have the
    same units.
    """

    with pytest.raises(UnitsError) as exc:
        g = Gaussian1D(1 * u.J, 1. * u.m, 0.1 * u.s)
    assert exc.value.args[0] == ("Parameters 'mean' and 'stddev' should "
                                 "have matching units")


@pytest.mark.xfail
def test_parameter_unit_linking_custom():

    # TODO: this is just an idea for an API, and open to feedback on the best
    # way to do this. The current approach is to provide a list of tuples,
    # where each tuple is a group of units that should have linked units. I
    # found it easier conceptually to keep this outside the parameter
    # definition, otherwise we have to navigate all parameters and construct a
    # tree.

    with pytest.raises(ModelDefinitionError) as exc:

        class TestModel(BaseTestModel):

            linked_units = [('a', 'b')]

            a = Parameter(unit=u.m)
            b = Parameter(unit=u.s)

    assert exc.value.args[0] == ("Parameters 'a' and 'b' should "
                                 "have matching units")

    # We now check what happens if units are set at initialization time

    class TestModel(BaseTestModel):

        linked_units = [('a', 'b')]

        a = Parameter()
        b = Parameter()
        c = Parameter()

    TestModel(1 * u.m, 2 * u.cm, 3 * u.kg)

    with pytest.raises(UnitsError) as exc:
        TestModel(1 * u.m, 2 * u.s, 3 * u.kg)
    assert exc.value.args[0] == ("Parameters 'a' and 'b' should "
                                 "have matching units")

    # It's also important that parameters don't appear more than once in the
    # list of linked units (otherwise one could combine the groups anyway).
    with pytest.raises(ModelDefinitionError) as exc:

        class TestModel(BaseTestModel):

            linked_units = [('a', 'b'), ('a', 'c')]

            a = Parameter(unit=u.m)
            b = Parameter(unit=u.km)
            c = Parameter(unit=u.mm)

    assert exc.value.args[0] == ("Parameter 'a' appears multiple times in 'linked_units'")

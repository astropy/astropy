# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests that relate to evaluating models with quantity parameters
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from numpy.testing import assert_allclose


from ..core import Fittable1DModel
from ..parameters import Parameter
from ..models import Gaussian1D
from ... import units as u
from ...units import UnitsError
from ...tests.helper import pytest, assert_quantity_allclose


@pytest.mark.xfail
def test_evaluate_with_quantities():
    """
    Test evaluation of a single model with Quantity parameters that does
    not explicitly require units.
    """

    # We create two models here - one with quantities, and one without. The one
    # without is used to create the reference values for comparison.

    g = Gaussian1D(1, 1, 0.1)
    gq = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    # We first check that calling the Gaussian with quantities returns the
    # expected result
    assert_quantity_allclose(gq(1 * u.m), g(1) * u.J)

    # Units have to be specified for the Gaussian with quantities - if not, an
    # error is raised
    with pytest.raises(UnitsError) as exc:
        gq(1)
    assert exc.value.args[0] == ("Units of input 'x', (dimensionless), could not be "
                                 "converted to required input units of m (length)")

    # However, zero is a special case
    assert_quantity_allclose(gq(0), g(0) * u.J)

    # We can also evaluate models with equivalent units
    assert_allclose(gq(0.0005 * u.km).value, g(0.5))

    # But not with incompatible units
    with pytest.raises(UnitsError) as exc:
        gq(3 * u.s)
    assert exc.value.args[0] == ("Units of input 'x', s (time), could not be "
                                 "converted to required input units of m (length)")

    # We also can't evaluate the model without quantities with a quantity
    with pytest.raises(UnitsError) as exc:
        g(3 * u.m)
    assert exc.value.args[0] == ("Input in m (length), could not be converted "
                                 "to required to dimensionless input")


@pytest.mark.xfail
def test_evaluate_with_quantities_and_equivalencies():
    """
    We now make sure that equivalencies are correctly taken into account
    """

    g = Gaussian1D(1 * u.Jy, 10 * u.nm, 2 * u.nm)

    # We haven't set the equivalencies yet, so this won't work
    with pytest.raises(UnitsError) as exc:
        g(30 * u.PHz)
    assert exc.value.args[0] == ("Units of input 'x', PHz (frequency), could "
                                 "not be converted to required input units of "
                                 "nm (length)")

    g.mean.equivalencies = u.spectral()
    g.stddev.equivalencies = u.spectral()

    # But it should now work
    assert_quantity_allclose(g(30 * u.PHz), g(9.993081933333332 * u.nm))


# We now have a series of tests to test the functionality of setting the
# output_units attribute on models to what we want the acceptable outputs to be.
# This can be set to the name of another parameter, or 'x' for the input values,
# and can also be set to a lambda function or a fixed unit. All output is
# converted to these units on-the-fly, and an error is raised if this is not
# possible. We test all the different cases below.


def test_output_units_link_parameter():
    """
    Test setting output_units to match one of the parameters
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

    # The output units always match the input units
    assert m(0).unit is u.m
    assert m(1 * u.m).unit is u.m
    assert m(27 * u.s).unit is u.m


@pytest.mark.xfail
def test_output_units_link_parameter_inconsistent():
    """
    Test setting output_units to match one of the parameters, in the case
    where the evaluate method returns inconsistent units.
    """

    class TestModelA(Fittable1DModel):
        a = Parameter()
        output_units = 'a'

        @staticmethod
        def evaluate(x, a):
            return x * a

    m = TestModelA(a=1 * u.m)

    assert m(5).unit is u.m
    with pytest.raises(UnitsError) as exc:
        m(3 * u.m)
    assert exc.value.args[0] == ("Units of output 'y', m2 (area), could not be "
                                 "converted to required output units of m (length)")


def test_output_units_link_input():
    """
    Test setting output_units to match the input.
    """

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


@pytest.mark.xfail
def test_output_units_link_input_inconsistent():
    """
    Test setting output_units to match the input, in the case where the evaluate
    method returns inconsistent units.
    """

    class TestModelB(Fittable1DModel):
        a = Parameter()
        output_units = 'x'

        @staticmethod
        def evaluate(x, a):
            # In this cfase only the input units should matter
            return a * x

    m = TestModelB(a=1 / u.s)
    with pytest.raises(UnitsError) as exc:
        m(3 * u.m)
    assert exc.value.args[0] == ("Units of output 'y', m / s (speed), could not be "
                                 "converted to required output units of m (length)")


def test_output_units_functional():
    """
    Test setting output_units to a custom function.
    """

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


@pytest.mark.xfail
def test_output_units_functional_inconsistent():
    """
    Test setting output_units to a custom function that returns units
    inconsistent with the output of the evaluate method.
    """

    class TestModelC(Fittable1DModel):
        a = Parameter()
        output_units = lambda a, x: a.unit * x.unit

        @staticmethod
        def evaluate(x, a):
            # In this case the output's units are some compound
            # involving both the input and parameter's units
            return a * x * u.s

    m = TestModelC(a=1 / u.s)
    with pytest.raises(UnitsError) as exc:
        m(3 * u.m)
    assert exc.value.args[0] == ("Units of output 'y', m (length), could not be "
                                 "converted to required output units of m / s "
                                 "(speed)")


def test_output_units_fixed_unit():
    """
    Test setting output_units to a fixed unit
    """

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


def test_output_units_fixed_unit_inconsistent():
    """
    Test setting output_units to a fixed unit that is inconsistent with
    evaluate method output.
    """

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


def test_output_units_instance_specific():
    """
    Test changing output_units on-the-fly for a specific instance
    """

    class TestModelD(Fittable1DModel):

        output_units = u.m

        @staticmethod
        def evaluate(x):
            # This is a no-op model that just always forces the output to be in
            # meters (if the input is a length)
            return 2 * x

    m1 = TestModelD()
    m2 = TestModelD()
    m2.output_units = u.s

    assert m1(1 * u.m) == 2 * u.m
    with pytest.raises(UnitsError) as exc:
        m2(1 * u.m)
    assert exc.value.args[0] == ("Units of input 'y', m (length), could not be "
                                 "converted to required input units of s (time)")
    assert m2(1 * u.s) == 2 * u.s


# We now have a series of similar tests for the input_units attribute


def test_input_units_link_parameter():
    """
    Test setting input_units to match one of the parameters
    """

    class TestModelA(Fittable1DModel):
        a = Parameter()
        input_units = 'a'
        @staticmethod
        def evaluate(x, a):
            return a * x

    m = TestModelA(a=1 * u.m)

    with pytest.raises(UnitsError) as exc:
        m(1)
    assert exc.value.args[0] == ("Units of input 'x', (dimensionless), could "
                                 "not be converted to required input units of m "
                                 "(length)")

    with pytest.raises(UnitsError) as exc:
        m(1 * u.s)
    assert exc.value.args[0] == ("Units of input 'x', s (time), could "
                                 "not be converted to required input units of m "
                                 "(length)")

    assert m(2 * u.m) == 2 * u.m ** 2

def test_input_units_link_input():
    """
    Test setting input_units to match the input essentially places no
    constraints since that would be true by definition.
    """

    class TestModelB(Fittable1DModel):
        a = Parameter()
        input_units = 'x'
        @staticmethod
        def evaluate(x, a):
            return a * x

    m = TestModelB(a=1 / u.s)
    assert m(1) == 1 / u.s
    assert m(3 * u.m) == 3 * u.m / u.s
    assert m(4 * u.s) == 4


def test_input_units_functional():
    """
    Test setting input_units to a custom function.
    """

    class TestModelC(Fittable1DModel):
        a = Parameter()
        input_units = lambda a: a.unit * u.m
        @staticmethod
        def evaluate(x, a):
            return a * x * u.s ** 2

    m = TestModelC(a=1 / u.s)
    assert m(3 * u.m / u.s) == 3 * u.m

    with pytest.raises(UnitsError) as exc:
        m(1 * u.m)
    assert exc.value.args[0] == ("Units of input 'x', m (length), could "
                                 "not be converted to required input units of "
                                 "m / s (speed)")


def test_input_units_fixed_unit():
    """
    Test setting input_units to a fixed unit
    """

    class TestModelD(Fittable1DModel):
        input_units = u.m
        @staticmethod
        def evaluate(x):
            return 2 * x

    m = TestModelD()

    with pytest.raises(UnitsError) as exc:
        m(1 * u.s)
    assert exc.value.args[0] == ("Units of input 'x', s (time), could "
                                 "not be converted to required input units of "
                                 "m (length)")

    assert m(2 * u.m) == 4 * u.m
    assert_quantity_allclose(m(300 * u.cm), 6 * u.m)


def test_input_units_instance_specific():
    """
    Test changing input_units on-the-fly for a specific instance
    """

    class TestModelD(Fittable1DModel):
        input_units = u.m
        @staticmethod
        def evaluate(x):
            return 2 * x

    m1 = TestModelD()
    m2 = TestModelD()
    m2.input_units = u.s

    assert m1(1 * u.m) == 2 * u.m
    with pytest.raises(UnitsError) as exc:
        m2(1 * u.m)
    assert exc.value.args[0] == ("Units of input 'x', m (length), could not be "
                                 "converted to required input units of s (time)")
    assert m2(1 * u.s) == 2 * u.s


# We now test that we are able to set equivalencies for input/output units. We
# don't need to repeat all the above tests, we simply choose the case where the
# input/output units are fixed.

@pytest.mark.xfail
def test_input_units_equivalencies():

    class TestModel(Fittable1DModel):
        input_units = u.m
        input_equivalencies = u.spectral()
        @staticmethod
        def evaluate(x):
            return 2 * x

    m = TestModel()
    assert_quantity_allclose(m(3 * u.nm), 6 * u.nm)

    # We now check it works with frequency. Note that we should make sure
    # we evaluate the model outside the context manager to make sure that we
    # are picking up the equivalency from input_equivalencies. The context
    # manager is only there because assert_quantity_allclose doesn't support
    # equivalencies.
    result = m(400 * u.GHz)
    with u.set_enabled_equivalencies(u.spectral()):
        assert_quantity_allclose(result, 800 * u.GHz)


@pytest.mark.xfail
def test_output_units_equivalencies():

    class TestModel(Fittable1DModel):
        output_units = u.m
        output_equivalencies = u.spectral()
        @staticmethod
        def evaluate(x):
            return 2 * x

    m = TestModel()
    assert_quantity_allclose(m(3 * u.nm), 6 * u.nm)

    # We now check it works with frequency. As before, we need to evaluate the
    # model outside the context manager.
    result = m(400 * u.GHz)

    # The output units should actually be meters since that is what was
    # requested - we're just checking that no error was raised above.
    assert result.unit is u.m

    # And we now check the actual accuracy of the result to be sure
    with u.set_enabled_equivalencies(u.spectral()):
        assert_quantity_allclose(result, 800 * u.GHz)

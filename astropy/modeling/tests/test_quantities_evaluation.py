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
def test_evaluate_with_quantities_with_equivalencies():
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


def test_evaluate_output_units():
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

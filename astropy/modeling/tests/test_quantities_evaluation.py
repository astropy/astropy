# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests that relate to evaluating models with quantity parameters
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
import pytest
from numpy.testing import assert_allclose


from ..core import Model
from ..models import Gaussian1D
from ... import units as u
from ...units import UnitsError
from ...tests.helper import assert_quantity_allclose


# We start off by taking some simple cases where the units are defined by
# whatever the model is initialized with, and we check that the model evaluation
# returns quantities.


def test_evaluate_with_quantities():
    """
    Test evaluation of a single model with Quantity parameters that do
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
    # TODO: determine what error message should be here
    # assert exc.value.args[0] == ("Units of input 'x', m (length), could not be "
    #                              "converted to required dimensionless input")


def test_evaluate_with_quantities_and_equivalencies():
    """
    We now make sure that equivalencies are correctly taken into account
    """

    g = Gaussian1D(1 * u.Jy, 10 * u.nm, 2 * u.nm)

    # We aren't setting the equivalencies, so this won't work
    with pytest.raises(UnitsError) as exc:
        g(30 * u.PHz)
    assert exc.value.args[0] == ("Units of input 'x', PHz (frequency), could "
                                 "not be converted to required input units of "
                                 "nm (length)")

    # But it should now work if we pass equivalencies when evaluating
    assert_quantity_allclose(g(30 * u.PHz, equivalencies={'x': u.spectral()}),
                             g(9.993081933333332 * u.nm))


class MyTestModel(Model):
    inputs = ('a', 'b')
    outputs = ('f',)

    def evaluate(self, a, b):
        print('a', a)
        print('b', b)
        return a * b


class TestInputUnits():

    def setup_method(self, method):
        self.model = MyTestModel()

    def test_evaluate(self):
        # We should be able to evaluate with anything
        assert_quantity_allclose(self.model(3, 5), 15)
        assert_quantity_allclose(self.model(4 * u.m, 5), 20 * u.m)
        assert_quantity_allclose(self.model(3 * u.deg, 5), 15 * u.deg)

    def test_input_units(self):

        self.model.input_units = {'a': u.deg}

        assert_quantity_allclose(self.model(3 * u.deg, 4), 12 * u.deg)
        assert_quantity_allclose(self.model(4 * u.rad, 2), 8 * u.rad)
        assert_quantity_allclose(self.model(4 * u.rad, 2 * u.s), 8 * u.rad * u.s)

        with pytest.raises(UnitsError) as exc:
            self.model(4 * u.s, 3)
        assert exc.value.args[0] == ("Units of input 'a', s (time), could not be "
                                     "converted to required input units of deg (angle)")

        with pytest.raises(UnitsError) as exc:
            self.model(3, 3)
        assert exc.value.args[0] == ("Units of input 'a', (dimensionless), could "
                                     "not be converted to required input units of deg (angle)")

    def test_input_units_allow_dimensionless(self):

        self.model.input_units = {'a': u.deg}
        self.model.input_units_allow_dimensionless = True

        assert_quantity_allclose(self.model(3 * u.deg, 4), 12 * u.deg)
        assert_quantity_allclose(self.model(4 * u.rad, 2), 8 * u.rad)

        with pytest.raises(UnitsError) as exc:
            self.model(4 * u.s, 3)
        assert exc.value.args[0] == ("Units of input 'a', s (time), could not be "
                                     "converted to required input units of deg (angle)")

        assert_quantity_allclose(self.model(3, 3), 9)

    def test_input_units_strict(self):

        self.model.input_units = {'a': u.deg}
        self.model.input_units_strict = True

        assert_quantity_allclose(self.model(3 * u.deg, 4), 12 * u.deg)

        result = self.model(np.pi * u.rad, 2)
        assert_quantity_allclose(result, 360 * u.deg)
        assert result.unit is u.deg

    def test_input_units_equivalencies(self):

        self.model.input_units = {'a': u.micron}

        with pytest.raises(UnitsError) as exc:
            self.model(3 * u.PHz, 3)
        assert exc.value.args[0] == ("Units of input 'a', PHz (frequency), could "
                                     "not be converted to required input units of "
                                     "micron (length)")

        self.model.input_units_equivalencies = {'a': u.spectral()}

        assert_quantity_allclose(self.model(3 * u.PHz, 3),
                                3 * (3 * u.PHz).to(u.micron, equivalencies=u.spectral()))

    def test_return_units(self):

        self.model.input_units = {'a': u.deg}
        self.model.return_units = {'f': u.rad}

        result = self.model(3 * u.deg, 4)

        assert_quantity_allclose(result, 12 * u.deg)
        assert result.unit is u.rad

    def test_return_units_scalar(self):

        # Check that return_units also works when giving a single unit since
        # there is only one output, so is unambiguous.

        self.model.input_units = {'a': u.deg}
        self.model.return_units = u.rad

        result = self.model(3 * u.deg, 4)

        assert_quantity_allclose(result, 12 * u.deg)
        assert result.unit is u.rad

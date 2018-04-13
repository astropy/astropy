# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests that relate to evaluating models with quantity parameters
"""


import numpy as np
import pytest
from numpy.testing import assert_allclose


from ..core import Model
from ..models import Gaussian1D, Shift, Scale, Pix2Sky_TAN
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
    assert exc.value.args[0] == ("Gaussian1D: Units of input 'x', (dimensionless), could not be "
                                 "converted to required input units of m (length)")

    # However, zero is a special case
    assert_quantity_allclose(gq(0), g(0) * u.J)

    # We can also evaluate models with equivalent units
    assert_allclose(gq(0.0005 * u.km).value, g(0.5))

    # But not with incompatible units
    with pytest.raises(UnitsError) as exc:
        gq(3 * u.s)
    assert exc.value.args[0] == ("Gaussian1D: Units of input 'x', s (time), could not be "
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
    assert exc.value.args[0] == ("Gaussian1D: Units of input 'x', PHz (frequency), could "
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

        self.model._input_units = {'a': u.deg}

        assert_quantity_allclose(self.model(3 * u.deg, 4), 12 * u.deg)
        assert_quantity_allclose(self.model(4 * u.rad, 2), 8 * u.rad)
        assert_quantity_allclose(self.model(4 * u.rad, 2 * u.s), 8 * u.rad * u.s)

        with pytest.raises(UnitsError) as exc:
            self.model(4 * u.s, 3)
        assert exc.value.args[0] == ("MyTestModel: Units of input 'a', s (time), could not be "
                                     "converted to required input units of deg (angle)")

        with pytest.raises(UnitsError) as exc:
            self.model(3, 3)
        assert exc.value.args[0] == ("MyTestModel: Units of input 'a', (dimensionless), could "
                                     "not be converted to required input units of deg (angle)")

    def test_input_units_allow_dimensionless(self):

        self.model._input_units = {'a': u.deg}
        self.model._input_units_allow_dimensionless = True

        assert_quantity_allclose(self.model(3 * u.deg, 4), 12 * u.deg)
        assert_quantity_allclose(self.model(4 * u.rad, 2), 8 * u.rad)

        with pytest.raises(UnitsError) as exc:
            self.model(4 * u.s, 3)
        assert exc.value.args[0] == ("MyTestModel: Units of input 'a', s (time), could not be "
                                     "converted to required input units of deg (angle)")

        assert_quantity_allclose(self.model(3, 3), 9)

    def test_input_units_strict(self):

        self.model._input_units = {'a': u.deg}
        self.model._input_units_strict = True

        assert_quantity_allclose(self.model(3 * u.deg, 4), 12 * u.deg)

        result = self.model(np.pi * u.rad, 2)
        assert_quantity_allclose(result, 360 * u.deg)
        assert result.unit is u.deg

    def test_input_units_equivalencies(self):

        self.model._input_units = {'a': u.micron}

        with pytest.raises(UnitsError) as exc:
            self.model(3 * u.PHz, 3)
        assert exc.value.args[0] == ("MyTestModel: Units of input 'a', PHz (frequency), could "
                                     "not be converted to required input units of "
                                     "micron (length)")

        self.model.input_units_equivalencies = {'a': u.spectral()}

        assert_quantity_allclose(self.model(3 * u.PHz, 3),
                                3 * (3 * u.PHz).to(u.micron, equivalencies=u.spectral()))

    def test_return_units(self):

        self.model._input_units = {'a': u.deg}
        self.model._return_units = {'f': u.rad}

        result = self.model(3 * u.deg, 4)

        assert_quantity_allclose(result, 12 * u.deg)
        assert result.unit is u.rad

    def test_return_units_scalar(self):

        # Check that return_units also works when giving a single unit since
        # there is only one output, so is unambiguous.

        self.model._input_units = {'a': u.deg}
        self.model._return_units = u.rad

        result = self.model(3 * u.deg, 4)

        assert_quantity_allclose(result, 12 * u.deg)
        assert result.unit is u.rad


def test_and_input_units():
    """
    Test units to first model in chain.
    """
    s1 = Shift(10 * u.deg)
    s2 = Shift(10 * u.deg)

    cs = s1 & s2

    out = cs(10 * u.arcsecond, 20 * u.arcsecond)

    assert_quantity_allclose(out[0], 10 * u.deg + 10 * u.arcsec)
    assert_quantity_allclose(out[1], 10 * u.deg + 20 * u.arcsec)


def test_plus_input_units():
    """
    Test units to first model in chain.
    """
    s1 = Shift(10 * u.deg)
    s2 = Shift(10 * u.deg)

    cs = s1 + s2

    out = cs(10 * u.arcsecond)

    assert_quantity_allclose(out, 20 * u.deg + 20 * u.arcsec)


def test_compound_input_units():
    """
    Test units to first model in chain.
    """
    s1 = Shift(10 * u.deg)
    s2 = Shift(10 * u.deg)

    cs = s1 | s2

    out = cs(10 * u.arcsecond)

    assert_quantity_allclose(out, 20 * u.deg + 10 * u.arcsec)


def test_compound_input_units_fail():
    """
    Test incompatible units to first model in chain.
    """
    s1 = Shift(10 * u.deg)
    s2 = Shift(10 * u.deg)

    cs = s1 | s2

    with pytest.raises(UnitsError):
        cs(10 * u.pix)


def test_compound_incompatible_units_fail():
    """
    Test incompatible model units in chain.
    """
    s1 = Shift(10 * u.pix)
    s2 = Shift(10 * u.deg)

    cs = s1 | s2

    with pytest.raises(UnitsError):
        cs(10 * u.pix)


def test_compound_pipe_equiv_call():
    """
    Check that equivalencies work when passed to evaluate, for a chained model
    (which has one input).
    """
    s1 = Shift(10 * u.deg)
    s2 = Shift(10 * u.deg)

    cs = s1 | s2

    out = cs(10 * u.pix, equivalencies={'x': u.pixel_scale(0.5 * u.deg / u.pix)})
    assert_quantity_allclose(out, 25 * u.deg)


def test_compound_and_equiv_call():
    """
    Check that equivalencies work when passed to evaluate, for a compsite model
    with two inputs.
    """
    s1 = Shift(10 * u.deg)
    s2 = Shift(10 * u.deg)

    cs = s1 & s2

    out = cs(10 * u.pix, 10 * u.pix, equivalencies={'x0': u.pixel_scale(0.5 * u.deg / u.pix),
                                                    'x1': u.pixel_scale(0.5 * u.deg / u.pix)})
    assert_quantity_allclose(out[0], 15 * u.deg)
    assert_quantity_allclose(out[1], 15 * u.deg)


def test_compound_input_units_equivalencies():
    """
    Test setting input_units_equivalencies on one of the models.
    """

    s1 = Shift(10 * u.deg)
    s1.input_units_equivalencies = {'x': u.pixel_scale(0.5 * u.deg / u.pix)}
    s2 = Shift(10 * u.deg)
    sp = Shift(10 * u.pix)

    cs = s1 | s2

    out = cs(10 * u.pix)
    assert_quantity_allclose(out, 25 * u.deg)

    cs = sp | s1

    out = cs(10 * u.pix)
    assert_quantity_allclose(out, 20 * u.deg)

    cs = s1 & s2
    cs = cs.rename('TestModel')
    out = cs(20 * u.pix, 10 * u.deg)
    assert_quantity_allclose(out, 20 * u.deg)

    with pytest.raises(UnitsError) as exc:
        out = cs(20 * u.pix, 10 * u.pix)
    assert exc.value.args[0] == "TestModel: Units of input 'x1', pix (unknown), could not be converted to required input units of deg (angle)"


def test_compound_input_units_strict():
    """
    Test setting input_units_strict on one of the models.
    """

    class ScaleDegrees(Scale):
        input_units = {'x': u.deg}

    s1 = ScaleDegrees(2)
    s2 = Scale(2)

    cs = s1 | s2

    out = cs(10 * u.arcsec)
    assert_quantity_allclose(out, 40 * u.arcsec)
    assert out.unit is u.deg  # important since this tests input_units_strict

    cs = s2 | s1

    out = cs(10 * u.arcsec)
    assert_quantity_allclose(out, 40 * u.arcsec)
    assert out.unit is u.deg  # important since this tests input_units_strict

    cs = s1 & s2

    out = cs(10 * u.arcsec, 10 * u.arcsec)
    assert_quantity_allclose(out, 20 * u.arcsec)
    assert out[0].unit is u.deg
    assert out[1].unit is u.arcsec


def test_compound_input_units_allow_dimensionless():
    """
    Test setting input_units_allow_dimensionless on one of the models.
    """

    class ScaleDegrees(Scale):
        input_units = {'x': u.deg}

    s1 = ScaleDegrees(2)
    s1._input_units_allow_dimensionless = True
    s2 = Scale(2)

    cs = s1 | s2
    cs = cs.rename('TestModel')
    out = cs(10)
    assert_quantity_allclose(out, 40 * u.one)

    out = cs(10 * u.arcsec)
    assert_quantity_allclose(out, 40 * u.arcsec)

    with pytest.raises(UnitsError) as exc:
        out = cs(10 * u.m)
    assert exc.value.args[0] == "TestModel: Units of input 'x', m (length), could not be converted to required input units of deg (angle)"

    s1._input_units_allow_dimensionless = False

    cs = s1 | s2
    cs = cs.rename('TestModel')

    with pytest.raises(UnitsError) as exc:
        out = cs(10)
    assert exc.value.args[0] == "TestModel: Units of input 'x', (dimensionless), could not be converted to required input units of deg (angle)"

    s1._input_units_allow_dimensionless = True

    cs = s2 | s1
    cs = cs.rename('TestModel')

    out = cs(10)
    assert_quantity_allclose(out, 40 * u.one)

    out = cs(10 * u.arcsec)
    assert_quantity_allclose(out, 40 * u.arcsec)

    with pytest.raises(UnitsError) as exc:
        out = cs(10 * u.m)
    assert exc.value.args[0] == "ScaleDegrees: Units of input 'x', m (length), could not be converted to required input units of deg (angle)"

    s1._input_units_allow_dimensionless = False

    cs = s2 | s1

    with pytest.raises(UnitsError) as exc:
        out = cs(10)
    assert exc.value.args[0] == "ScaleDegrees: Units of input 'x', (dimensionless), could not be converted to required input units of deg (angle)"

    s1._input_units_allow_dimensionless = True

    s1 = ScaleDegrees(2)
    s1._input_units_allow_dimensionless = True
    s2 = ScaleDegrees(2)
    s2._input_units_allow_dimensionless = False

    cs = s1 & s2
    cs = cs.rename('TestModel')

    out = cs(10, 10 * u.arcsec)
    assert_quantity_allclose(out[0], 20 * u.one)
    assert_quantity_allclose(out[1], 20 * u.arcsec)

    with pytest.raises(UnitsError) as exc:
        out = cs(10, 10)
    assert exc.value.args[0] == "TestModel: Units of input 'x1', (dimensionless), could not be converted to required input units of deg (angle)"


def test_compound_return_units():
    """
    Test that return_units on the first model in the chain is respected for the
    input to the second.
    """

    class PassModel(Model):

        inputs = ('x', 'y')
        outputs = ('x', 'y')

        @property
        def input_units(self):
            """ Input units. """
            return {'x': u.deg, 'y': u.deg}

        @property
        def return_units(self):
            """ Output units. """
            return {'x': u.deg, 'y': u.deg}

        def evaluate(self, x, y):
            return x.value, y.value


    cs = Pix2Sky_TAN() | PassModel()

    assert_quantity_allclose(cs(0*u.deg, 0*u.deg), (0, 90)*u.deg)

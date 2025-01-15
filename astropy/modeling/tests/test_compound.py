# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name, pointless-statement

import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import astropy.units as u
from astropy.modeling.core import CompoundModel, Model, ModelDefinitionError
from astropy.modeling.fitting import DogBoxLSQFitter, LevMarLSQFitter
from astropy.modeling.models import (
    Chebyshev1D,
    Chebyshev2D,
    Const1D,
    Gaussian1D,
    Gaussian2D,
    Identity,
    Legendre1D,
    Legendre2D,
    Linear1D,
    Mapping,
    Polynomial1D,
    Polynomial2D,
    Rotation2D,
    Scale,
    Shift,
    Tabular1D,
    fix_inputs,
)
from astropy.modeling.parameters import Parameter
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_SCIPY


@pytest.mark.parametrize(
    ("expr", "result"),
    [
        (lambda x, y: x + y, [5.0, 5.0]),
        (lambda x, y: x - y, [-1.0, -1.0]),
        (lambda x, y: x * y, [6.0, 6.0]),
        (lambda x, y: x / y, [2.0 / 3.0, 2.0 / 3.0]),
        (lambda x, y: x**y, [8.0, 8.0]),
    ],
)
def test_model_set(expr, result):
    s = expr(Const1D((2, 2), n_models=2), Const1D((3, 3), n_models=2))
    out = s(0, model_set_axis=False)
    assert_array_equal(out, result)


@pytest.mark.parametrize(
    ("expr", "result"),
    [
        (lambda x, y: x + y, [5.0, 5.0]),
        (lambda x, y: x - y, [-1.0, -1.0]),
        (lambda x, y: x * y, [6.0, 6.0]),
        (lambda x, y: x / y, [2.0 / 3.0, 2.0 / 3.0]),
        (lambda x, y: x**y, [8.0, 8.0]),
    ],
)
def test_model_set_raises_value_error(expr, result):
    """Check that creating model sets with components whose _n_models are
    different raise a value error
    """
    MESSAGE = r"Both operands must have equal values for .*"
    with pytest.raises(ValueError, match=MESSAGE):
        expr(Const1D((2, 2), n_models=2), Const1D(3, n_models=1))


@pytest.mark.parametrize(
    ("expr", "result"),
    [
        (lambda x, y: x + y, 5.0),
        (lambda x, y: x - y, -1.0),
        (lambda x, y: x * y, 6.0),
        (lambda x, y: x / y, 2.0 / 3.0),
        (lambda x, y: x**y, 8.0),
    ],
)
def test_two_model_instance_arithmetic_1d(expr, result):
    """
    Like test_two_model_class_arithmetic_1d, but creates a new model from two
    model *instances* with fixed parameters.
    """

    s = expr(Const1D(2), Const1D(3))

    assert isinstance(s, CompoundModel)
    assert s.n_inputs == 1
    assert s.n_outputs == 1

    out = s(0)
    assert out == result

    assert isinstance(out, float)


def test_simple_two_model_compose_1d():
    """
    Shift and Scale are two of the simplest models to test model composition
    with.
    """

    S1 = Shift(2) | Scale(3)  # First shift then scale
    assert isinstance(S1, CompoundModel)
    assert S1.n_inputs == 1
    assert S1.n_outputs == 1
    assert S1(1) == 9.0

    S2 = Scale(2) | Shift(3)  # First scale then shift
    assert isinstance(S2, CompoundModel)
    assert S2.n_inputs == 1
    assert S2.n_outputs == 1
    assert S2(1) == 5.0

    # Test with array inputs
    assert_array_equal(S2([1, 2, 3]), [5.0, 7.0, 9.0])


def test_simple_two_model_compose_2d():
    """
    A simple example consisting of two rotations.
    """

    r1 = Rotation2D(45) | Rotation2D(45)

    assert isinstance(r1, CompoundModel)
    assert r1.n_inputs == 2
    assert r1.n_outputs == 2
    assert_allclose(r1(0, 1), (-1, 0), atol=1e-10)

    r2 = Rotation2D(90) | Rotation2D(90)  # Rotate twice by 90 degrees
    assert_allclose(r2(0, 1), (0, -1), atol=1e-10)

    # Compose R with itself to produce 4 rotations
    r3 = r1 | r1

    assert_allclose(r3(0, 1), (0, -1), atol=1e-10)


def test_n_submodels():
    """
    Test that CompoundModel.n_submodels properly returns the number
    of components.
    """

    g2 = Gaussian1D() + Gaussian1D()
    assert g2.n_submodels == 2

    g3 = g2 + Gaussian1D()
    assert g3.n_submodels == 3

    g5 = g3 | g2
    assert g5.n_submodels == 5

    g7 = g5 / g2
    assert g7.n_submodels == 7


def test_expression_formatting():
    """
    Test that the expression strings from compound models are formatted
    correctly.
    """

    # For the purposes of this test it doesn't matter a great deal what
    # model(s) are used in the expression, I don't think
    G = Gaussian1D(1, 1, 1)
    G2 = Gaussian2D(1, 2, 3, 4, 5, 6)

    M = G + G
    assert M._format_expression() == "[0] + [1]"

    M = G + G + G
    assert M._format_expression() == "[0] + [1] + [2]"

    M = G + G * G
    assert M._format_expression() == "[0] + [1] * [2]"

    M = G * G + G
    assert M._format_expression() == "[0] * [1] + [2]"

    M = G + G * G + G
    assert M._format_expression() == "[0] + [1] * [2] + [3]"

    M = (G + G) * (G + G)
    assert M._format_expression() == "([0] + [1]) * ([2] + [3])"

    # This example uses parentheses in the expression, but those won't be
    # preserved in the expression formatting since they technically aren't
    # necessary, and there's no way to know that they were originally
    # parenthesized (short of some deep, and probably not worthwhile
    # introspection)
    M = (G * G) + (G * G)
    assert M._format_expression() == "[0] * [1] + [2] * [3]"

    M = G**G
    assert M._format_expression() == "[0] ** [1]"

    M = G + G**G
    assert M._format_expression() == "[0] + [1] ** [2]"

    M = (G + G) ** G
    assert M._format_expression() == "([0] + [1]) ** [2]"

    M = G + G | G
    assert M._format_expression() == "[0] + [1] | [2]"

    M = G + (G | G)
    assert M._format_expression() == "[0] + ([1] | [2])"

    M = G & G | G2
    assert M._format_expression() == "[0] & [1] | [2]"

    M = G & (G | G)
    assert M._format_expression() == "[0] & ([1] | [2])"


def test_basic_compound_inverse():
    """
    Test basic inversion of compound models in the limited sense supported for
    models made from compositions and joins only.
    """

    t = (Shift(2) & Shift(3)) | (Scale(2) & Scale(3)) | Rotation2D(90)
    assert_allclose(t.inverse(*t(0, 1)), (0, 1))


@pytest.mark.parametrize(
    "model",
    [
        Shift(0) + Shift(0) | Shift(0),
        Shift(0) - Shift(0) | Shift(0),
        Shift(0) * Shift(0) | Shift(0),
        Shift(0) / Shift(0) | Shift(0),
        Shift(0) ** Shift(0) | Shift(0),
        Gaussian1D(1, 2, 3) | Gaussian1D(4, 5, 6),
    ],
)
def test_compound_unsupported_inverse(model):
    """
    Ensure inverses aren't supported in cases where it shouldn't be.
    """

    MESSAGE = r"No analytical or user-supplied inverse transform .*"
    with pytest.raises(NotImplementedError, match=MESSAGE):
        model.inverse


def test_mapping_basic_permutations():
    """
    Tests a couple basic examples of the Mapping model--specifically examples
    that merely permute the outputs.
    """

    x, y = Rotation2D(90)(1, 2)

    rs = Rotation2D(90) | Mapping((1, 0))
    x_prime, y_prime = rs(1, 2)
    assert_allclose((x, y), (y_prime, x_prime))

    # A more complicated permutation
    m = Rotation2D(90) & Scale(2)
    x, y, z = m(1, 2, 3)

    ms = m | Mapping((2, 0, 1))
    x_prime, y_prime, z_prime = ms(1, 2, 3)
    assert_allclose((x, y, z), (y_prime, z_prime, x_prime))


def test_mapping_inverse():
    """Tests inverting a compound model that includes a `Mapping`."""

    rs1 = Rotation2D(12.1) & Scale(13.2)
    rs2 = Rotation2D(14.3) & Scale(15.4)

    # Rotates 2 of the coordinates and scales the third--then rotates on a
    # different axis and scales on the axis of rotation.  No physical meaning
    # here just a simple test
    m = rs1 | Mapping([2, 0, 1]) | rs2

    assert_allclose((0, 1, 2), m.inverse(*m(0, 1, 2)), atol=1e-08)


def test_identity_input():
    """
    Test a case where an Identity (or Mapping) model is the first in a chain
    of composite models and thus is responsible for handling input broadcasting
    properly.

    Regression test for https://github.com/astropy/astropy/pull/3362
    """

    ident1 = Identity(1)
    shift = Shift(1)
    rotation = Rotation2D(angle=90)
    model = ident1 & shift | rotation
    assert_allclose(model(1, 2), [-3.0, 1.0])


def test_invalid_operands():
    """
    Test that certain operators do not work with models whose inputs/outputs do
    not match up correctly.
    """

    MESSAGE = r"Unsupported operands for |:.*"
    with pytest.raises(ModelDefinitionError, match=MESSAGE):
        Rotation2D(90) | Gaussian1D(1, 0, 0.1)

    MESSAGE = r"Both operands must match numbers of inputs and outputs"
    with pytest.raises(ModelDefinitionError, match=MESSAGE):
        Rotation2D(90) + Gaussian1D(1, 0, 0.1)


@pytest.mark.parametrize("poly", [Chebyshev2D(1, 2), Polynomial2D(2), Legendre2D(1, 2)])
def test_compound_with_polynomials_2d(poly):
    """
    Tests that polynomials are scaled when used in compound models.
    Issue #3699
    """
    poly.parameters = [1, 2, 3, 4, 1, 2]
    shift = Shift(3)
    model = poly | shift
    x, y = np.mgrid[:20, :37]
    result_compound = model(x, y)
    result = shift(poly(x, y))
    assert_allclose(result, result_compound)


def test_fix_inputs():
    g1 = Gaussian2D(1, 0, 0, 1, 2)
    g2 = Gaussian2D(1.5, 0.5, -0.2, 0.5, 0.3)
    sg1_1 = fix_inputs(g1, {1: 0})
    assert_allclose(sg1_1(0), g1(0, 0))
    assert_allclose(sg1_1([0, 1, 3]), g1([0, 1, 3], [0, 0, 0]))
    sg1_2 = fix_inputs(g1, {"x": 1})
    assert_allclose(sg1_2(1.5), g1(1, 1.5))
    gg1 = g1 & g2
    sgg1_1 = fix_inputs(gg1, {1: 0.1, 3: 0.2})
    assert_allclose(sgg1_1(0, 0), gg1(0, 0.1, 0, 0.2))
    sgg1_2 = fix_inputs(gg1, {"x0": -0.1, 2: 0.1})
    assert_allclose(sgg1_2(1, 1), gg1(-0.1, 1, 0.1, 1))
    assert_allclose(sgg1_2(y0=1, y1=1), gg1(-0.1, 1, 0.1, 1))


def test_fix_inputs_invalid():
    g1 = Gaussian2D(1, 0, 0, 1, 2)

    MESSAGE = r"Substitution key .* not among possible input choices"
    with pytest.raises(ValueError, match=MESSAGE):
        fix_inputs(g1, {"x0": 0, 0: 0})

    with pytest.raises(ValueError, match=MESSAGE):
        fix_inputs(g1, {3: 2})

    with pytest.raises(ValueError, match=MESSAGE):
        fix_inputs(g1, {np.int32(3): 2})

    with pytest.raises(ValueError, match=MESSAGE):
        fix_inputs(g1, {np.int64(3): 2})

    with pytest.raises(ValueError, match=MESSAGE):
        fix_inputs(g1, {"w": 2})

    MESSAGE = r'Expected a dictionary for second argument of "fix_inputs"'
    with pytest.raises(ValueError, match=MESSAGE):
        fix_inputs(g1, (0, 1))

    MESSAGE = r".*Illegal operator: ', '#'.*"
    with pytest.raises(ModelDefinitionError, match=MESSAGE):
        CompoundModel("#", g1, g1)

    MESSAGE = r"Too many input arguments - expected 1, got 2"
    with pytest.raises(ValueError, match=MESSAGE):
        gg1 = fix_inputs(g1, {0: 1})
        gg1(2, y=2)

    with pytest.raises(ValueError, match=MESSAGE):
        gg1 = fix_inputs(g1, {np.int32(0): 1})
        gg1(2, y=2)

    with pytest.raises(ValueError, match=MESSAGE):
        gg1 = fix_inputs(g1, {np.int64(0): 1})
        gg1(2, y=2)


def test_fix_inputs_with_bounding_box():
    g1 = Gaussian2D(1, 0, 0, 1, 1)
    g2 = Gaussian2D(1, 0, 0, 1, 1)
    assert_allclose(g1.bounding_box, ((-5.5, 5.5), (-5.5, 5.5)))

    gg1 = g1 & g2
    gg1.bounding_box = ((-5.5, 5.5), (-5.4, 5.4), (-5.3, 5.3), (-5.2, 5.2))
    assert gg1.bounding_box == ((-5.5, 5.5), (-5.4, 5.4), (-5.3, 5.3), (-5.2, 5.2))

    sg = fix_inputs(gg1, {0: 0, 2: 0})
    assert sg.bounding_box == ((-5.5, 5.5), (-5.3, 5.3))

    g1 = Gaussian1D(10, 3, 1)
    g = g1 & g1
    g.bounding_box = ((1, 4), (6, 8))
    gf = fix_inputs(g, {0: 1})
    assert gf.bounding_box == (1, 4)


def test_indexing_on_instance():
    """Test indexing on compound model instances."""

    m = Gaussian1D(1, 0, 0.1) + Const1D(2)
    assert isinstance(m[0], Gaussian1D)
    assert isinstance(m[1], Const1D)
    assert m.param_names == ("amplitude_0", "mean_0", "stddev_0", "amplitude_1")

    # Test parameter equivalence
    assert m[0].amplitude == 1 == m.amplitude_0
    assert m[0].mean == 0 == m.mean_0
    assert m[0].stddev == 0.1 == m.stddev_0
    assert m[1].amplitude == 2 == m.amplitude_1

    # Test that parameter value updates are symmetric between the compound
    # model and the submodel returned by indexing
    const = m[1]
    m.amplitude_1 = 42
    assert const.amplitude == 42
    const.amplitude = 137
    assert m.amplitude_1 == 137

    # Similar couple of tests, but now where the compound model was created
    # from model instances
    g = Gaussian1D(1, 2, 3, name="g")
    p = Polynomial1D(2, name="p")
    m = g + p
    assert m[0].name == "g"
    assert m[1].name == "p"
    assert m["g"].name == "g"
    assert m["p"].name == "p"

    poly = m[1]
    m.c0_1 = 12345
    assert poly.c0 == 12345
    poly.c1 = 6789
    assert m.c1_1 == 6789

    # Test negative indexing
    assert isinstance(m[-1], Polynomial1D)
    assert isinstance(m[-2], Gaussian1D)

    MESSAGE = r"list index out of range"
    with pytest.raises(IndexError, match=MESSAGE):
        m[42]

    MESSAGE = r"No component with name 'foobar' found"
    with pytest.raises(IndexError, match=MESSAGE):
        m["foobar"]

    # Confirm index-by-name works with fix_inputs
    g = Gaussian2D(1, 2, 3, 4, 5, name="g")
    m = fix_inputs(g, {0: 1})
    assert m["g"].name == "g"

    # Test string slicing
    A = Const1D(1.1, name="A")
    B = Const1D(2.1, name="B")
    C = Const1D(3.1, name="C")
    M = A + B * C
    assert_allclose(M["B":"C"](1), 6.510000000000001)


class _ConstraintsTestA(Model):
    stddev = Parameter(default=0, min=0, max=0.3)
    mean = Parameter(default=0, fixed=True)

    @staticmethod
    def evaluate(stddev, mean):
        return stddev, mean


class _ConstraintsTestB(Model):
    mean = Parameter(default=0, fixed=True)

    @staticmethod
    def evaluate(mean):
        return mean


def test_inherit_constraints():
    """
    Various tests for copying of constraint values between compound models and
    their members.

    Regression test for https://github.com/astropy/astropy/issues/3481
    """
    model = Gaussian1D(bounds={"stddev": (0, 0.3)}, fixed={"mean": True}) + Gaussian1D(
        fixed={"mean": True}
    )

    # Lots of assertions in this test as there are multiple interfaces to
    # parameter constraints

    assert "stddev_0" in model.bounds
    assert model.bounds["stddev_0"] == (0, 0.3)
    assert model.stddev_0.bounds == (0, 0.3)
    assert "mean_0" in model.fixed
    assert model.fixed["mean_0"] is True
    assert model.mean_0.fixed is True
    assert "mean_1" in model.fixed
    assert model.fixed["mean_1"] is True
    assert model.mean_1.fixed is True

    assert model.stddev_0 is model[0].stddev
    # Great, all the constraints were inherited properly
    # Now what about if we update them through the sub-models?
    model.stddev_0.bounds = (0, 0.4)
    assert model[0].stddev.bounds == (0, 0.4)
    assert model[0].bounds["stddev"] == (0, 0.4)

    model.stddev_0.bounds = (0.1, 0.5)
    assert model[0].stddev.bounds == (0.1, 0.5)
    assert model[0].bounds["stddev"] == (0.1, 0.5)

    model[1].mean.fixed = False
    assert model.mean_1.fixed is False
    assert model[1].mean.fixed is False

    # Now turn off syncing of constraints
    assert model.bounds["stddev_0"] == (0.1, 0.5)
    model.sync_constraints = False
    model[0].stddev.bounds = (0, 0.2)
    assert model.bounds["stddev_0"] == (0.1, 0.5)
    model.sync_constraints = True
    assert model.bounds["stddev_0"] == (0, 0.2)


def test_compound_custom_inverse():
    """
    Test that a compound model with a custom inverse has that inverse applied
    when the inverse of another model, of which it is a component, is computed.
    Regression test for https://github.com/astropy/astropy/issues/3542
    """

    poly = Polynomial1D(1, c0=1, c1=2)
    scale = Scale(1)
    shift = Shift(1)

    model1 = poly | scale
    model1.inverse = poly

    # model1 now has a custom inverse (the polynomial itself, ignoring the
    # trivial scale factor)
    model2 = shift | model1

    assert_allclose(model2.inverse(1), (poly | shift.inverse)(1))

    # Make sure an inverse is not allowed if the models were combined with the
    # wrong operator, or if one of the models doesn't have an inverse defined
    MESSAGE = (
        r"No analytical or user-supplied inverse transform has been implemented for"
        r" this model"
    )
    with pytest.raises(NotImplementedError, match=MESSAGE):
        (shift + model1).inverse

    with pytest.raises(NotImplementedError, match=MESSAGE):
        (model1 & poly).inverse


def test_pickle_compound():
    """
    Regression test for
    https://github.com/astropy/astropy/issues/3867#issuecomment-114547228
    """

    # Test pickling a compound model instance
    g1 = Gaussian1D(1.0, 0.0, 0.1)
    g2 = Gaussian1D([2.0, 3.0], [0.0, 0.0], [0.2, 0.3])
    m = g1 + g2
    m2 = pickle.loads(pickle.dumps(m))
    assert m.param_names == m2.param_names
    assert m.__class__.__name__ == m2.__class__.__name__
    assert np.all(m.parameters == m2.parameters)
    assert np.all(m(0) == m2(0))


def test_update_parameters():
    offx = Shift(1)
    scl = Scale(2)
    m = offx | scl
    assert m(1) == 4

    offx.offset = 42
    assert m(1) == 86

    m.factor_1 = 100
    assert m(1) == 4300
    m2 = m | offx
    assert m2(1) == 4342


def test_name():
    offx = Shift(1)
    scl = Scale(2)
    m = offx | scl
    scl.name = "scale"
    assert m.submodel_names == ("None_0", "scale")
    assert m.name is None
    m.name = "M"
    assert m.name == "M"
    m1 = m.rename("M1")
    assert m.name == "M1"
    assert m1.name == "M1"


def test_name_index():
    g1 = Gaussian1D(1, 1, 1)
    g2 = Gaussian1D(1, 2, 1)
    g = g1 + g2

    MESSAGE = r"No component with name 'bozo' found"
    with pytest.raises(IndexError, match=MESSAGE):
        g["bozo"]
    g1.name = "bozo"
    assert g["bozo"].mean == 1
    g2.name = "bozo"
    MESSAGE = r"Multiple components found using 'bozo' as name.*"
    with pytest.raises(IndexError, match=MESSAGE):
        g["bozo"]


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_tabular_in_compound():
    """
    Issue #7411 - evaluate should not change the shape of the output.
    """
    t = Tabular1D(points=([1, 5, 7],), lookup_table=[12, 15, 19], bounds_error=False)
    rot = Rotation2D(2)
    p = Polynomial1D(1)
    x = np.arange(12).reshape((3, 4))
    # Create a compound model which does not execute Tabular.__call__,
    # but model.evaluate and is followed by a Rotation2D which
    # checks the exact shapes.
    model = p & t | rot
    x1, y1 = model(x, x)
    assert x1.ndim == 2
    assert y1.ndim == 2


def test_bounding_box():
    g = Gaussian2D() + Gaussian2D(2, 0.5, 0.1, 2, 3, 0)
    g.bounding_box = ((0, 1), (0, 0.5))
    y, x = np.mgrid[0:10, 0:10]
    y = y / 3.0
    x = x / 3.0
    val = g(x, y, with_bounding_box=True)

    # fmt: off
    compare = np.array(
        [
            [2.93738984, 2.93792011, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [2.87857153, 2.88188761, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [2.70492922, 2.71529265, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [2.45969972, 2.47912103, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan, np.nan, np.nan]
        ]
    )
    # fmt: on

    mask = ~np.isnan(val)
    assert_allclose(val[mask], compare[mask])
    val2 = g(x + 2, y + 2, with_bounding_box=True)
    assert np.isnan(val2).sum() == 100
    # val3 = g(.1, .1, with_bounding_box=True)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_bounding_box_with_units():
    points = np.arange(5) * u.pix
    lt = np.arange(5) * u.AA
    t = Tabular1D(points, lt)

    assert t(1 * u.pix, with_bounding_box=True) == 1.0 * u.AA


@pytest.mark.parametrize("poly", [Chebyshev1D(5), Legendre1D(5), Polynomial1D(5)])
def test_compound_with_polynomials_1d(poly):
    """
    Tests that polynomials are offset when used in compound models.
    Issue #3699
    """
    poly.parameters = [1, 2, 3, 4, 1, 2]
    shift = Shift(3)
    model = poly | shift
    x = np.linspace(-5, 5, 10)
    result_compound = model(x)
    result = shift(poly(x))
    assert_allclose(result, result_compound)
    assert model.param_names == (
        "c0_0",
        "c1_0",
        "c2_0",
        "c3_0",
        "c4_0",
        "c5_0",
        "offset_1",
    )


def test_replace_submodel():
    """
    Replace a model in a Compound model
    """
    S1 = Shift(2, name="shift2") | Scale(3, name="scale3")  # First shift then scale
    S2 = Scale(2, name="scale2") | Shift(3, name="shift3")  # First scale then shift

    m = S1 & S2
    assert m(1, 2) == (9, 7)

    m2 = m.replace_submodel("scale3", Scale(4, name="scale4"))
    assert m2(1, 2) == (12, 7)
    assert m(1, 2) == (9, 7)
    # Check the inverse has been updated
    assert m2.inverse(12, 7) == (1, 2)

    # Produce the same result by replacing a single model with a compound
    m3 = m.replace_submodel("shift2", Shift(2) | Scale(2))
    assert m(1, 2) == (9, 7)
    assert m3(1, 2) == (18, 7)
    # Check the inverse has been updated
    assert m3.inverse(18, 7) == (1, 2)

    # Test with arithmetic model compunding operator
    m = S1 + S2
    assert m(1) == 14
    m2 = m.replace_submodel("scale2", Scale(4, name="scale4"))
    assert m2(1) == 16

    # Test with fix_inputs()
    R = fix_inputs(Rotation2D(angle=90, name="rotate"), {0: 1})
    m4 = S1 | R
    assert_allclose(m4(0), (-6, 1))

    m5 = m4.replace_submodel("rotate", Rotation2D(180))
    assert_allclose(m5(0), (-1, -6))

    # Check we get a value error when model name doesn't exist
    MESSAGE = r"No submodels found named not_there"
    with pytest.raises(ValueError, match=MESSAGE):
        m2 = m.replace_submodel("not_there", Scale(2))

    # And now a model set
    P = Polynomial1D(degree=1, n_models=2, name="poly")
    S = Shift([1, 2], n_models=2)
    m = P | S
    assert_array_equal(m([0, 1]), (1, 2))
    MESSAGE = r"New and old models must have equal values for n_models"
    with pytest.raises(ValueError, match=MESSAGE):
        m2 = m.replace_submodel("poly", Polynomial1D(degree=1, c0=1))
    m2 = m.replace_submodel("poly", Polynomial1D(degree=1, c0=[1, 2], n_models=2))
    assert_array_equal(m2([0, 1]), (2, 4))

    # Ensure previous _user_inverse doesn't stick around
    S1 = Shift(1)
    S2 = Shift(2)
    S3 = Shift(3, name="S3")

    S23 = S2 | S3
    S23.inverse = Shift(-4.9)
    m = S1 & S23

    # This should delete the S23._user_inverse
    m2 = m.replace_submodel("S3", Shift(4))
    assert m2(1, 2) == (2, 8)
    assert m2.inverse(2, 8) == (1, 2)


@pytest.mark.parametrize(
    "expr",
    [
        lambda m1, m2: m1 + m2,
        lambda m1, m2: m1 - m2,
        lambda m1, m2: m1 * m2,
        lambda m1, m2: m1 / m2,
    ],
)
def test_compound_evaluate(expr):
    """
    Tests that compound evaluate function produces the same
    result as the models with the operator applied
    """
    x = np.linspace(-5, 5, 10)
    # Some evaluate functions assume that inputs are numpy arrays or quantities including Const1D
    p1 = np.array([1, 2, 3, 4, 1, 2])
    p2 = np.array([1, 0, 0.5])

    model1 = Polynomial1D(5)
    model2 = Gaussian1D(2, 1, 5)
    compound = expr(model1, model2)

    assert_array_equal(
        compound.evaluate(x, *p1, *p2),
        expr(model1.evaluate(x, *p1), model2.evaluate(x, *p2)),
    )


def test_compound_evaluate_power():
    """
    Tests that compound evaluate function produces the same
    result as the models with the power operator applied
    """
    x = np.linspace(-5, 5, 10)
    p1 = np.array([1, 0, 0.2])
    p2 = np.array([3])

    model1 = Gaussian1D(2, 1, 5)
    model2 = Const1D(2)
    compound = model1**model2

    assert_array_equal(
        compound.evaluate(x, *p1, *p2),
        model1.evaluate(x, *p1) ** model2.evaluate(x, *p2),
    )


def test_compound_evaluate_double_shift():
    x = np.linspace(-5, 5, 10)
    y = np.linspace(-5, 5, 10)

    m1 = Gaussian2D(1, 0, 0, 1, 1, 1)
    m2 = Shift(1)
    m3 = Shift(2)
    m = Gaussian2D(1, 0, 0, 1, 1, 1) & Shift(1) & Shift(2)
    assert_array_equal(
        m.evaluate(x, y, x - 10, y + 20, 1, 0, 0, 1, 1, 1, 1, 2),
        [
            m1.evaluate(x, y, 1, 0, 0, 1, 1, 1),
            m2.evaluate(x - 10, 1),
            m3.evaluate(y + 20, 2),
        ],
    )


@pytest.mark.parametrize(
    "expr",
    [
        lambda m1, m2: m1 + m2,
        lambda m1, m2: m1 - m2,
        lambda m1, m2: m1 * m2,
        lambda m1, m2: m1 / m2,
    ],
)
def test_compound_evaluate_named_param(expr):
    """
    Tests that compound evaluate function produces the same
    result as the models with the operator applied
    """
    x = np.linspace(-5, 5, 10)
    p1 = np.array([1, 0, 0.2])
    p2 = np.array([3, 0.5, 0.5])

    model1 = Gaussian1D(2, 1, 5)
    model2 = Gaussian1D(2, 1, 5)
    compound = expr(model1, model2)

    assert_array_equal(
        compound.evaluate(x, *p2, amplitude_0=p1[0], mean_0=p1[1], stddev_0=p1[2]),
        expr(model1.evaluate(x, *p1), model2.evaluate(x, *p2)),
    )


def test_compound_evaluate_name_param_power():
    """
    Tests that compound evaluate function produces the same
    result as the models with the power operator applied
    """
    x = np.linspace(-5, 5, 10)
    p1 = np.array([1, 0, 0.2])
    p2 = np.array([3])

    model1 = Gaussian1D(2, 1, 5)
    model2 = Const1D(2)
    compound = model1**model2

    assert_array_equal(
        compound.evaluate(x, *p2, amplitude_0=p1[0], mean_0=p1[1], stddev_0=p1[2]),
        model1.evaluate(x, *p1) ** model2.evaluate(x, *p2),
    )


def test_compound_evaluate_and():
    """
    Tests that compound evaluate function produces the same
    result as the models with the operator applied
    """
    x = np.linspace(-5, 5, 10)
    p1 = np.array([1, 0.1, 0.5])
    p2 = np.array([3])

    model1 = Gaussian1D()
    model2 = Shift()
    compound = model1 & model2

    assert_array_equal(
        compound.evaluate(x, x, *p1, p2),
        [model1.evaluate(x, *p1), model2.evaluate(x, p2)],
    )


def test_compound_evaluate_or():
    """
    Tests that compound evaluate function produces the same
    result as the models with the operator applied
    """
    x = np.linspace(-5, 5, 10)
    p1 = np.array([0.5])
    p2_amplitude = np.array([3])
    p2_mean = np.array([0])
    p2_std = np.array([0.1])

    model1 = Shift(0.5)
    model2 = Gaussian1D(1, 0, 0.5)
    compound = model1 | model2

    assert_array_equal(
        compound.evaluate(x, p1, p2_amplitude, p2_mean, p2_std),
        model2.evaluate(model1.evaluate(x, p1), p2_amplitude, p2_mean, p2_std),
    )


def test_compound_evaluate_fix_inputs_by_keyword():
    """
    Tests that compound evaluate function produces the same
    result as the models fix_inputs operator is applied
    when using the keyword
    """
    y, x = np.mgrid[:10, :10]

    model_params = [3, 0, 0.1, 1, 0.5, 0]

    model = Gaussian2D(1, 2, 0, 0.5)
    compound = fix_inputs(model, {"x": x + 5})

    assert_array_equal(
        compound.evaluate(x, y, *model_params),
        model.evaluate(x + 5, y, *model_params),
    )


def test_compound_evaluate_fix_inputs_by_position():
    """
    Tests that compound evaluate function produces the same
    result as the models fix_inputs operator is applied
    when using the input index
    """
    y, x = np.mgrid[:10, :10]

    model_params = [3, 0, 0.1, 1, 0.5, 0]

    model = Gaussian2D(1, 2, 0, 0.5)
    compound = fix_inputs(model, {0: x + 5})

    assert_array_equal(
        compound.evaluate(x, y, *model_params),
        model.evaluate(x + 5, y, *model_params),
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fit_multiplied_compound_model_with_mixed_units():
    """
    Regression test for issue #12320
    """

    fitter = LevMarLSQFitter()
    x = np.linspace(0, 1, 101) * u.s
    y = np.linspace(5, 10, 101) * u.m * u.kg / u.s

    m1 = Linear1D(slope=5 * u.m / u.s / u.s, intercept=1.0 * u.m / u.s)
    m2 = Linear1D(slope=0.0 * u.kg / u.s, intercept=10.0 * u.kg)
    truth = m1 * m2

    # We need to fix some of the parameters to avoid degeneracies
    truth.slope_1.fixed = True
    truth.intercept_0.fixed = True

    fit = fitter(truth, x, y)

    unfit_output = truth(x)
    fit_output = fit(x)

    assert unfit_output.unit == fit_output.unit == (u.kg * u.m / u.s)

    # The unfit model is 10 kg m^2 / s for x=0s and goes up to 60 kg m^2 / s
    # for x=1s, whereas the actual data being fit goes from 5 to 10, so we need
    # to correct this.
    assert_allclose(unfit_output / 10 + 4 * u.kg * u.m / u.s, fit_output)

    assert_quantity_allclose(fit.slope_0, 1 * u.m / u.s / u.s)
    assert_quantity_allclose(fit.intercept_0, 1.0 * u.m / u.s)
    assert_quantity_allclose(fit.slope_1, 0 * u.kg / u.s)
    assert_quantity_allclose(fit.intercept_1, 5 * u.kg)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fit_multiplied_recursive_compound_model_with_mixed_units():
    """
    Regression test for issue #12320
    """

    fitter = LevMarLSQFitter()

    x = np.linspace(0, 1, 101) * u.s
    y = np.linspace(5, 10, 101) * u.m * u.m * u.kg / u.s

    m1 = Linear1D(slope=5 * u.m / u.s / u.s, intercept=1.0 * u.m / u.s)
    m2 = Linear1D(slope=0.0 * u.kg / u.s, intercept=10.0 * u.kg)
    m3 = Linear1D(slope=0.0 * u.m / u.s, intercept=10.0 * u.m)
    truth = m1 * m2 * m3

    # We need to fix some of the parameters to avoid degeneracies
    truth.slope_1.fixed = True
    truth.slope_2.fixed = True
    truth.intercept_0.fixed = True
    truth.intercept_1.fixed = True

    fit = fitter(truth, x, y)

    unfit_output = truth(x)
    fit_output = fit(x)

    assert unfit_output.unit == fit_output.unit == (u.kg * u.m * u.m / u.s)

    # The unfit model is 100 kg m^2 / s for x=0s and goes up to 600 kg m^2 / s
    # for x=1s, whereas the actual data being fit goes from 5 to 10, so we need
    # to correct this.
    assert_allclose(unfit_output / 100 + 4 * u.kg * u.m * u.m / u.s, fit_output)

    assert_quantity_allclose(fit.slope_0, 1 * u.m / u.s / u.s)
    assert_quantity_allclose(fit.intercept_0, 1.0 * u.m / u.s)
    assert_quantity_allclose(fit.slope_1, 0 * u.kg / u.s)
    assert_quantity_allclose(fit.intercept_1, 10 * u.kg)
    assert_quantity_allclose(fit.slope_2, 0 * u.m / u.s)
    assert_quantity_allclose(fit.intercept_2, 0.5 * u.m)

    x = np.linspace(0, 1, 101) * u.s
    y = np.linspace(5, 10, 101) * u.m * u.m * u.kg * u.kg / u.s

    m1 = Linear1D(slope=5 * u.m / u.s / u.s, intercept=1.0 * u.m / u.s)
    m2 = Linear1D(slope=0.0 * u.kg / u.s, intercept=10.0 * u.kg)
    m3 = Linear1D(slope=0.0 * u.m / u.s, intercept=10.0 * u.m)
    m4 = Linear1D(slope=0.0 * u.kg / u.s, intercept=10.0 * u.kg)
    m11 = m1 * m2
    m22 = m3 * m4
    truth = m11 * m22

    # We need to fix some of the parameters to avoid degeneracies
    truth.slope_1.fixed = True
    truth.slope_2.fixed = True
    truth.slope_3.fixed = True
    truth.intercept_0.fixed = True
    truth.intercept_1.fixed = True
    truth.intercept_2.fixed = True

    fit = fitter(truth, x, y)

    unfit_output = truth(x)
    fit_output = fit(x)

    assert unfit_output.unit == fit_output.unit == (u.kg * u.kg * u.m * u.m / u.s)

    # The unfit model is 1000 kg m^2 / s for x=0s and goes up to 6000 kg m^2 / s
    # for x=1s, whereas the actual data being fit goes from 5 to 10, so we need
    # to correct this.
    assert_allclose(unfit_output / 1000 + 4 * u.kg * u.kg * u.m * u.m / u.s, fit_output)

    assert_quantity_allclose(fit.slope_0, 1 * u.m / u.s / u.s)
    assert_quantity_allclose(fit.intercept_0, 1.0 * u.m / u.s)
    assert_quantity_allclose(fit.slope_1, 0 * u.kg / u.s)
    assert_quantity_allclose(fit.intercept_1, 10 * u.kg)
    assert_quantity_allclose(fit.slope_2, 0 * u.m / u.s)
    assert_quantity_allclose(fit.intercept_2, 10 * u.m)
    assert_quantity_allclose(fit.slope_3, 0 * u.kg / u.s)
    assert_quantity_allclose(fit.intercept_3, 0.05 * u.kg)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fit_divided_compound_model_with_mixed_units():
    """
    Regression test for issue #12320
    """

    fitter = LevMarLSQFitter()
    x = np.linspace(0, 1, 101) * u.s
    y = np.linspace(5, 10, 101) * u.kg * u.m / u.s

    m1 = Linear1D(slope=5 * u.kg * u.m / u.s, intercept=1.0 * u.kg * u.m)
    m2 = Linear1D(slope=0.0 * u.s / u.s, intercept=10.0 * u.s)
    truth = m1 / m2

    # We need to fix some of the parameters to avoid degeneracies
    truth.slope_1.fixed = True
    truth.intercept_0.fixed = True

    fit = fitter(truth, x, y)

    unfit_output = truth(x)
    fit_output = fit(x)

    assert unfit_output.unit == fit_output.unit == (u.kg * u.m / u.s)

    # The unfit model is 0.1 kg m / s for x=0s and goes up to 0.5 kg m / s for
    # x=1s, whereas the actual data being fit goes from 5 to 10, so we need
    # to correct this.
    assert_allclose(unfit_output * 10 + 4 * u.kg * u.m / u.s, fit_output, rtol=1e-4)

    assert_quantity_allclose(fit.slope_0, 1 * u.kg * u.m / u.s)
    assert_quantity_allclose(fit.intercept_0, 1.0 * u.kg * u.m)
    assert_quantity_allclose(fit.slope_1, 0 * u.s / u.s)
    assert_quantity_allclose(fit.intercept_1, 0.2 * u.s)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fit_mixed_recursive_compound_model_with_mixed_units():
    """
    Regression test for issue #12320
    """

    fitter = LevMarLSQFitter()

    x = np.linspace(0, 1, 101) * u.s
    y = np.linspace(5, 10, 101) * u.kg * u.m * u.m / u.s

    m1 = Linear1D(slope=5 * u.kg * u.m / u.s, intercept=1.0 * u.kg * u.m)
    m2 = Linear1D(slope=0.0 * u.s / u.s, intercept=10.0 * u.s)
    m3 = Linear1D(slope=0.0 * u.m / u.s, intercept=10.0 * u.m)
    truth = m1 / m2 * m3

    # We need to fix some of the parameters to avoid degeneracies
    truth.slope_1.fixed = True
    truth.slope_2.fixed = True
    truth.intercept_0.fixed = True
    truth.intercept_1.fixed = True

    fit = fitter(truth, x, y)

    unfit_output = truth(x)
    fit_output = fit(x)

    assert unfit_output.unit == fit_output.unit == (u.kg * u.m * u.m / u.s)

    # The unfit model is 1 kg m^2 / s for x=0s and goes up to 6 kg m^2 / s for
    # x=1s, whereas the actual data being fit goes from 5 to 10, so we need
    # to correct this.
    assert_allclose(unfit_output + 4 * u.kg * u.m * u.m / u.s, fit_output)

    assert_quantity_allclose(fit.slope_0, 1 * u.kg * u.m / u.s)
    assert_quantity_allclose(fit.intercept_0, 1.0 * u.kg * u.m)
    assert_quantity_allclose(fit.slope_1, 0 * u.s / u.s)
    assert_quantity_allclose(fit.intercept_1, 10 * u.s)
    assert_quantity_allclose(fit.slope_2, 0 * u.m / u.s)
    assert_quantity_allclose(fit.intercept_2, 50 * u.m)

    x = np.linspace(0, 1, 101) * u.s
    y = np.linspace(5, 10, 101) * u.kg * u.kg * u.m * u.m / u.s

    m1 = Linear1D(slope=5 * u.kg * u.m / u.s, intercept=1.0 * u.kg * u.m)
    m2 = Linear1D(slope=0.0 * u.s / u.s, intercept=10.0 * u.s)
    m3 = Linear1D(slope=0.0 * u.m / u.s, intercept=10.0 * u.m)
    m4 = Linear1D(slope=0.0 * u.kg / u.s, intercept=10.0 * u.kg)
    m11 = m1 / m2
    m22 = m3 * m4
    truth = m11 * m22

    # We need to fix some of the parameters to avoid degeneracies
    truth.slope_1.fixed = True
    truth.slope_2.fixed = True
    truth.slope_3.fixed = True
    truth.intercept_0.fixed = True
    truth.intercept_1.fixed = True
    truth.intercept_2.fixed = True

    fit = fitter(truth, x, y)

    unfit_output = truth(x)
    fit_output = fit(x)

    assert unfit_output.unit == fit_output.unit == (u.kg * u.kg * u.m * u.m / u.s)

    # The unfit model is 10 kg^2 m^2 / s for x=0s and goes up to 60 kg^2 m^2 / s
    # for x=1s, whereas the actual data being fit goes from 5 to 10, so we need
    # to correct this.

    assert_allclose(unfit_output / 10 + 4 * u.kg * u.kg * u.m * u.m / u.s, fit_output)

    assert_quantity_allclose(fit.slope_0, 1 * u.kg * u.m / u.s)
    assert_quantity_allclose(fit.intercept_0, 1.0 * u.kg * u.m)
    assert_quantity_allclose(fit.slope_1, 0 * u.s / u.s)
    assert_quantity_allclose(fit.intercept_1, 10 * u.s)
    assert_quantity_allclose(fit.slope_2, 0 * u.m / u.s)
    assert_quantity_allclose(fit.intercept_2, 10 * u.m)
    assert_quantity_allclose(fit.slope_3, 0 * u.kg / u.s)
    assert_quantity_allclose(fit.intercept_3, 5 * u.kg)


def numerical_partial_deriv(model, *inputs, param_idx, delta=1e-5):
    """
    Evaluate the central difference approximation of the derivative for param_idx.

    Parameters
    ----------
    model
        The model to evaluate
    inputs
        The inputs to the model
    param_idx
        The index of the parameter to compute the partial derivative for.
    delta
        The step size with which to compute the central difference.
    """
    param = model.parameters

    param_down = param.copy()
    param_down[param_idx] = param[param_idx] - delta
    param_up = param.copy()
    param_up[param_idx] = param[param_idx] + delta

    up = model.evaluate(*inputs, *param_up)
    down = model.evaluate(*inputs, *param_down)

    return (up - down) / (2 * delta)


@pytest.mark.parametrize(
    "model",
    [
        pytest.param(
            m, id=m._format_expression(format_leaf=lambda i, l: type(l).__name__)
        )
        for m in [
            Gaussian1D(5, 2, 3) + Linear1D(2, 3),
            Gaussian1D(5, 2, 3) - Linear1D(2, 3),
            Polynomial1D(2) * Gaussian1D(),
            Polynomial1D(2) / Gaussian1D(),
            Polynomial1D(2) + Gaussian1D(),
            Polynomial2D(2) + Gaussian2D(),
        ]
    ],
)
@pytest.mark.parametrize("input_ndim", (1, 2))
def test_compound_fit_deriv(model, input_ndim):
    """
    Given some compound models compare the numerical derivatives to analytical ones.
    """

    x = np.linspace(1, 5, num=10)
    y = np.linspace(1, 5, num=10)

    if input_ndim == 2:
        x = x.reshape((5, 2))
        y = y.reshape((5, 2))

    inputs = (x,) if model.n_inputs == 1 else (x, y)

    numerical = [
        numerical_partial_deriv(model, *inputs, param_idx=i)
        for i in range(len(model.parameters))
    ]
    analytical = model.fit_deriv(*inputs, *model.parameters)

    numerical = np.asarray(numerical)
    analytical = np.asarray(analytical)

    # Reshape output to ravel all but the first dimension since some models do this
    numerical = numerical.reshape((numerical.shape[0], -1))
    analytical = analytical.reshape((analytical.shape[0], -1))

    assert_allclose(numerical, analytical)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fit_compound_polynomial2d():
    """
    Regression test for a bug that caused compound models with Polynomial2D
    to not be fittable due to a bug in CompoundModel.fit_deriv
    """

    # Generate fake data
    rng = np.random.default_rng(0)
    y, x = np.mgrid[:128, :128]
    z = 2.0 * x**2 - 0.5 * x**2 + 1.5 * x * y - 1.0
    z += rng.normal(0.0, 0.1, z.shape) * 50000.0
    z += Gaussian2D(amplitude=50000, x_mean=60, y_mean=60, x_stddev=5, y_stddev=5)(x, y)

    # Fit the data using astropy.modeling
    p_init = Polynomial2D(degree=2) + Gaussian2D(amplitude=50000, x_mean=60, y_mean=60)
    fit_p = DogBoxLSQFitter()

    # We just make sure the fitting works, as it previously crashed
    fit_p(p_init, x, y, z)

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import pickle

from copy import deepcopy

import numpy as np

from numpy.testing.utils import assert_allclose, assert_array_equal

from astropy.utils import minversion
from astropy.modeling.core import Model, ModelDefinitionError
from astropy.modeling.parameters import Parameter
from astropy.modeling.models import (Const1D, Shift, Scale, Rotation2D, Gaussian1D,
                                     Gaussian2D, Polynomial1D, Polynomial2D,
                                     Chebyshev2D, Legendre2D, Chebyshev1D, Legendre1D,
                                     Identity, Mapping,
                                     Tabular1D)
import astropy.units as u
from ..core import CompoundModel


try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

HAS_SCIPY_14 = HAS_SCIPY and minversion(scipy, "0.14")


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, [5.0, 5.0]),
                          (lambda x, y: x - y, [-1.0, -1.0]),
                          (lambda x, y: x * y, [6.0, 6.0]),
                          (lambda x, y: x / y, [2.0 / 3.0, 2.0 / 3.0]),
                          (lambda x, y: x ** y, [8.0, 8.0])])
def test_model_set(expr, result):
    s = expr(Const1D((2, 2), n_models=2), Const1D((3, 3), n_models=2))
    out = s(0, model_set_axis=False)
    assert_array_equal(out, result)


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, [5.0, 5.0]),
                          (lambda x, y: x - y, [-1.0, -1.0]),
                          (lambda x, y: x * y, [6.0, 6.0]),
                          (lambda x, y: x / y, [2.0 / 3.0, 2.0 / 3.0]),
                          (lambda x, y: x ** y, [8.0, 8.0])])
def test_model_set_raises_value_error(expr, result):
    """Check that creating model sets with components whose _n_models are
       different raise a value error
    """
    with pytest.raises(ValueError):
        s = expr(Const1D((2, 2), n_models=2), Const1D(3, n_models=1))


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, 5.0),
                          (lambda x, y: x - y, -1.0),
                          (lambda x, y: x * y, 6.0),
                          (lambda x, y: x / y, 2.0 / 3.0),
                          (lambda x, y: x ** y, 8.0)])
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
    assert M._format_expression() == '[0] + [1]'

    M = G + G + G
    assert M._format_expression() == '[0] + [1] + [2]'

    M = G + G * G
    assert M._format_expression() == '[0] + [1] * [2]'

    M = G * G + G
    assert M._format_expression() == '[0] * [1] + [2]'

    M = G + G * G + G
    assert M._format_expression() == '[0] + [1] * [2] + [3]'

    M = (G + G) * (G + G)
    assert M._format_expression() == '([0] + [1]) * ([2] + [3])'

    # This example uses parentheses in the expression, but those won't be
    # preserved in the expression formatting since they technically aren't
    # necessary, and there's no way to know that they were originally
    # parenthesized (short of some deep, and probably not worthwhile
    # introspection)
    M = (G * G) + (G * G)
    assert M._format_expression() == '[0] * [1] + [2] * [3]'

    M = G ** G
    assert M._format_expression() == '[0] ** [1]'

    M = G + G ** G
    assert M._format_expression() == '[0] + [1] ** [2]'

    M = (G + G) ** G
    assert M._format_expression() == '([0] + [1]) ** [2]'

    M = G + G | G
    assert M._format_expression() == '[0] + [1] | [2]'

    M = G + (G | G)
    assert M._format_expression() == '[0] + ([1] | [2])'

    M = G & G | G2
    assert M._format_expression() == '[0] & [1] | [2]'

    M = G & (G | G)
    assert M._format_expression() == '[0] & ([1] | [2])'


def test_basic_compound_inverse():
    """
    Test basic inversion of compound models in the limited sense supported for
    models made from compositions and joins only.
    """

    t = (Shift(2) & Shift(3)) | (Scale(2) & Scale(3)) | Rotation2D(90)
    assert_allclose(t.inverse(*t(0, 1)), (0, 1))


@pytest.mark.parametrize('model', [
    Shift(0) + Shift(0) | Shift(0),
    Shift(0) - Shift(0) | Shift(0),
    Shift(0) * Shift(0) | Shift(0),
    Shift(0) / Shift(0) | Shift(0),
    Shift(0) ** Shift(0) | Shift(0),
    Gaussian1D(1, 2, 3) | Gaussian1D(4, 5, 6)])
def test_compound_unsupported_inverse(model):
    """
    Ensure inverses aren't supported in cases where it shouldn't be.
    """

    with pytest.raises(NotImplementedError):
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

    with pytest.raises(ModelDefinitionError):
        Rotation2D(90) | Gaussian1D(1, 0, 0.1)

    with pytest.raises(ModelDefinitionError):
        Rotation2D(90) + Gaussian1D(1, 0, 0.1)


@pytest.mark.parametrize('poly', [Chebyshev2D(1, 2), Polynomial2D(2), Legendre2D(1, 2),
                                  Chebyshev1D(5), Legendre1D(5), Polynomial1D(5)])
def test_compound_with_polynomials(poly):
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


def test_indexing_on_instance():
    """Test indexing on compound model instances."""

    m = Gaussian1D(1, 0, 0.1) + Const1D(2)
    assert isinstance(m[0], Gaussian1D)
    assert isinstance(m[1], Const1D)
    # assert isinstance(m['Gaussian1D'], Gaussian1D)
    # assert isinstance(m['Const1D'], Const1D)

    # Test parameter equivalence
    assert m[0].amplitude == 1
    assert m[0].mean == 0
    assert m[0].stddev == 0.1
    assert m[1].amplitude == 2

    # Similar couple of tests, but now where the compound model was created
    # from model instances
    g = Gaussian1D(1, 2, 3, name='g')
    p = Polynomial1D(2, name='p')
    m = g + p
    assert m[0].name == 'g'
    assert m[1].name == 'p'
    assert m['g'].name == 'g'
    assert m['p'].name == 'p'

    # Test negative indexing
    assert isinstance(m[-1], Polynomial1D)
    assert isinstance(m[-2], Gaussian1D)

    with pytest.raises(IndexError):
        m[42]

    with pytest.raises(IndexError):
        m['foobar']

    # Test string slicing
    A = Const1D(1.1, name='A')
    B = Const1D(2.1, name='B')
    C = Const1D(3.1, name='C')
    M = A + B * C
    assert_allclose(M['B':'C'](1), 6.510000000000001)


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
    model = (Gaussian1D(bounds={'stddev': (0, 0.3)}, fixed={'mean': True}) +
             Gaussian1D(fixed={'mean': True}))
    # We have to copy the model before modifying it, otherwise the test fails
    # if it is run twice in a row, because the state of the model instance
    # would be preserved from one run to the next.
    model = deepcopy(model)
    model.map_parameters()
    # Lots of assertions in this test as there are multiple interfaces to
    # parameter constraints

    assert 'stddev_0' in model.bounds
    assert model.bounds['stddev_0'] == (0, 0.3)
    assert model.stddev_0.bounds == (0, 0.3)
    assert 'mean_0' in model.fixed
    assert model.fixed['mean_0'] is True
    assert model.mean_0.fixed is True
    assert 'mean_1' in model.fixed
    assert model.fixed['mean_1'] is True
    assert model.mean_1.fixed is True

    assert model.stddev_0 is model[0].stddev
    # Great, all the constraints were inherited properly
    # Now what about if we update them through the sub-models?
    model.stddev_0.bounds = (0, 0.4)
    assert model[0].stddev.bounds == (0, 0.4)
    assert model[0].bounds['stddev'] == (0, 0.4)

    model.stddev_0.bounds = (0.1, 0.5)
    assert model[0].stddev.bounds == (0.1, 0.5)
    assert model[0].bounds['stddev'] == (0.1, 0.5)

    model[1].mean.fixed = False
    assert model.mean_1.fixed is False
    assert model[1].mean.fixed is False


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
    with pytest.raises(NotImplementedError):
        (shift + model1).inverse

    with pytest.raises(NotImplementedError):
        (model1 & poly).inverse


@pytest.mark.skipif(str("sys.version_info < (2, 7, 3)"))
def test_pickle_compound():
    """
    Regression test for
    https://github.com/astropy/astropy/issues/3867#issuecomment-114547228
    """

    # Test pickling a compound model instance
    g1 = Gaussian1D(1.0, 0.0, 0.1)
    g2 = Gaussian1D([2.0, 3.0], [0.0, 0.0], [0.2, 0.3])
    m = g1 + g2
    m.map_parameters()
    m2 = pickle.loads(pickle.dumps(m))
    assert m.param_names == m2.param_names
    assert m.__class__.__name__ == m2.__class__.__name__
    assert np.all(m.parameters == m2.parameters)
    assert np.all(m(0) == m2(0))


def test_update_parameters():
    offx = Shift(1)
    scl = Scale(2)
    m = offx | scl
    m.map_parameters()
    assert(m(1) == 4)

    offx.offset = 42
    assert(m(1) == 86)

    m.factor_1 = 100
    assert(m(1) == 4300)
    m2 = m | offx
    assert(m2(1) == 4342)


def test_name():
    offx = Shift(1)
    scl = Scale(2)
    m = offx | scl
    scl.name = "scale"
    assert m.submodel_names == ('None_0', 'scale')
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
    with pytest.raises(IndexError):
        g['bozo']
    g1.name = 'bozo'
    assert g['bozo'].mean == 1
    g2.name = 'bozo'
    with pytest.raises(IndexError):
        g['bozo']


@pytest.mark.skipif("not HAS_SCIPY_14")
def test_tabular_in_compound():
    """
    Issue #7411 - evaluate should not change the shape of the output.
    """
    t = Tabular1D(points=([1, 5, 7],), lookup_table=[12, 15, 19],
                  bounds_error=False)
    rot = Rotation2D(2)
    p = Polynomial1D(1)
    x = np.arange(12).reshape((3, 4))
    # Create a compound model which does ot execute Tabular.__call__,
    # but model.evaluate and is followed by a Rotation2D which
    # checks the exact shapes.
    model = p & t | rot
    x1, y1 = model(x, x)
    assert x1.ndim == 2
    assert y1.ndim == 2


def test_bounding_box():
    g = Gaussian2D() + Gaussian2D(2, .5, .1, 2, 3, 0)
    g.bounding_box = ((0, 1), (0, .5))
    y, x = np.mgrid[0:10, 0:10]
    y = y / 3.
    x = x / 3.
    val = g(x, y, with_bounding_box=True)
    compare = np.array([
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
         np.nan, np.nan, np.nan, np.nan, np.nan]])
    mask = ~np.isnan(val)
    assert_allclose(val[mask], compare[mask])
    val2 = g(x+2, y+2, with_bounding_box=True)
    assert(np.isnan(val2).sum() == 100)
    val3 = g(.1, .1, with_bounding_box=True)


@pytest.mark.skipif("not HAS_SCIPY_14")
def test_bounding_box_with_units():
    points = np.arange(5) * u.pix
    lt = np.arange(5) * u.AA
    t = Tabular1D(points, lt)

    assert(t(1 * u.pix, with_bounding_box=True) == 1. * u.AA)

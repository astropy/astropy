# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import inspect
from copy import deepcopy

import numpy as np

from numpy.testing.utils import (assert_allclose, assert_array_equal,
                                 assert_almost_equal)

from ...extern.six.moves import cPickle as pickle
from ...tests.helper import pytest

from ..core import Model, ModelDefinitionError, set_compound_model
from ..parameters import Parameter
from ..models import (Const1D, Shift, Scale, Rotation2D, Gaussian1D,
                      Gaussian2D, Polynomial1D, Polynomial2D,
                      Chebyshev2D, Legendre2D, Chebyshev1D, Legendre1D,
                      AffineTransformation2D, Identity, Mapping)
from ..altcompound import _AltCompoundModel


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, [5.0, 5.0]),
                          (lambda x, y: x - y, [-1.0, -1.0]),
                          (lambda x, y: x * y, [6.0, 6.0]),
                          (lambda x, y: x / y, [2.0 / 3.0, 2.0 / 3.0]),
                          (lambda x, y: x ** y, [8.0, 8.0])])
def test_model_set(expr, result):
    set_compound_model('lite')
    s = expr(Const1D((2, 2), n_models=2), Const1D((3, 3), n_models=2))
    out = s(0, model_set_axis=False)
    set_compound_model('regular')
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
    #set_compound_model('lite')
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

    set_compound_model('lite')
    s = expr(Const1D(2), Const1D(3))

    assert isinstance(s, _AltCompoundModel)
    assert s.n_inputs == 1
    assert s.n_outputs == 1

    out = s(0)
    assert out == result
    assert isinstance(out, float)
    set_compound_model('regular')


def test_simple_two_model_compose_1d():
    """
    Shift and Scale are two of the simplest models to test model composition
    with.
    """

    set_compound_model('lite')
    S1 = Shift(2) | Scale(3)  # First shift then scale
    assert isinstance(S1, _AltCompoundModel)
    assert S1.n_inputs == 1
    assert S1.n_outputs == 1
    assert S1(1) == 9.0

    S2 = Scale(2) | Shift(3)  # First scale then shift
    assert isinstance(S2, _AltCompoundModel)
    assert S2.n_inputs == 1
    assert S2.n_outputs == 1
    assert S2(1) == 5.0

    set_compound_model('regular')
    # Test with array inputs
    assert_array_equal(S2([1, 2, 3]), [5.0, 7.0, 9.0])


def test_simple_two_model_compose_2d():
    """
    A simple example consisting of two rotations.
    """

    set_compound_model('lite')
    r1 = Rotation2D(45) | Rotation2D(45)
    assert isinstance(r1, _AltCompoundModel)
    assert r1.n_inputs == 2
    assert r1.n_outputs == 2
    assert_allclose(r1(0, 1), (-1, 0), atol=1e-10)

    r2 = Rotation2D(90) | Rotation2D(90)  # Rotate twice by 90 degrees
    assert_allclose(r2(0, 1), (0, -1), atol=1e-10)

    # Compose R with itself to produce 4 rotations
    r3 = r1 | r1

    set_compound_model('regular')
    assert_allclose(r3(0, 1), (0, -1), atol=1e-10)

def xtest_n_submodels():
    """
    Test that CompoundModel.n_submodels properly returns the number
    of components.
    """
    set_compound_model('lite')
    g2 = Gaussian1D() + Gaussian1D()
    assert g2.n_submodels() == 2

    g3 = g2 + Gaussian1D()
    assert g3.n_submodels() == 3

    g5 = g3 | g2
    assert g5.n_submodels() == 5

    g7 = g5 / g2
    assert g7.n_submodels() == 7

    # make sure it works as class method
    p = Polynomial1D + Polynomial1D

    set_compound_model('regular')
    assert p.n_submodels() == 2

def xtest_expression_formatting():
    """
    Test that the expression strings from compound models are formatted
    correctly.
    """

    set_compound_model('lite')
    # For the purposes of this test it doesn't matter a great deal what
    # model(s) are used in the expression, I don't think
    G = Gaussian1D
    G2 = Gaussian2D

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
    set_compound_model('regular')



def test_basic_compound_inverse():
    """
    Test basic inversion of compound models in the limited sense supported for
    models made from compositions and joins only.
    """

    set_compound_model('lite')
    t = (Shift(2) & Shift(3)) | (Scale(2) & Scale(3)) | Rotation2D(90)
    set_compound_model('regular')
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

    #set_compound_model('lite')
    with pytest.raises(NotImplementedError):
        model.inverse


def test_mapping_basic_permutations():
    """
    Tests a couple basic examples of the Mapping model--specifically examples
    that merely permute the outputs.
    """

    set_compound_model('lite')
    x, y = Rotation2D(90)(1, 2)

    rs = Rotation2D(90) | Mapping((1, 0))
    x_prime, y_prime = rs(1, 2)
    assert_allclose((x, y), (y_prime, x_prime))

    # A more complicated permutation
    m = Rotation2D(90) & Scale(2)
    x, y, z = m(1, 2, 3)

    ms = m | Mapping((2, 0, 1))
    x_prime, y_prime, z_prime = ms(1, 2, 3)
    set_compound_model('regular')
    assert_allclose((x, y, z), (y_prime, z_prime, x_prime))


def test_mapping_inverse():
    """Tests inverting a compound model that includes a `Mapping`."""

    set_compound_model('lite')
    rs1 = Rotation2D(12.1) & Scale(13.2)
    rs2 = Rotation2D(14.3) & Scale(15.4)

    # Rotates 2 of the coordinates and scales the third--then rotates on a
    # different axis and scales on the axis of rotation.  No physical meaning
    # here just a simple test
    m = rs1 | Mapping([2, 0, 1]) | rs2

    set_compound_model('regular')
    assert_allclose((0, 1, 2), m.inverse(*m(0, 1, 2)), atol=1e-08)


def test_identity_input():
    """
    Test a case where an Identity (or Mapping) model is the first in a chain
    of composite models and thus is responsible for handling input broadcasting
    properly.

    Regression test for https://github.com/astropy/astropy/pull/3362
    """

    set_compound_model('lite')
    ident1 = Identity(1)
    shift = Shift(1)
    rotation = Rotation2D(angle=90)
    model = ident1 & shift | rotation
    assert_allclose(model(1, 2), [-3.0, 1.0])

    set_compound_model('regular')


def test_invalid_operands():
    """
    Test that certain operators do not work with models whose inputs/outputs do
    not match up correctly.
    """

    set_compound_model('lite')
    with pytest.raises(ModelDefinitionError):
        Rotation2D(90) | Gaussian1D(1, 0, 0.1)

    with pytest.raises(ModelDefinitionError):
        Rotation2D(90) + Gaussian1D(1, 0, 0.1)
    set_compound_model('regular')

@pytest.mark.parametrize('poly', [Chebyshev2D(1, 2), Polynomial2D(2), Legendre2D(1, 2),
                                  Chebyshev1D(5), Legendre1D(5), Polynomial1D(5)])
def test_compound_with_polynomials(poly):
    """
    Tests that polynomials are scaled when used in compound models.
    Issue #3699
    """
    set_compound_model('lite')
    poly.parameters = [1, 2, 3, 4, 1, 2]
    shift = Shift(3)
    model = poly | shift
    x, y = np.mgrid[:20, :37]
    result_compound = model(x, y)
    result = shift(poly(x, y))
    set_compound_model('regular')
    assert_allclose(result, result_compound)



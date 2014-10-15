# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from numpy.testing.utils import (assert_allclose, assert_array_equal,
                                 assert_almost_equal)

from ...tests.helper import pytest

from ..core import Model, Identity, Mapping
from ..models import (Const1D, Shift, Scale, Rotation2D, Gaussian1D,
                      Polynomial1D, Polynomial2D, AffineTransformation2D)


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, 5.0),
                          (lambda x, y: x - y, -1.0),
                          (lambda x, y: x * y, 6.0),
                          (lambda x, y: x / y, 2.0 / 3.0),
                          (lambda x, y: x ** y, 8.0)])
def test_two_model_class_arithmetic_1d(expr, result):
    # Const1D is perhaps the simplest model to test basic arithmetic with.
    # TODO: Should define more tests later on for more complicated
    # combinations of models

    S = expr(Const1D, Const1D)

    assert issubclass(S, Model)
    assert S.n_inputs == 1
    assert S.n_outputs == 1

    # Initialize an instance of the model, providing values for the two
    # "amplitude" parameters
    s = S(2, 3)

    # It shouldn't matter what input we evaluate on since this is a constant
    # function
    out = s(0)
    assert out == result
    assert isinstance(out, float)


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

    assert isinstance(s, Model)
    assert s.n_inputs == 1
    assert s.n_outputs == 1

    out = s(0)
    assert out == result
    assert isinstance(out, float)


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, 5.0),
                          (lambda x, y: x - y, -1.0),
                          (lambda x, y: x * y, 6.0),
                          (lambda x, y: x / y, 2.0 / 3.0),
                          (lambda x, y: x ** y, 8.0)])
def test_two_model_mixed_arithmetic_1d(expr, result):
    """
    Like test_two_model_class_arithmetic_1d, but creates a new model from an
    expression of one model class with one model instance (and vice-versa).
    """

    S1 = expr(Const1D, Const1D(3))
    S2 = expr(Const1D(2), Const1D)

    for cls in (S1, S2):
        assert issubclass(cls, Model)
        assert cls.n_inputs == 1
        assert cls.n_outputs == 1

    # Takes only one argument--the one unspecified amplitude
    s1 = S1(2)
    s2 = S2(3)

    for out in (s1(0), s2(0)):
        assert out == result
        assert isinstance(out, float)


def test_simple_two_model_class_compose_1d():
    """
    Shift and Scale are two of the simplest models to test model composition
    with.
    """

    S1 = Shift | Scale  # First shift then scale
    assert issubclass(S1, Model)
    assert S1.n_inputs == 1
    assert S1.n_outputs == 1

    s1 = S1(2, 3)  # Shift by 2 and scale by 3
    assert s1(1) == 9.0

    S2 = Scale | Shift  # First scale then shift
    assert issubclass(S2, Model)
    assert S2.n_inputs == 1
    assert S2.n_outputs == 1

    s2 = S2(2, 3)  # Scale by 2 then shift by 3
    assert s2(1) == 5.0

    # Test with array inputs
    assert_array_equal(s2([1, 2, 3]), [5.0, 7.0, 9.0])


def test_simple_two_model_class_compose_2d():
    """
    A simple example consisting of two rotations.
    """

    R = Rotation2D | Rotation2D
    assert issubclass(R, Model)
    assert R.n_inputs == 2
    assert R.n_outputs == 2

    r1 = R(45, 45)  # Rotate twice by 45 degrees
    assert_allclose(r1(0, 1), (-1, 0), atol=1e-10)

    r2 = R(90, 90)  # Rotate twice by 90 degrees
    assert_allclose(r2(0, 1), (0, -1), atol=1e-10)

    # Compose R with itself to produce 4 rotations
    R2 = R | R

    r3 = R2(45, 45, 45, 45)
    assert_allclose(r3(0, 1), (0, -1), atol=1e-10)


class TestCompositeLegacy(object):
    """
    Tests inspired by the original _CompositeModel tests in test_core.py,
    this implements the equivalent tests implemented in the new framework.

    Note: These aren't *exactly* the same as the original tests, as they used
    overly trivial models (polynomials with all coeffs 0).
    """

    def setup_class(self):
        self.y, self.x = np.mgrid[:5, :5]

    def test_single_array_input(self):
        p1 = Polynomial1D(3, c0=1, c1=2, c2=3, c3=4)
        p2 = Polynomial1D(3, c0=2, c1=3, c2=4, c3=5)
        m = p1 | p2
        assert_almost_equal(p2(p1(self.x)), m(self.x))

    def test_labeledinput_1(self):
        # Note: No actual use of LabeledInput in this test; this just uses the
        # same name for symmetry with the old tests
        p1 = Polynomial1D(3, c0=1, c1=2, c2=3, c3=4)
        p2 = Polynomial2D(3, c0_0=1, c2_0=2, c0_1=3, c2_1=4)
        m = p2 | p1
        assert_almost_equal(p1(p2(self.x, self.y)), m(self.x, self.y))

    def test_labledinput_2(self):
        rot = Rotation2D(angle=23.4)
        offx = Shift(-2)
        offy = Shift(1.2)
        m = rot | (offx & Identity(1)) | (Identity(1) & offy)

        x, y = rot(self.x, self.y)
        x = offx(x)
        y = offy(y)

        assert_almost_equal(x, m(self.x, self.y)[0])
        assert_almost_equal(y, m(self.x, self.y)[1])

        a = np.deg2rad(23.4)
        # For kicks
        matrix = [[np.cos(a), -np.sin(a)],
                  [np.sin(a), np.cos(a)]]
        x, y = AffineTransformation2D(matrix, [-2, 1.2])(self.x, self.y)
        assert_almost_equal(x, m(self.x, self.y)[0])
        assert_almost_equal(y, m(self.x, self.y)[1])

    def test_multiple_input(self):
        """
        Despite the name, this actually tests inverting composite models,
        which is not yet supported in the new framework (but should be).
        """

        rot = Rotation2D(-60)
        m = rot | rot
        xx, yy = m(self.x, self.y)
        x0, y0 = m.inverse(xx, yy)
        assert_almost_equal(x0, self.x)
        assert_almost_equal(y0, self.y)


def test_expression_formatting():
    """
    Test that the expression strings from compound models are formatted
    correctly.
    """

    # For the purposes of this test it doesn't matter a great deal what
    # model(s) are used in the expression, I don't think
    G = Gaussian1D

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

    M = G + (G | G )
    assert M._format_expression() == '[0] + ([1] | [2])'

    M = G & G | G
    assert M._format_expression() == '[0] & [1] | [2]'

    M = G & (G | G)
    assert M._format_expression() == '[0] & ([1] | [2])'


def test_indexing_on_class():
    """
    Test indexing on compound model class objects, including cases where the
    submodels are classes, as well as instances, or both.
    """

    g = Gaussian1D(1, 2, 3, name='g')
    p = Polynomial1D(2, name='p')

    M = Gaussian1D + Const1D
    assert M[0] is Gaussian1D
    assert M[1] is Const1D
    assert M['Gaussian1D'] is M[0]
    assert M['Const1D'] is M[1]

    M = Gaussian1D + p
    assert M[0] is Gaussian1D
    assert M[1] is p
    assert M['Gaussian1D'] is M[0]
    assert M['p'] is M[1]

    M = g + p
    assert M[0] is g
    assert M[1] is p
    assert M['g'] is M[0]
    assert M['p'] is M[1]

    # Test negative indexing
    assert M[-1] is p
    assert M[-2] is g

    with pytest.raises(IndexError):
        M[42]

    with pytest.raises(IndexError):
        M['foobar']


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

    RS = Rotation2D | Mapping((1, 0))
    x_prime, y_prime = RS(90)(1, 2)
    assert_allclose((x, y), (y_prime, x_prime))

    # A more complicated permutation
    M = Rotation2D & Scale
    m = M(90, 2)
    x, y, z = m(1, 2, 3)

    MS = M | Mapping((2, 0, 1))
    ms = MS(90, 2)
    x_prime, y_prime, z_prime = ms(1, 2, 3)
    assert_allclose((x, y, z), (y_prime, z_prime, x_prime))


def test_mapping_inverse():
    """Tests inverting a compound model that includes a `Mapping`."""

    RS = Rotation2D & Scale

    # Rotates 2 of the coordinates and scales the third--then rotates on a
    # different axis and scales on the axis of rotation.  No physical meaning
    # here just a simple test
    M = RS | Mapping([2, 0, 1]) | RS

    m = M(12.1, 13.2, 14.3, 15.4)

    assert_allclose((0, 1, 2), m.inverse(*m(0, 1, 2)), atol=1e-08)

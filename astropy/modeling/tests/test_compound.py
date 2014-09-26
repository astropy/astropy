# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from numpy.testing.utils import assert_allclose, assert_array_equal

from ...tests.helper import pytest

from ..core import Model
from ..models import Const1D, Shift, Scale, Rotation2D, Gaussian1D


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

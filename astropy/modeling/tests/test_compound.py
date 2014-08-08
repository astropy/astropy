# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ...tests.helper import pytest

from ..core import Model
from ..models import Const1D, Shift, Scale


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
    assert s(0) == result


def test_two_model_class_compose_1d():
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

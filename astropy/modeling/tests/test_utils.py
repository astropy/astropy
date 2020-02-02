# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name
import operator

import numpy as np

from astropy.tests.helper import catch_warnings
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.modeling.utils import ExpressionTree as ET, ellipse_extent
from astropy.modeling.models import Ellipse2D


def test_traverse_postorder_duplicate_subtrees():
    """
    Regression test for a bug in `ExpressionTree.traverse_postorder`
    where given an expression like ``(1 + 2) + (1 + 2)`` where the two proper
    subtrees are actually the same object.
    """
    with catch_warnings(AstropyDeprecationWarning):
        subtree = ET('+', ET(1), ET(2))
        tree = ET('+', subtree, subtree)
    traversal = [n.value for n in tree.traverse_postorder()]
    assert traversal == [1, 2, '+', 1, 2, '+', '+']


def test_tree_evaluate_subexpression():
    """Test evaluating a subexpression from an expression tree."""

    operators = {'+': operator.add, '-': operator.sub, '*': operator.mul,
                 '/': operator.truediv, '**': operator.pow}
    # The full expression represented by this tree is:
    # 1.0 + 2 - 3 * 4 / 5 ** 6 (= 2.999232 if you must know)
    with catch_warnings(AstropyDeprecationWarning):
        tree = ET('+', ET(1.0), ET('-', ET(2.0),
                                   ET('*', ET(3.0), ET('/', ET(4.0),
                                                       ET('**', ET(5.0), ET(6.0))))))

    def test_slice(start, stop, expected):
        assert np.allclose(tree.evaluate(operators, start=start, stop=stop),
                           expected)

    assert tree.evaluate(operators) == (1.0 + 2.0 - 3.0 * 4.0 / 5.0 ** 6.0)
    test_slice(0, 5, (1.0 + 2.0 - 3.0 * 4.0 / 5.0))
    test_slice(0, 4, (1.0 + 2.0 - 3.0 * 4.0))
    test_slice(0, 3, (1.0 + 2.0 - 3.0))
    test_slice(0, 2, (1.0 + 2.0))
    test_slice(0, 1, 1.0)

    test_slice(1, 6, (2.0 - 3.0 * 4.0 / 5.0 ** 6.0))
    test_slice(1, 5, (2.0 - 3.0 * 4.0 / 5.0))
    test_slice(1, 4, (2.0 - 3.0 * 4.0))
    test_slice(1, 3, (2.0 - 3.0))
    test_slice(1, 2, 2.0)

    test_slice(2, 6, (3.0 * 4.0 / 5.0 ** 6.0))
    test_slice(2, 5, (3.0 * 4.0 / 5.0))
    test_slice(2, 4, (3.0 * 4.0))
    test_slice(2, 3, 3.0)

    test_slice(3, 6, (4.0 / 5.0 ** 6.0))
    test_slice(3, 5, (4.0 / 5.0))
    test_slice(3, 4, 4.0)

    test_slice(4, 6, (5.0 ** 6.0))
    test_slice(4, 5, 5.0)

    test_slice(5, 6, 6.0)


def test_ellipse_extent():
    # Test this properly bounds the ellipse

    imshape = (100, 100)
    coords = y, x = np.indices(imshape)

    amplitude = 1
    x0 = 50
    y0 = 50
    a = 30
    b = 10
    theta = np.pi / 4

    model = Ellipse2D(amplitude, x0, y0, a, b, theta)

    dx, dy = ellipse_extent(a, b, theta)

    limits = ((y0 - dy, y0 + dy), (x0 - dx, x0 + dx))

    model.bounding_box = limits

    actual = model.render(coords=coords)

    expected = model(x, y)

    # Check that the full ellipse is captured
    np.testing.assert_allclose(expected, actual, atol=0, rtol=1)

    # Check the bounding_box isn't too large
    limits = np.array(limits).flatten()
    for i in [0, 1]:
        s = actual.sum(axis=i)
        diff = np.abs(limits[2 * i] - np.where(s > 0)[0][0])
        assert diff < 1

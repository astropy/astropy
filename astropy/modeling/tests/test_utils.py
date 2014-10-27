# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import operator

import numpy as np

from ..utils import ExpressionTree as ET


def test_traverse_postorder_duplicate_subtrees():
    """
    Regression test for a bug in `ExpressionTree.traverse_postorder`
    where given an expression like ``(1 + 2) + (1 + 2)`` where the two proper
    subtrees are actually the same object.
    """

    subtree = ET('+', ET(1), ET(2))
    tree = ET('+', subtree, subtree)
    traversal = [n.value for n in tree.traverse_postorder()]
    assert traversal == [1, 2, '+', 1, 2, '+', '+']


# TODO: It might prove useful to implement a simple expression parser to build
# trees; this would be easy and might find use elsewhere
def test_tree_evaluate_subexpression():
    """Test evaluating a subexpression from an expression tree."""

    operators = {'+': operator.add, '-': operator.sub, '*': operator.mul,
                 '/': operator.truediv, '**': operator.pow}
    # The full expression represented by this tree is:
    # 1.0 + 2 - 3 * 4 / 5 ** 6 (= 2.999232 if you must know)
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

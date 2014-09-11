# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

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

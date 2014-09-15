# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides utility functions for the models package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import textwrap

from collections import deque

import numpy as np

from ..extern.six.moves import xrange, zip_longest


__all__ = ['ExpressionTree', 'check_broadcast', 'poly_map_domain', 'comb']


class ExpressionTree(object):
    __slots__ = ['left', 'right', 'value']

    def __init__(self, value, left=None, right=None):
        self.value = value
        self.left = left

        # Two subtrees can't be the same *object* or else traverse_postorder
        # breaks, so we just always copy the right subtree to subvert that.
        if right is not None and left is right:
            right = right.copy()

        self.right = right

    @property
    def isleaf(self):
        return self.left is None and self.right is None

    def traverse_preorder(self):
        stack = deque([self])
        while stack:
            node = stack.pop()
            yield node

            if node.right is not None:
                stack.append(node.right)
            if node.left is not None:
                stack.append(node.left)

    def traverse_inorder(self):
        stack = deque()
        node = self
        while stack or node is not None:
            if node is not None:
                stack.append(node)
                node = node.left
            else:
                node = stack.pop()
                yield node
                node = node.right

    def traverse_postorder(self):
        stack = deque([self])
        last = None
        while stack:
            node = stack[-1]
            if last is None or node is last.left or node is last.right:
                if node.left is not None:
                    stack.append(node.left)
                elif node.right is not None:
                    stack.append(node.right)
            elif node.left is last and node.right is not None:
                stack.append(node.right)
            else:
                yield stack.pop()
            last = node

    def evaluate(self, operators, getter=None):
        """Evaluate the expression represented by this tree.

        ``Operators`` should be a dictionary mapping operator names ('tensor',
        'product', etc.) to a function that implements that operator for the
        correct number of operands.

        If given, ``getter`` is a function evaluated on each *leaf* node's
        value before applying the operator between them.  This could be used,
        for example, to operate on an attribute of the node values rather than
        directly on the node values.
        """

        operands = deque()

        for node in self.traverse_postorder():
            if node.isleaf:
                # For a "tree" containing just a single operator at the root
                operands.append(getter(node.value))
            else:
                operator = operators[node.value]
                right = operands.pop()
                left = operands.pop()
                operands.append(operator(left, right))

        return operands.pop()

    def copy(self):
        # Hopefully this won't blow the stack for any practical case; if such a
        # case arises that this won't work then I suppose we can find an
        # iterative approach.

        children = []
        for child in (self.left, self.right):
            if isinstance(child, ExpressionTree):
                children.append(child.copy())
            else:
                children.append(child)

        return self.__class__(self.value, left=children[0], right=children[1])

    def format_expression(self, operator_precedence, format_leaf=None):
        leaf_idx = 0
        operands = deque()

        if format_leaf is None:
            format_leaf = lambda i, l: '[{0}]'.format(i)

        for node in self.traverse_postorder():
            if node.isleaf:
                operands.append(format_leaf(leaf_idx, node))
                leaf_idx += 1
                continue

            oper_order = operator_precedence[node.value]
            right = operands.pop()
            left = operands.pop()

            if (node.left is not None and not node.left.isleaf and
                    operator_precedence[node.left.value] < oper_order):
                left = '({0})'.format(left)
            if (node.right is not None and not node.right.isleaf and
                    operator_precedence[node.right.value] < oper_order):
                right = '({0})'.format(right)

            operands.append(' '.join((left, node.value, right)))

        return ''.join(operands)

    # TODO: This could still use a lot of improvement; in particular the trees
    # it outputs are often too wide, and could be made more compactly.  More
    # formatting control would be useful too.
    def format_tree_ascii(self, format_leaf=lambda i, l: str(l)):
        """
        Format the tree using an ASCII character representation.

        Parameters
        ----------
        format_leaf : callable
            A function of a single argument which, given a node value,
            returns a string representing that node in the tree display.
        """

        stack = deque()
        leaf_idx = 0

        for node in self.traverse_postorder():
            if node.isleaf:
                text = format_leaf(leaf_idx, node.value)
                stack.append((len(text), [text]))
                leaf_idx += 1
                continue

            right_width, right = stack.pop()
            left_width, left = stack.pop()

            if left_width > right_width:
                right = [r.center(left_width) for r in right]
                child_width = left_width
            elif right_width > left_width:
                left = [l.center(right_width) for l in left]
                child_width = right_width
            else:
                child_width = left_width  # without loss of generality

            root = '[{0}]'.format(node.value)
            spine = '/   \\'
            fill = ' ' * len(spine)
            width = 2 * child_width + len(fill)

            offset = child_width - (child_width % 2)
            spine = spine.center(width - 2 - offset, '_').center(width)
            root = root.center(width)

            lines = [root, spine]
            for l, r in zip_longest(left, right,
                                    fillvalue=' ' * child_width):
                lines.append(l + fill + r)

            stack.append((width, lines))

        return textwrap.dedent('\n'.join(stack[0][1]))


def make_binary_operator_eval(oper, f, g):
    """
    Given a binary operator (as a callable of two arguments) ``oper`` and
    two callables ``f`` and ``g`` which accept the same arguments,
    returns a *new* function that takes the same arguments as ``f`` and ``g``,
    but passes the outputs of ``f`` and ``g`` in the given ``oper``.

    ``f`` and ``g`` are assumed to return tuples (which may be 1-tuples).  The
    given operator is applied element-wise to tuple outputs).

    Example
    -------

    >>> from operator import add
    >>> def prod(x, y):
    ...     return (x * y,)
    ...
    >>> sum_of_prod = make_binary_operator_eval(add, prod, prod)
    >>> sum_of_prod(3, 5)
    (30,)
    """

    return lambda *args: tuple(oper(x, y)
                               for x, y in zip(f(*args), g(*args)))


class IncompatibleShapeError(ValueError):
    def __init__(self, shape_a, shape_a_idx, shape_b, shape_b_idx):
        super(IncompatibleShapeError, self).__init__(
                shape_a, shape_a_idx, shape_b, shape_b_idx)


def check_broadcast(*shapes):
    """
    Determines whether two or more Numpy arrays can be broadcast with each
    other based on their shape tuple alone.

    Parameters
    ----------
    *shapes : tuple
        All shapes to include in the comparison.  If only one shape is given it
        is passed through unmodified.  If no shapes are given returns an empty
        `tuple`.

    Returns
    -------
    broadcast : `tuple`
        If all shapes are mutually broadcastable, returns a tuple of the full
        broadcast shape.
    """

    if len(shapes) == 0:
        return ()
    elif len(shapes) == 1:
        return shapes[0]

    reversed_shapes = (reversed(shape) for shape in shapes)

    full_shape = []

    for dims in zip_longest(*reversed_shapes, fillvalue=1):
        max_dim = 1
        max_dim_idx = None
        for idx, dim in enumerate(dims):
            if dim == 1:
                continue

            if max_dim == 1:
                # The first dimension of size greater than 1
                max_dim = dim
                max_dim_idx = idx
            elif dim != max_dim:
                raise IncompatibleShapeError(
                    shapes[max_dim_idx], max_dim_idx, shapes[idx], idx)

        full_shape.append(max_dim)

    return tuple(full_shape[::-1])


def poly_map_domain(oldx, domain, window):
    """
    Map domain into window by shifting and scaling.

    Parameters
    ----------
    oldx : array
          original coordinates
    domain : list or tuple of length 2
          function domain
    window : list or tuple of length 2
          range into which to map the domain
    """
    domain = np.array(domain, dtype=np.float64)
    window = np.array(window, dtype=np.float64)
    scl = (window[1] - window[0]) / (domain[1] - domain[0])
    off = (window[0] * domain[1] - window[1] * domain[0]) / (domain[1] - domain[0])
    return off + scl * oldx


def comb(N, k):
    """
    The number of combinations of N things taken k at a time.

    Parameters
    ----------
    N : int, array
        Number of things.
    k : int, array
        Number of elements taken.

    """
    if (k > N) or (N < 0) or (k < 0):
        return 0
    val = 1
    for j in xrange(min(k, N - k)):
        val = (val * (N - j)) / (j + 1)
    return val


def array_repr_oneline(array):
    """
    Represents a multi-dimensional Numpy array flattened onto a single line.
    """

    r = np.array2string(array, separator=',', suppress_small=True)
    return ' '.join(l.strip() for l in r.splitlines())


def combine_labels(left, right):
    """
    For use with the join operator &: Combine left input/output labels with
    right input/output labels.

    If none of the labels conflict then this just returns a sum of tuples.
    However if *any* of the labels conflict, this appends '0' to the left-hand
    labels and '1' to the right-hand labels so there is no ambiguity).
    """

    if set(left).intersection(right):
        left = tuple(l + '0' for l in left)
        right = tuple(r + '1' for r in right)

    return left + right


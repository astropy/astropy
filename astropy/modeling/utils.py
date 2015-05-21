# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides utility functions for the models package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from collections import deque, MutableMapping

import numpy as np

from ..extern import six
from ..extern.six.moves import xrange, zip_longest


__all__ = ['ExpressionTree', 'AliasDict', 'check_broadcast',
           'poly_map_domain', 'comb']


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

    def evaluate(self, operators, getter=None, start=0, stop=None):
        """Evaluate the expression represented by this tree.

        ``Operators`` should be a dictionary mapping operator names ('tensor',
        'product', etc.) to a function that implements that operator for the
        correct number of operands.

        If given, ``getter`` is a function evaluated on each *leaf* node's
        value before applying the operator between them.  This could be used,
        for example, to operate on an attribute of the node values rather than
        directly on the node values.  The ``getter`` is passed both the index
        of the leaf (a count starting at 0 that is incremented after each leaf
        is found) and the leaf node itself.

        The ``start`` and ``stop`` arguments allow evaluating a sub-expression
        within the expression tree.

        TODO: Document this better.
        """

        stack = deque()

        if getter is None:
            getter = lambda idx, value: value

        if start is None:
            start = 0

        leaf_idx = 0
        for node in self.traverse_postorder():
            if node.isleaf:
                # For a "tree" containing just a single operator at the root
                # Also push the index of this leaf onto the stack, which will
                # prove useful for evaluating subexpressions
                stack.append((getter(leaf_idx, node.value), leaf_idx))
                leaf_idx += 1
            else:
                operator = operators[node.value]

                if len(stack) < 2:
                    # Skip this operator if there are not enough operands on
                    # the stack; this can happen if some operands were skipped
                    # when evaluating a sub-expression
                    continue

                right = stack.pop()
                left = stack.pop()
                operands = []

                for operand in (left, right):
                    # idx is the leaf index; -1 if not a leaf node
                    if operand[-1] == -1:
                        operands.append(operand)
                    else:
                        operand, idx = operand
                        if start <= idx and (stop is None or idx < stop):
                            operands.append((operand, idx))

                if len(operands) == 2:
                    # evaluate the operator with the given operands and place
                    # the result on the stack (with -1 for the "leaf index"
                    # since this result is not a leaf node
                    left, right = operands
                    stack.append((operator(left[0], right[0]), -1))
                elif len(operands) == 0:
                    # Just push the left one back on the stack
                    # TODO: Explain and/or refactor this better
                    # This is here because even if both operands were "skipped"
                    # due to being outside the (start, stop) range, we've only
                    # skipped one operator.  But there should be at least 2
                    # operators involving these operands, so we push the one
                    # from the left back onto the stack so that the next
                    # operator will be skipped as well.  Should probably come
                    # up with an easier to follow way to write this algorithm
                    stack.append(left)
                else:
                    # one or more of the operands was not included in the
                    # sub-expression slice, so don't evaluate the operator;
                    # instead place left over operands (if any) back on the
                    # stack for later use
                    stack.extend(operands)

        return stack.pop()[0]

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


class AliasDict(MutableMapping):
    """
    Creates a `dict` like object that wraps an existing `dict` or other
    `MutableMapping`, along with a `dict` of *key aliases* that translate
    between specific keys in this dict to different keys in the underlying
    dict.

    In other words, keys that do not have an associated alias are accessed and
    stored like a normal `dict`.  However, a key that has an alias is accessed
    and stored to the "parent" dict via the alias.

    Parameters
    ----------
    parent : dict-like
        The parent `dict` that aliased keys and accessed from and stored to.

    aliases : dict-like
        Maps keys in this dict to their associated keys in the parent dict.

    Examples
    --------

    >>> parent = {'a': 1, 'b': 2, 'c': 3}
    >>> aliases = {'foo': 'a', 'bar': 'c'}
    >>> alias_dict = AliasDict(parent, aliases)
    >>> alias_dict['foo']
    1
    >>> alias_dict['bar']
    3

    Keys in the original parent dict are not visible if they were not
    aliased::

    >>> alias_dict['b']
    Traceback (most recent call last):
    ...
    KeyError: 'b'

    Likewise, updates to aliased keys are reflected back in the parent dict::

    >>> alias_dict['foo'] = 42
    >>> alias_dict['foo']
    42
    >>> parent['a']
    42

    However, updates/insertions to keys that are *not* aliased are not
    reflected in the parent dict::

    >>> alias_dict['qux'] = 99
    >>> alias_dict['qux']
    99
    >>> 'qux' in parent
    False

    In particular, updates on the `AliasDict` to a key that is equal to
    one of the aliased keys in the parent dict does *not* update the parent
    dict.  For example, ``alias_dict`` aliases ``'foo'`` to ``'a'``.  But
    assigning to a key ``'a'`` on the `AliasDict` does not impact the
    parent::

    >>> alias_dict['a'] = 'nope'
    >>> alias_dict['a']
    'nope'
    >>> parent['a']
    42
    """

    _store_type = dict
    """
    Subclasses may override this to use other mapping types as the underlying
    storage, for example an `OrderedDict`.  However, even in this case
    additional work may be needed to get things like the ordering right.
    """

    def __init__(self, parent, aliases):
        self._parent = parent
        self._store = self._store_type()
        self._aliases = dict(aliases)

    def __getitem__(self, key):
        if key in self._aliases:
            try:
                return self._parent[self._aliases[key]]
            except KeyError:
                raise KeyError(key)

        return self._store[key]

    def __setitem__(self, key, value):
        if key in self._aliases:
            self._parent[self._aliases[key]] = value
        else:
            self._store[key] = value

    def __delitem__(self, key):
        if key in self._aliases:
            try:
                del self._parent[self._aliases[key]]
            except KeyError:
                raise KeyError(key)
        else:
            del self._store[key]

    def __iter__(self):
        """
        First iterates over keys from the parent dict (if the aliased keys are
        present in the parent), followed by any keys in the local store.
        """

        for key, alias in six.iteritems(self._aliases):
            if alias in self._parent:
                yield key

        for key in self._store:
            yield key

    def __len__(self):
        # TODO:
        # This could be done more efficiently, but at present the use case for
        # it is narrow if non-existent.
        return len(list(iter(self)))

    def __repr__(self):
        # repr() just like any other dict--this should look transparent
        store_copy = self._store_type()
        for key, alias in six.iteritems(self._aliases):
            if alias in self._parent:
                store_copy[key] = self._parent[alias]

        store_copy.update(self._store)

        return repr(store_copy)


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

    return lambda inputs, params: \
            tuple(oper(x, y) for x, y in zip(f(inputs, params),
                                             g(inputs, params)))


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

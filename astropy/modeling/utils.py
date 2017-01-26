# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides utility functions for the models package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from collections import deque, MutableMapping

import numpy as np

from ..extern import six
from ..extern.six.moves import range, zip_longest, zip

from ..utils import isiterable, check_broadcast
from ..utils.compat.funcsigs import signature


__all__ = ['ExpressionTree', 'AliasDict', 'check_broadcast',
           'poly_map_domain', 'comb', 'ellipse_extent']


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

    def __getstate__(self):
        # For some reason the default pickle protocol on Python 2 does not just
        # do this.  On Python 3 it's not a problem.
        return dict((slot, getattr(self, slot)) for slot in self.__slots__)

    def __setstate__(self, state):
        for slot, value in state.items():
            setattr(self, slot, value)

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


class _BoundingBox(tuple):
    """
    Base class for models with custom bounding box templates (methods that
    return an actual bounding box tuple given some adjustable parameters--see
    for example `~astropy.modeling.models.Gaussian1D.bounding_box`).

    On these classes the ``bounding_box`` property still returns a `tuple`
    giving the default bounding box for that instance of the model.  But that
    tuple may also be a subclass of this class that is callable, and allows
    a new tuple to be returned using a user-supplied value for any adjustable
    parameters to the bounding box.
    """

    _model = None

    def __new__(cls, input_, _model=None):
        self = super(_BoundingBox, cls).__new__(cls, input_)
        if _model is not None:
            # Bind this _BoundingBox (most likely a subclass) to a Model
            # instance so that its __call__ can access the model
            self._model = _model

        return self

    def __call__(self, *args, **kwargs):
        raise NotImplementedError(
            "This bounding box is fixed by the model and does not have "
            "adjustable parameters.")

    @classmethod
    def validate(cls, model, bounding_box):
        """
        Validate a given bounding box sequence against the given model (which
        may be either a subclass of `~astropy.modeling.Model` or an instance
        thereof, so long as the ``.inputs`` attribute is defined.

        Currently this just checks that the bounding_box is either a 2-tuple
        of lower and upper bounds for 1-D models, or an N-tuple of 2-tuples
        for N-D models.

        This also returns a normalized version of the bounding_box input to
        ensure it is always an N-tuple (even for the 1-D case).
        """

        nd = model.n_inputs

        if nd == 1:
            msg = ("Bounding box for {0} model must be a sequence of length "
                   "2 consisting of a lower and upper bound, or a 1-tuple "
                   "containing such a sequence as its sole element.").format(
                       model.name)

            assert (isiterable(bounding_box) and
                        np.shape(bounding_box) in ((2,), (1, 2))), msg

            if len(bounding_box) == 1:
                return cls((tuple(bounding_box[0]),))
            else:
                return cls(tuple(bounding_box))
        else:
            msg = ("Bounding box for {0} model must be a sequence of length "
                   "{1} (the number of model inputs) consisting of pairs of "
                   "lower and upper bounds for those inputs on which to "
                   "evaluate the model.").format(model.name, nd)

            assert (isiterable(bounding_box) and
                        np.shape(bounding_box) == (nd, 2)), msg

            return cls(tuple(bounds) for bounds in bounding_box)


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
    for j in range(min(k, N - k)):
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


def ellipse_extent(a, b, theta):
    """
    Calculates the extent of a box encapsulating a rotated 2D ellipse.

    Parameters
    ----------
    a : float
        Major axis.
    b : float
        Minor axis.
    theta : float
        Rotation angle in radians.

    Returns
    -------
    offsets : tuple
        The absolute value of the offset distances from the ellipse center that
        define its bounding box region, ``(dx, dy)``.

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.modeling.models import Ellipse2D
        from astropy.modeling.utils import ellipse_extent, render_model

        amplitude = 1
        x0 = 50
        y0 = 50
        a = 30
        b = 10
        theta = np.pi/4

        model = Ellipse2D(amplitude, x0, y0, a, b, theta)

        dx, dy = ellipse_extent(a, b, theta)

        limits = [x0 - dx, x0 + dx, y0 - dy, y0 + dy]

        model.bounding_box = limits

        image = render_model(model)

        plt.imshow(image, cmap='binary', interpolation='nearest', alpha=.5,
                  extent = limits)
        plt.show()
    """

    t = np.arctan2(-b * np.tan(theta), a)
    dx = a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)

    t = np.arctan2(b, a * np.tan(theta))
    dy = b * np.sin(t) * np.cos(theta) + a * np.cos(t) * np.sin(theta)

    return np.abs([dx, dy])


def get_inputs_and_params(func):
    """
    Given a callable, determine the input variables and the
    parameters.

    Parameters
    ----------
    func : callable

    Returns
    -------
    inputs, params : tuple
        Each entry is a list of inspect.Parameter objects
    """
    sig = signature(func)

    inputs = []
    params = []
    for param in sig.parameters.values():
        if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
            raise ValueError("Signature must not have *args or **kwargs")
        if param.default == param.empty:
            inputs.append(param)
        else:
            params.append(param)

    return inputs, params

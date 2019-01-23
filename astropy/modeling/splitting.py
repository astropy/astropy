"""
Functions for splitting `~astropy.modeling.CompoundModel` objects into smaller
models.
"""

from functools import wraps
from collections import defaultdict

from astropy.modeling import separable
from astropy.modeling.core import BINARY_OPERATORS, _model_oper, _CompoundModel, _CompoundModelMeta
from astropy.modeling.utils import ExpressionTree

OPERATORS = dict((oper, _model_oper(oper)) for oper in BINARY_OPERATORS)


__all__ = ['remove_input_frame']


def tree_or_compound(func):
    """
    A decorator which processes the first argument of the function to be an
    ExpressionTree.

    If the first argument is a CompoundModel then it will be converted into an
    ExpressionTree and then back to a CompoundModel before being returned.

    """
    @wraps(func)
    def compound_wrapper(tree, *args, **kwargs):
        re_tree = False
        if isinstance(tree, _CompoundModel) or isinstance(tree, _CompoundModelMeta):
            tree = tree._tree
            re_tree = True
        elif not isinstance(tree, ExpressionTree):
            raise TypeError("The tree argument must be either a"
                            " CompoundModel or an ExpressionTree")

        out_trees = func(tree, *args, **kwargs)

        if re_tree:
            return re_model_trees(out_trees)

        return out_trees
    return compound_wrapper


@tree_or_compound
def remove_input_frame(tree, inp, remove_coupled_trees=False):
    """
    Given a tree, remove the smallest subtree needed to remove the input from
    the tree.

    This method traverses the expression tree until it finds a tree with only
    the given input or a tree which has non-separable inputs and takes the
    given input. It then removes this tree and returns a list of all the other
    trees in the original tree.

    Parameters
    ----------
    tree : `astropy.modelling.utils.ExpressionTree`, `astropy.modeling.CompoundModel`
        The tree to analyse the inputs of.

    inp : `str`
        The name of the input to be removed.

    remove_coupled_trees : (optional) `bool`
        If `True` remove the subtree that contains input even if it has other
        non-separable inputs as well. Defaults to `False`.

    Returns
    -------
    new_trees : `list`
        A list of all the trees. Can have the `&` operator applied to
        reconstruct a `CompoundModel`.
    """
    # If the input is not found, noop
    if inp not in tree.inputs:
        return [tree]

    if tree.value != "&":
        # If the input is not a "&" then we have to respect remove_coupled_trees
        sep = tree_is_separable(tree)
        if all(~sep):
            if not remove_coupled_trees:
                return [tree]
        # Otherwise, we know this tree has the input, so we just drop it.
        return []

    # map the input names of this tree to the names of it's children
    inp_map = make_forward_input_map(tree)

    # Map the names of the inputs of tree to the two subtrees of tree
    input_tree_map = make_tree_input_map(tree)

    new_trees = []
    for stree in (tree.left, tree.right):
        # Given the two halves of this tree, and one input we want to remove,
        # let us check if we can drop one half of this tree or if we need to
        # split one of these subtrees.

        # Get the names of the inputs for this subtree
        tree_inputs = tuple(input_tree_map[stree])
        # If the input we want to drop is not in this subtree, we keep it and
        # check the next subtree.
        if inp not in tree_inputs:
            new_trees.append(stree)
            continue

        # If the input is in this subtree, and this subtree only takes one
        # input (which is the input we want to drop) drop this subtree and then
        # check the next one. (In the case where we have dropped this tree we
        # would expect the next iteration (if there is one) to keep the other
        # half of the tree, i.e. meet the above if criteria).
        elif len(tree_inputs) == 1:
            continue

        # Finally if the subtree with the input we want to drop has more than
        # one input, then we need to determine how to handle it based on if the
        # inputs are separable or not.
        else:
            sep = tree_is_separable(stree)
            # If they are separable, then we either keep the whole tree and
            # move on, or we drop the subtree if we are removing whole coupled
            # trees
            if all(~sep):
                if not remove_coupled_trees:
                    new_trees.append(stree)
                    continue
            # Otherwise we split up the subtree by recusing down the tree.
            else:
                new_trees += remove_input_frame(stree, inp_map[inp],
                                                remove_coupled_trees)
    return new_trees


def tree_is_separable(tree):
    """
    Given a tree, convert it to a `CompoundModel` and then return the
    separability of the inputs.
    """
    return separable.is_separable(tree.evaluate(OPERATORS))


def make_tree_input_map(tree):
    """
    Construct a mapping of tree to input names.

    This function iterates over all the inputs of a model and determines which
    side of the tree (left or right) the input corresponds to. It returns a
    mapping of tree to set of all inputs for that side of the tree (which is
    also a tree).

    Parameters
    ----------
    tree : `astropy.modelling.utils.ExpressionTree`
        The tree to analyse the inputs of.

    Returns
    -------
    tree_input_map : `dict`
       A mapping of tree to a set of input names.
    """
    tree_input_map = defaultdict(set)
    for i, inp in enumerate(tree.inputs):

        if tree.isleaf or tree.value != "&":
            # If the tree is a leaf node then the inputs map to the original tree.
            return {tree: set(tree.inputs)}

        # If this input number is less than the number of inputs in the left
        # hand side of the tree then the input maps to the LHS of the tree. If
        # not it maps to the right.
        if i < len(tree.left.inputs):
            tree_input_map[tree.left].add(inp)
        else:
            tree_input_map[tree.right].add(inp)
    return tree_input_map


def make_forward_input_map(tree):
    """
    Given a tree, generate a mapping of inputs to the tree, to inputs of it's
    branches.
    """
    inp_map = {}
    assert tree.value == "&"
    for i, ori_inp in enumerate(tree.inputs):
        if i < len(tree.left.inputs):
            inp_map[ori_inp] = tree.left.inputs[i]
        else:
            inp_map[ori_inp] = tree.right.inputs[i - len(tree.left.inputs)]
    return inp_map


def re_model_trees(trees):
    """
    Given a list of trees return a `CompoundModel` by applying the `&`
    operator.

    Parameters
    ----------
    tree : `list`
        A list of `astropy.modelling.utils.ExpressionTree` objects.

    Returns
    -------
    model : `astropy.modelling.CompoundModel`
        A model.
    """
    left = trees.pop(0).evaluate(OPERATORS)
    for right in trees:
        left = left & right.evaluate(OPERATORS)
    return left

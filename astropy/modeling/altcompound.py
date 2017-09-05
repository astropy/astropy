# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import division, print_function, absolute_import


"""
This module provides an alternate implementation of compound models that
is lighter weight than the default implementation. This is primarily done
to minimize the memory and construction time demands when compound models
become very complex.

The motivation for this arose when JWST NIRSPEC compound models for its
MOS mode were involving nearly 100 constituent model elements and on the
order of 300 parameters. Since the default compound model implementation
maps parameters from the consistuent components (which themselves are often
compound models), this led to many thousands of parameter instances since
the total tends to go with the square of the number of model instances.
This led to enormous memory usage and very long instantiation times.

The current implementation is already very complex and relies heavily
on metaclass tools making it difficult for mere mortals to understand
and modify. In principle it is possible to optimize the current
implementation but it is very hard to see how to avoid embedding operator
expressions in strings to avoid the many instantiations. Besides being
complex and difficult to add, it adds clumsiness to the user interface
(though admittedly, the alternate solution has its own clumsy aspect).
So this led to a completely different approach that seems much simpler
for certain use cases (it does forgo the ability to fit compound models
for example, though adding such capability without affecting the
performance issues isn't out of the question in the following approach)

This is an initial implementation that implements only very basic
functionality. Some of the functionality that the default implementation
provides will never be provided by this implementation; other functionality
is deferred depending on the experience in using this initial implementation.
Note that this does not affect the underlying basic model classes at all.
This model class does not inherit from the basic model class! (So some
solution to making generic tests for models is needed)

Using this alternate version of compound models requires calling a function
in core to make this one the default. If this is used, *is is higly recommended
that the mode be set back to the default at the end of the code constructing
compound models so that other code depending on the default behavior is
not affected!*

As an example of how to do this:

from astropy.modeling.core import set_compound_model
prevcm = set_compound_model('lite')
compound_model = Gaussian1D(1., 0.5, 0.1) + Gaussian1D(2, 0.7, 0.2)
set_compound_model(prevcm) # the default model type is 'regular'

Things currently supported:

- evaluation
- inverse evaluation (if possible or provided)

Things not currently supported (but will be if adopted):

- named access to constituent models (And thus their parameters)
- good string representation and such
- conversion to and from the default compound model scheme
- pickling
- possibly more efficient evaluation as done in the style of the current
  compound models (this implementation walks the tree every time)
- support for units, when available
- and other things I've overlooked at this moment...

Things that will never be supported:

- Compound models of model classes (as opposed to instances)
- Transparent and automatic mapping of parameters of constituent
  models to the compound model level.
- Fitting (the focus of this almost entirely on evaluation)
"""


import operator
from .core import ModelDefinitionError
from .utils import combine_labels


def _alt_model_oper(oper, **kwargs):
    """
    This is an alternate version of compound models intended to use
    much less memory than the default.
    """
    return lambda left, right: _AltCompoundModel(oper, left, right, **kwargs)


class _AltCompoundModel(object):
    '''
    Lightweight compound model implementation: see altcompound.py documentation
    for details.
    '''

    def __init__(self, op, left, right, name=None, inverse=None):
        self.op = op
        self.left = left
        self.right = right
        self.bounding_box = None
        self._has_inverse = False # may be set to True in following code
        if op in ['+', '-', '*', '/', '**']:
            if (left.n_inputs != right.n_inputs) or \
               (left.n_outputs != right.n_outputs):
                raise ModelDefinitionError(
                    'Both operands must match numbers of inputs and outputs')
            else:
                self.n_inputs = left.n_inputs
                self.n_outputs = left.n_outputs
                self.inputs = left.inputs
                self.outputs = left.outputs
        elif op == '&':
            self.n_inputs = left.n_inputs + right.n_inputs
            self.n_outputs = left.n_outputs + right.n_outputs
            self.inputs = combine_labels(left.inputs, right.inputs)
            self.outputs = combine_labels(left.outputs, right.outputs)
            if inverse is None and self.both_inverses_exist():
                self._has_inverse = True
                self._inverse = _AltCompoundModel('&',
                    self.left.inverse, self.right.inverse, inverse=self)
        elif op == '|':
            if left.n_outputs != right.n_inputs:
                print('>>>>>>>>', left.n_outputs, right.n_inputs)
                print(left.name, right.name)
                print(left)
                print(right)
                raise ModelDefinitionError(
                    'left operand number of outputs must match right operand number of inputs')
            self.n_inputs = left.n_inputs
            self.n_outputs = right.n_outputs
            self.inputs = left.inputs
            self.outputs = right.outputs
            if inverse is None and self.both_inverses_exist():
                self._has_inverse = True
                self._inverse = _AltCompoundModel('|',
                    self.right.inverse, self.left.inverse, inverse=self)
        elif op == '%':
            if type(right) != type({}):
                raise ValueError('expecting dictionary for right side of "%" operator')
            else:
                # dict keys must match either possible indices for model on left side,
                # or names for inputs.
                self.n_inputs = left.n_inputs - len(right)
                self.outputs = left.outputs
                newinputs = list(left.inputs)                
                keys = right.keys()
                input_ind = []
                for key in keys:
                    if type(key) == type(0):
                        if key >= left.n_inputs or key < 0:
                            raise ValueError(
                                'substitution key integer value not among possible input choices')
                        else:
                            input_ind.append(key)
                    elif type(key) == type(''):
                        if key not in left.inputs:
                            raise ValueError(
                                'Substitution key string not among possible input choices')
                        # check to see it doesn't match positional specification
                        ind = left.inputs.index(key)
                        if ind in input_ind:
                            raise ValueError("Duplicate specification of same input (index/name)")
                        else:
                            input_ind.append(ind)
                # Remove substituted inputs
                input_ind.sort()
                input_ind.reverse()
                for ind in input_ind:
                    del newinputs[ind]
                self.inputs = tuple(newinputs)

        else:
            raise ModelDefinitionError('Illegal operator: ', self.op)
        if inverse is not None:
            self._inverse = inverse
            self._has_inverse = True
        self.name = name

    def both_inverses_exist(self):
        '''
        if both members of this compound model have inverses return True
        '''
        try:
            linv = self.left.inverse
            rinv = self.right.inverse
        except NotImplementedError:
            return False
        if isinstance(self.left, _AltCompoundModel):
            if  not self.left.has_inverse():
                return False
        if isinstance(self.right, _AltCompoundModel):
            if not self.right.has_inverse():
                return False
        return  True

    def __call__(self, *args, **kw):
        op = self.op
        if op != '%':
            if op != '&':
                leftval = self.left(*args, **kw)
                if op != '|':
                    rightval = self.right(*args, **kw)
            else:
                leftval = self.left(*(args[:self.left.n_inputs]), **kw)
                rightval = self.right(*(args[self.left.n_inputs:]), **kw)
            if op == '+':
                return binary_operation(operator.add, leftval, rightval)
            elif op == '-':
                return binary_operation(operator.sub, leftval, rightval)
            elif op == '*':
                return binary_operation(operator.mul, leftval, rightval)
            elif op == '/':
                return binary_operation(operator.truediv, leftval, rightval)
            elif op == '**':
                return binary_operation(operator.pow, leftval, rightval)
            elif op == '&':
                if not isinstance(leftval, tuple):
                    leftval = (leftval,)
                if not isinstance(rightval, tuple):
                    rightval = (rightval,)
                return leftval + rightval
            elif op == '|':
                if isinstance(leftval, tuple):
                    return self.right(*leftval, **kw)
                else:
                    return self.right(leftval, **kw)
        elif op == '%':
            subs = self.right
            newargs = list(args)
            subinds = []
            subvals = []
            for key in subs.keys():
                if type(key) == type(0):
                    subinds.append(key)
                elif type(key) == type(''):
                    ind = self.left.inputs.index(key)
                    subinds.append(ind)
                subvals.append(subs[key])
            # Turn inputs specified in kw into positional indices.
            # Names for compound inputs do not propagate to sub models.
            kwind = []
            kwval = []
            for kwkey in list(kw.keys()):
                if kwkey in self.inputs:
                    ind = self.inputs.index(kwkey)
                    if ind < len(args):
                        raise ValueError("Keyword argument duplicates positional value supplied")
                    kwind.append(ind)
                    kwval.append(kw[kwkey])
                    del kw[kwkey]
            # Build new argument list
            # Append keyword specified args first
            if kwind:
                kwargs = list(zip(kwind, kwval))
                kwargs.sort()
                kwindsorted, kwvalsorted = list(zip(*kwargs))
                newargs = newargs + list(kwvalsorted)
            if subinds:
                subargs = list(zip(subinds, subvals))
                subargs.sort()
                subindsorted, subvalsorted = list(zip(*subargs))
            # The substitutions must be inserted in order
            for ind, val in subargs:
                newargs.insert(ind, val)
            return self.left(*newargs, **kw)
        else:
            raise ModelDefinitionError('unrecognized operator')

    def __getitem__(self, index):
        tdict = {}
        leaflist = []
        make_subtree_dict(self, '', tdict, leaflist)
        if isinstance(index, slice):
            if index.step:
                raise ValueError('Steps in slices not supported for compound models')
            # Following won't work for negative indices
            if index.start:
                start = index.start
            else:
                start = 0
            if index.stop:
                stop = index.stop
            else:
                stop = len(leaflist) - 1
            if start < 0:
                start = len(leaflist) + start
            if stop < 0:
                stop = len(leaflist) +  stop
            # now search for matching node:
            for key in tdict:
                node, leftind, rightind = tdict[key]
                if leftind == start and rightind == stop:
                    return node
            raise IndexError("No appropriate subtree matches slice")
        elif isinstance(index, type(0)):
            return leaflist[index]
        else:
            raise TypeError('index must be integer or slice')

    def rename(self, name):
        self.name = name
        return self

    @property
    def isleaf(self):
        return False

    def has_inverse(self):
        return self._has_inverse

    @property
    def bounding_box(self):
        return self._bounding_box

    @bounding_box.setter
    def bounding_box(self, bounding_box):
        self._bounding_box = bounding_box

    @property
    def inverse(self):
        '''
        '''
        if self.has_inverse():
            return self._inverse
        else:
            raise NotImplementedError("Inverse function not provided")

    @inverse.setter
    def inverse(self, invmodel):
        if not (isinstance(invmodel, Model) or isinstance(invmodel, _AltCompoundModel)):
            raise ValueError("Attempt to assign non model to inverse")
        self._has_inverse = True 
        self._inverse = invmodel

    __add__ =     _alt_model_oper('+')
    __sub__ =     _alt_model_oper('-')
    __mul__ =     _alt_model_oper('*')
    __truediv__ = _alt_model_oper('/')
    __pow__ =     _alt_model_oper('**')
    __or__ =      _alt_model_oper('|')
    __and__ =     _alt_model_oper('&')
    __mod__ =     _alt_model_oper('%')


def binary_operation(binoperator, left, right):
    '''
    Perform binary operation. Operands may be matching tuples of operands.
    '''
    if isinstance(left, tuple) and isinstance(right, tuple):
        return tuple([binoperator(item[0], item[1]) for item in zip(left, right)])
    else:
        return binoperator(left, right)

def tree_to_list(tree):
    '''
    Walk a tree and return a list of the leaves in order of traversal.
    '''
    if not hasattr(tree, 'isleaf'):
        return [tree]
    else:
        return tree_to_list(tree.left, get_oper) + tree_to_list(tree.right, get_oper)

def make_subtree_dict(tree, nodepath, tdict, leaflist):
    '''
    Traverse a tree noting each node by a key that indicates all the left/right choices
    necessary to reach that node. Each key will reference a tuple that contains:

    - reference to the compound model for that node.
    - left most index contained within that subtree
       (relative to all indices for the whole tree)
    - right most index contained within that subtree
    '''
    # if this is a leaf, just append it to the leaflist
    if not hasattr(tree, 'isleaf'):
        leaflist.append(tree)
    else:
        leftmostind = len(leaflist)
        make_subtree_dict(tree.left, nodepath+'l', tdict, leaflist)
        make_subtree_dict(tree.right, nodepath+'r', tdict, leaflist)
        rightmostind = len(leaflist)-1
        tdict[nodepath] = (tree, leftmostind, rightmostind)



from .core import Model

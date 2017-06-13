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
        self._has_inverse = False # may be set to True in following code
        if op in ['+', '-', '*', '/', '**']:
            if (left.n_inputs != right.n_inputs) or \
               (left.n_outputs != right.n_outputs):
                raise ModelDefinitionError(
                    'Both operands must match numbers of inputs and outputs')
            else:
                self.n_inputs = left.n_inputs
                self.n_outputs = left.n_outputs
        elif op == '&':
            self.n_inputs = left.n_inputs + right.n_inputs
            self.n_outputs = left.n_outputs + right.n_outputs
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
            if inverse is None and self.both_inverses_exist():
                self._has_inverse = True
                self._inverse = _AltCompoundModel('|',
                    self.right.inverse, self.left.inverse, inverse=self)
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
        else:
            raise ModelDefinitionError('unrecognized operator')

    def __getitem__(self, index):
        if isinstance(index, slice):
            model_list = tree_to_list(self, get_oper=True)
            # now figure out which list indicies correspond to the requested ones
            if index.step:
                raise ValueError('Steps in slices not supported for compound models')
            # Following won't work for negative indices
            if index.start:
                lstart = index.start * 2
            else:
                lstart = None
            if index.stop:
                lstop = indes.stop * 2
            else:
                lstop = None
            sliced_model_list = model_list[slice(lstart, lstop)]
            # Now create new compound model from sliced list
            if len(sliced_model_list) < 1:
                raise ValueError('must have at least one element in slice') 
            if len(sliced_model_list) == 1:
                return sliced_model_list[0]
            else:
                temp = _AltCompoundModel(sliced_model_list[1],
                                         sliced_model_list[0],
                                         sliced_model_list[2])
                sliced_model_list = sliced_model_list[3:]
                while sliced_model_list:
                    temp = _AltCompoundModel( 
                        sliced_model_list[0],
                        temp,
                        sliced_model_list[1])
                    sliced_model_list = sliced_model_list[2:]
                return temp
        elif isinstance(index, type(0)):
            model_list = tree_to_list(self)
        else:
            raise TypeError('index must be integer')
        return model_list[index]



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
        raise NotImplementedError("no bounding box defined")

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
        self._inverse = invmodel

    __add__ =     _alt_model_oper('+')
    __sub__ =     _alt_model_oper('-')
    __mul__ =     _alt_model_oper('*')
    __truediv__ = _alt_model_oper('/')
    __pow__ =     _alt_model_oper('**')
    __or__ =      _alt_model_oper('|')
    __and__ =     _alt_model_oper('&')


def binary_operation(binoperator, left, right):
    '''
    Perform binary operation. Operands may be matching tuples of operands.
    '''
    if isinstance(left, tuple) and isinstance(right, tuple):
        return tuple([binoperator(item[0], item[1]) for item in zip(left, right)])
    else:
        return binoperator(left, right)

def tree_to_list(tree, get_oper=False):
    '''
    Walk a tree and return a list of the leaves in order of traversal.
    '''
    if not hasattr(tree,'isleaf'):
        return [tree]
    else:
        if not get_oper:
            return tree_to_list(tree.left, get_oper) + tree_to_list(tree.right, get_oper)
        else:
            return tree_to_list(tree.left, get_oper) + [tree.op] + tree_to_list(tree.right, get_oper)
 


from .core import Model

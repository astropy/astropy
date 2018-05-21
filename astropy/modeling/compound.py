# Licensed under a 3-clause BSD style license - see LICENSE.rst

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
from collections import deque, OrderedDict
import numpy as np
from . import core
from ..units import Quantity, UnitsError
from ..utils import indent
from .utils import combine_labels, _BoundingBox
from .parameters import _tofloat, Parameter, InputParameterError


def _model_oper(oper, **kwargs):
    """
    This is an alternate version of compound models intended to use
    much less memory than the default.
    """
    return lambda left, right: CompoundModel(oper, left, right, **kwargs)


class CompoundModel(object):
    '''
    Lightweight compound model implementation: see altcompound.py documentation
    for details.
    '''

    def __init__(self, op, left, right, name=None, inverse=None):
        self.__dict__['param_names'] = None
        self._n_submodels = None
        self.op = op
        self.left = left
        self.right = right
        self._bounding_box = None
        self._user_bounding_box = None
        self._leaflist = None
        self._parameters = None
        self._parameters_ = None
        self._param_metrics = None
        self._has_inverse = False  # may be set to True in following code
        if op != '%' and len(left) != len(right):
            raise ValueError(
                'Both operands must have equal values for n_models')
        else:
            self._n_models = len(left)
        if op in ['+', '-', '*', '/', '**']:
            if (left.n_inputs != right.n_inputs) or \
               (left.n_outputs != right.n_outputs):
                raise core.ModelDefinitionError(
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
                self._inverse = CompoundModel('&',
                                              self.left.inverse,
                                              self.right.inverse,
                                              inverse=self)
        elif op == '|':
            if left.n_outputs != right.n_inputs:
                print('>>>>>>>>', left.n_outputs, right.n_inputs)
                print(left.name, right.name)
                print(left)
                print(right)
                raise core.ModelDefinitionError(
                    'left operand number of outputs must'
                    'match right operand number of inputs')
            self.n_inputs = left.n_inputs
            self.n_outputs = right.n_outputs
            self.inputs = left.inputs
            self.outputs = right.outputs
            if inverse is None and self.both_inverses_exist():
                self._has_inverse = True
                self._inverse = CompoundModel('|',
                                              self.right.inverse,
                                              self.left.inverse,
                                              inverse=self)
        elif op == '%':
            if not isinstance(right, dict):
                raise ValueError('expecting dictionary for right side of "%"'
                                 ' operator')
            else:
                # Dict keys must match either possible indices
                # for model on left side,
                # or names for inputs.
                self.n_inputs = left.n_inputs - len(right)
                self.outputs = left.outputs
                self.n_outputs = left.n_outputs
                newinputs = list(left.inputs)
                keys = right.keys()
                input_ind = []
                for key in keys:
                    if isinstance(key, int):
                        if key >= left.n_inputs or key < 0:
                            raise ValueError(
                                'substitution key integer value '
                                'not among possible input choices')
                        else:
                            input_ind.append(key)
                    elif isinstance(key, str):
                        if key not in left.inputs:
                            raise ValueError(
                                'Substitution key string not among possible '
                                'input choices')
                        # Check to see it doesn't match positional
                        # specification.
                        ind = left.inputs.index(key)
                        if ind in input_ind:
                            raise ValueError("Duplicate specification of "
                                             "same input (index/name)")
                        else:
                            input_ind.append(ind)
                # Remove substituted inputs
                input_ind.sort()
                input_ind.reverse()
                for ind in input_ind:
                    del newinputs[ind]
                self.inputs = tuple(newinputs)

        else:
            raise core.ModelDefinitionError('Illegal operator: ', self.op)
        if inverse is not None:
            self._inverse = inverse
            self._has_inverse = True
        self.name = name
        self.fittable = False
        self.linear = False
        self.eqcons = False
        self.ineqcons = False

    def __len__(self):
        return self._n_models

    @property
    def n_submodels(self):
        if self._n_submodels is None:
            tdict = {}
            leaflist = []
            make_subtree_dict(self, '', tdict, leaflist)
            self._n_submodels = len(leaflist)
        return self._n_submodels

    def both_inverses_exist(self):
        '''
        if both members of this compound model have inverses return True
        '''
        try:
            linv = self.left.inverse
            rinv = self.right.inverse
        except NotImplementedError:
            return False
        if isinstance(self.left, CompoundModel):
            if not self.left.has_inverse():
                return False
        if isinstance(self.right, CompoundModel):
            if not self.right.has_inverse():
                return False
        return True

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
                if isinstance(key, int):
                    subinds.append(key)
                elif isinstance(key, str):
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
                        raise ValueError("Keyword argument duplicates "
                                        "positional value supplied")
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
            raise core.ModelDefinitionError('unrecognized operator')

    def _make_leaflist(self):
        tdict = {}
        leaflist = []
        make_subtree_dict(self, '', tdict, leaflist)
        self._leaflist = leaflist
        self._tdict = tdict

    def __getitem__(self, index):
        if self._leaflist is None:
            self._make_leaflist()
        leaflist = self._leaflist
        tdict = self._tdict
        if isinstance(index, slice):
            if index.step:
                raise ValueError('Steps in slices not supported '
                                 'for compound models')
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
                stop = len(leaflist) + stop
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

    def __setattr__(self, attr, value):
        # TODO should eliminate the duplicate code here with that in core
        if self.param_names is not None and attr in self.param_names:
            param = self.__dict__[attr]
            value = _tofloat(value)
            if param.unit is None:
                if isinstance(value, Quantity):
                    param._unit = value.unit
                    param.value = value.value
                else:
                    param.value = value
            else:
                if not isinstance(value, Quantity):
                    raise UnitsError("The '{0}' parameter should be given "
                                     "as a Quantity because it was originally "
                                     "initialized as a "
                                     "Quantity".format(param.name))
                else:
                    param._unit = value.unit
                    param.value = value.value
        else:
            super().__setattr__(attr, value)

    def traverse_postorder(self):
        stack = deque([self])
        stacked = deque([])
        while stack:
            node = stack[-1]
            if not isinstance(node, CompoundModel):
                yield stack.pop()
            elif node not in stacked:
                stacked.append(node)
                stack.append(node.right)
                stack.append(node.left)
            else:
                yield stack.pop()

    def _format_expression(self, format_leaf=None):
        leaf_idx = 0
        operands = deque()

        if format_leaf is None:
            format_leaf = lambda i, l: '[{0}]'.format(i)

        for node in self.traverse_postorder():
            if not isinstance(node, CompoundModel):
                operands.append(format_leaf(leaf_idx, node))
                leaf_idx += 1
                continue

            oper_order = OPERATOR_PRECEDENCE[node.op]
            right = operands.pop()
            left = operands.pop()

            if isinstance(node, CompoundModel):
                if (isinstance(node.left, CompoundModel) and
                        OPERATOR_PRECEDENCE[node.left.op] < oper_order):
                    left = '({0})'.format(left)
                if (isinstance(node.right, CompoundModel) and
                        OPERATOR_PRECEDENCE[node.right.op] < oper_order):
                    right = '({0})'.format(right)

            operands.append(' '.join((left, node.op, right)))

        return ''.join(operands)

    def _format_repr(self, keywords=[]):
        """
        Internal implementation of ``__repr__``.

        This is separated out for ease of use by subclasses that wish to
        override the default ``__repr__`` while keeping the same basic
        formatting.
        """

        # For the sake of familiarity start the output with the standard class
        # __repr__

        parts = []
        try:
            default_keywords = [
                ('Name', 'CompoundModel'),
                ('Inputs', self.inputs),
                ('Outputs', self.outputs),
            ]

            if self.param_names:
                default_keywords.append(('Fittable parameters',
                                         self.param_names))

            for keyword, value in default_keywords + keywords:
                if value is not None:
                    parts.append('{0}: {1}'.format(keyword, value))

            return '\n'.join(parts)
        except Exception:
            # If any of the above formatting fails fall back on the basic repr
            # (this is particularly useful in debugging)
            return parts[0]

    def _format_components(self):
        return '\n\n'.join('[{0}]: {1!r}'.format(idx, m)
                                 for idx, m in enumerate(self._leaflist))

    def __repr__(self):
        if self._parameters_ is None:
            return ("CompoundModel--call map_parameters()"+
                    " method to get full __repr__")
        expression = self._format_expression()
        components = self._format_components()
        keywords = [
            ('Expression', expression),
            ('Components', '\n' + indent(components))
        ]

        return self._format_repr(keywords=keywords)

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

    @property
    def parameters(self):
        """
        A flattened array of all parameter values in all parameter sets.

        Fittable parameters maintain this list and fitters modify it.
        """

        # Currently the sequence of a model's parameters must be contiguous
        # within the _parameters array (which may be a view of a larger array,
        # for example when taking a sub-expression of a compound model), so
        # the assumption here is reliable:
        if self.param_names is None:
            raise RuntimeError("Compound model parameter interface is not "
                               "supported\n"
                               "until the .map_parameters() method is called.")

        self._parameters_to_array()
        start = self._param_metrics[self.param_names[0]]['slice'].start
        stop = self._param_metrics[self.param_names[-1]]['slice'].stop

        return self._parameters[start:stop]

    @parameters.setter
    def parameters(self, value):
        """
        Assigning to this attribute updates the parameters array rather than
        replacing it.
        """

        if self.param_names is None:
            raise RuntimeError("Compound model parameter interface is not "
                               "supported\n"
                               "until the .map_parameters() method is called.")

        start = self._param_metrics[self.param_names[0]]['slice'].start
        stop = self._param_metrics[self.param_names[-1]]['slice'].stop

        try:
            value = np.array(value).flatten()
            self._parameters[start:stop] = value
        except ValueError as e:
            raise InputParameterError(
                "Input parameter values not compatible with the model "
                "parameters array: {0}".format(e))
        self._array_to_parameters()

    @inverse.setter
    def inverse(self, invmodel):
        if not (isinstance(invmodel, core.Model) or
                isinstance(invmodel, CompoundModel)):
            raise ValueError("Attempt to assign non model to inverse")
        self._has_inverse = True
        self._inverse = invmodel

    __add__ =     _model_oper('+')
    __sub__ =     _model_oper('-')
    __mul__ =     _model_oper('*')
    __truediv__ = _model_oper('/')
    __pow__ =     _model_oper('**')
    __or__ =      _model_oper('|')
    __and__ =     _model_oper('&')
    __mod__ =     _model_oper('%')

    def map_parameters(self, namestyle=None):
        """
        Map all the constituent model parameters to the compound object,
        renaming as necessary by appending a suffix number.

        This can be an expensive operation, particularly for a complex
        expression tree.

        All the corresponding parameter attributes are created that one
        expects for the Model class.

        The parameter objects that the attributes point to are the same
        objects as in the constiutent models. Changes made to parameter
        values to either are seen by both.

        Prior to calling this, none of the associated attributes will
        exist. This method must be called to make the model usable by
        fitting engines.

        If oldnames=True, then parameters are named as in the original
        implementation of compound models.
        """
        if self._parameters is not None:
            # do nothing
            return
        if self._leaflist is None:
            self._make_leaflist()
        self._parameters_ = OrderedDict()
        self.param_names = []
        for lindex, leaf in enumerate(self._leaflist):
            for param_name in leaf.param_names:
                param = getattr(leaf, param_name)
                new_param_name = "{}_{}".format(param_name, lindex)
                self.__dict__[new_param_name] = param
                self._parameters_[new_param_name] = param
                self.param_names.append(new_param_name)
        self._param_metrics = {}
        self._initialize_slices()
        self._initialize_constraints()

    def _initialize_slices(self):
        # TODO eliminate redundant code with core here and next two methods
        param_metrics = self._param_metrics
        total_size = 0

        for name in self.param_names:
            param = getattr(self, name)
            value = param.value
            param_size = np.size(value)
            param_shape = np.shape(value)
            param_slice = slice(total_size, total_size + param_size)
            param_metrics[name] = {}
            param_metrics[name]['slice'] = param_slice
            param_metrics[name]['shape'] = param_shape
            total_size += param_size
        self._parameters = np.empty(total_size, dtype=np.float64)

    def _initialize_constraints(self):

        self._constraints = {}
        for constraint in Parameter.constraints:
            self._constraints[constraint] = {}
            # Update with default parameter constraints
            for param_name in self.param_names:
                param = getattr(self, param_name)
                # Parameters don't have all constraint types
                value = getattr(param, constraint)
                if value is not None:
                    self._constraints[constraint][param_name] = value

    def _parameters_to_array(self):
        # Now set the parameter values (this will also fill
        # self._parameters)
        # TODO: This is a bit ugly, but easier to deal with than how this was
        # done previously.  There's still lots of opportunity for refactoring
        # though, in particular once we move the _get/set_model_value methods
        # out of Parameter and into Model (renaming them
        # _get/set_parameter_value)
        param_metrics = self._param_metrics
        for name in self.param_names:
            param = getattr(self, name)
            value = param.value
            if not isinstance(value, np.ndarray):
                value = np.array([value])
            self._parameters[param_metrics[name]['slice']] = value.ravel()

    def _array_to_parameters(self):
        param_metrics = self._param_metrics
        for name in self.param_names:
            param = getattr(self, name)
            param.value = self._parameters[param_metrics[name]['slice']]
            param.shape = param_metrics[name]['shape']
    
    @property
    def _has_units(self):
        # Returns True if any of the parameters have units
        for param in self.param_names:
            if getattr(self, param).unit is not None:
                return True
        else:
            return False

    @property
    def bounding_box(self):
        r"""
        A `tuple` of length `n_inputs` defining the bounding box limits, or
        `None` for no bounding box.

        The default limits are given by a ``bounding_box`` property or method
        defined in the class body of a specific model.  If not defined then
        this property just raises `NotImplementedError` by default (but may be
        assigned a custom value by a user).  ``bounding_box`` can be set
        manually to an array-like object of shape ``(model.n_inputs, 2)``. For
        further usage, see :ref:`bounding-boxes`

        The limits are ordered according to the `numpy` indexing
        convention, and are the reverse of the model input order,
        e.g. for inputs ``('x', 'y', 'z')``, ``bounding_box`` is defined:

        * for 1D: ``(x_low, x_high)``
        * for 2D: ``((y_low, y_high), (x_low, x_high))``
        * for 3D: ``((z_low, z_high), (y_low, y_high), (x_low, x_high))``

        Examples
        --------

        Setting the ``bounding_box`` limits for a 1D and 2D model:

        >>> from astropy.modeling.models import Gaussian1D, Gaussian2D
        >>> model_1d = Gaussian1D()
        >>> model_2d = Gaussian2D(x_stddev=1, y_stddev=1)
        >>> model_1d.bounding_box = (-5, 5)
        >>> model_2d.bounding_box = ((-6, 6), (-5, 5))

        Setting the bounding_box limits for a user-defined 3D `custom_model`:

        >>> from astropy.modeling.models import custom_model
        >>> def const3d(x, y, z, amp=1):
        ...    return amp
        ...
        >>> Const3D = custom_model(const3d)
        >>> model_3d = Const3D()
        >>> model_3d.bounding_box = ((-6, 6), (-5, 5), (-4, 4))

        To reset ``bounding_box`` to its default limits just delete the
        user-defined value--this will reset it back to the default defined
        on the class:

        >>> del model_1d.bounding_box

        To disable the bounding box entirely (including the default),
        set ``bounding_box`` to `None`:

        >>> model_1d.bounding_box = None
        >>> model_1d.bounding_box  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "astropy\modeling\core.py", line 980, in bounding_box
            "No bounding box is defined for this model (note: the "
        NotImplementedError: No bounding box is defined for this model (note:
        the bounding box was explicitly disabled for this model; use `del
        model.bounding_box` to restore the default bounding box, if one is
        defined for this model).
        """

        if self._user_bounding_box is not None:
            if self._user_bounding_box is NotImplemented:
                raise NotImplementedError(
                    "No bounding box is defined for this model (note: the "
                    "bounding box was explicitly disabled for this model; "
                    "use `del model.bounding_box` to restore the default "
                    "bounding box, if one is defined for this model).")
            return self._user_bounding_box
        elif self._bounding_box is None:
            raise NotImplementedError(
                    "No bounding box is defined for this model.")
        elif isinstance(self._bounding_box, _BoundingBox):
            # This typically implies a hard-coded bounding box.  This will
            # probably be rare, but it is an option
            return self._bounding_box
        elif isinstance(self._bounding_box, types.MethodType):
            return self._bounding_box()
        else:
            # The only other allowed possibility is that it's a _BoundingBox
            # subclass, so we call it with its default arguments and return an
            # instance of it (that can be called to recompute the bounding box
            # with any optional parameters)
            # (In other words, in this case self._bounding_box is a *class*)
            bounding_box = self._bounding_box((), _model=self)()
            return self._bounding_box(bounding_box, _model=self)

    @bounding_box.setter
    def bounding_box(self, bounding_box):
        """
        Assigns the bounding box limits.
        """

        if bounding_box is None:
            cls = None
            # We use this to explicitly set an unimplemented bounding box (as
            # opposed to no user bounding box defined)
            bounding_box = NotImplemented
        elif (isinstance(self._bounding_box, type) and
                issubclass(self._bounding_box, _BoundingBox)):
            cls = self._bounding_box
        else:
            cls = _BoundingBox

        if cls is not None:
            try:
                bounding_box = cls.validate(self, bounding_box)
            except ValueError as exc:
                raise ValueError(exc.args[0])

        self._user_bounding_box = bounding_box

    @bounding_box.deleter
    def bounding_box(self):
        self._user_bounding_box = None

    @property
    def has_user_bounding_box(self):
        """
        A flag indicating whether or not a custom bounding_box has been
        assigned to this model by a user, via assignment to
        ``model.bounding_box``.
        """

        return self._user_bounding_box is not None

    @property
    def fixed(self):
        """
        A `dict` mapping parameter names to their fixed constraint.
        """

        return self._constraints['fixed']

    @property
    def tied(self):
        """
        A `dict` mapping parameter names to their tied constraint.
        """

        return self._constraints['tied']

    @property
    def bounds(self):
        """
        A `dict` mapping parameter names to their upper and lower bounds as
        ``(min, max)`` tuples.
        """

        return self._constraints['bounds']

    def render(self, out=None, coords=None):
        """
        Evaluate a model at fixed positions, respecting the ``bounding_box``.

        The key difference relative to evaluating the model directly is that
        this method is limited to a bounding box if the `Model.bounding_box`
        attribute is set.

        Parameters
        ----------
        out : `numpy.ndarray`, optional
            An array that the evaluated model will be added to.  If this is not
            given (or given as ``None``), a new array will be created.
        coords : array-like, optional
            An array to be used to translate from the model's input coordinates
            to the ``out`` array. It should have the property that
            ``self(coords)`` yields the same shape as ``out``.  If ``out`` is
            not specified, ``coords`` will be used to determine the shape of
            the returned array. If this is not provided (or None), the model
            will be evaluated on a grid determined by `Model.bounding_box`.

        Returns
        -------
        out : `numpy.ndarray`
            The model added to ``out`` if  ``out`` is not ``None``, or else a
            new array from evaluating the model over ``coords``.
            If ``out`` and ``coords`` are both `None`, the returned array is
            limited to the `Model.bounding_box` limits. If
            `Model.bounding_box` is `None`, ``arr`` or ``coords`` must be
            passed.

        Raises
        ------
        ValueError
            If ``coords`` are not given and the the `Model.bounding_box` of
            this model is not set.

        Examples
        --------
        :ref:`bounding-boxes`
        """

        try:
            bbox = self.bounding_box
        except NotImplementedError:
            bbox = None

        ndim = self.n_inputs

        if (coords is None) and (out is None) and (bbox is None):
            raise ValueError('If no bounding_box is set, '
                             'coords or out must be input.')

        # for consistent indexing
        if ndim == 1:
            if coords is not None:
                coords = [coords]
            if bbox is not None:
                bbox = [bbox]

        if coords is not None:
            coords = np.asanyarray(coords, dtype=float)
            # Check dimensions match out and model
            assert len(coords) == ndim
            if out is not None:
                if coords[0].shape != out.shape:
                    raise ValueError('inconsistent shape of the output.')
            else:
                out = np.zeros(coords[0].shape)

        if out is not None:
            out = np.asanyarray(out, dtype=float)
            if out.ndim != ndim:
                raise ValueError('the array and model must have the same '
                                 'number of dimensions.')

        if bbox is not None:
            # Assures position is at center pixel, important when usin
            # add_array.
            pd = np.array([(np.mean(bb), np.ceil((bb[1] - bb[0]) / 2))
                           for bb in bbox]).astype(int).T
            pos, delta = pd

            if coords is not None:
                sub_shape = tuple(delta * 2 + 1)
                sub_coords = np.array([extract_array(c, sub_shape, pos)
                                       for c in coords])
            else:
                limits = [slice(p - d, p + d + 1, 1) for p, d in pd.T]
                sub_coords = np.mgrid[limits]

            sub_coords = sub_coords[::-1]

            if out is None:
                out = self(*sub_coords)
            else:
                try:
                    out = add_array(out, self(*sub_coords), pos)
                except ValueError:
                    raise ValueError(
                        'The `bounding_box` is larger than the input out in '
                        'one or more dimensions. Set '
                        '`model.bounding_box = None`.')
        else:
            if coords is None:
                im_shape = out.shape
                limits = [slice(i) for i in im_shape]
                coords = np.mgrid[limits]

            coords = coords[::-1]

            out += self(*coords)

        return out


def binary_operation(binoperator, left, right):
    '''
    Perform binary operation. Operands may be matching tuples of operands.
    '''
    if isinstance(left, tuple) and isinstance(right, tuple):
        return tuple([binoperator(item[0], item[1])
                      for item in zip(left, right)])
    else:
        return binoperator(left, right)


def tree_to_list(tree):
    '''
    Walk a tree and return a list of the leaves in order of traversal.
    '''
    if not hasattr(tree, 'isleaf'):
        return [tree]
    else:
        return tree_to_list(tree.left, get_oper) + \
               tree_to_list(tree.right, get_oper)


def make_subtree_dict(tree, nodepath, tdict, leaflist):
    '''
    Traverse a tree noting each node by a key that indicates all the
    left/right choices necessary to reach that node. Each key will
    reference a tuple that contains:

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

_ORDER_OF_OPERATORS = [('|',), ('&',), ('+', '-'), ('*', '/'), ('**',)]
OPERATOR_PRECEDENCE = {}
for idx, ops in enumerate(_ORDER_OF_OPERATORS):
    for op in ops:
        OPERATOR_PRECEDENCE[op] = idx
del idx, op, ops

try:
    import asdf_compound
except ImportError:
    pass

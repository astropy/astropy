# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module defines base classes for all models.  The base class of all
models is `~astropy.modeling.Model`. `~astropy.modeling.FittableModel` is
the base class for all fittable models. Fittable models can be linear or
nonlinear in a regression analysis sense.

All models provide a `__call__` method which performs the transformation in
a purely mathematical way, i.e. the models are unitless.  Model instances can
represent either a single model, or a "model set" representing multiple copies
of the same type of model, but with potentially different values of the
parameters in each model making up the set.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import abc
import copy
import inspect
import functools
import operator
import sys
import types
import warnings

from collections import defaultdict
from itertools import chain, islice

import numpy as np

from ..utils import indent, isinstancemethod, metadata
from ..extern import six
from ..extern.six.moves import copyreg
from ..table import Table
from ..utils import (deprecated, sharedmethod, find_current_module,
                     InheritDocstrings, OrderedDescriptorContainer)
from ..utils.codegen import make_function_with_signature
from ..utils.compat import ignored
from ..utils.compat.funcsigs import signature
from ..utils.compat.odict import OrderedDict
from ..utils.exceptions import AstropyDeprecationWarning
from .utils import (array_repr_oneline, check_broadcast, combine_labels,
                    make_binary_operator_eval, ExpressionTree,
                    IncompatibleShapeError, AliasDict, get_inputs_and_params,
                    _BoundingBox)
from ..nddata.utils import add_array, extract_array

from .parameters import Parameter, InputParameterError


__all__ = ['Model', 'FittableModel', 'Fittable1DModel', 'Fittable2DModel',
           'custom_model', 'ModelDefinitionError']


class ModelDefinitionError(TypeError):
    """Used for incorrect models definitions"""


def _model_oper(oper, **kwargs):
    """
    Returns a function that evaluates a given Python arithmetic operator
    between two models.  The operator should be given as a string, like ``'+'``
    or ``'**'``.

    Any additional keyword arguments passed in are passed to
    `_CompoundModelMeta._from_operator`.
    """

    # Note: Originally this used functools.partial, but that won't work when
    # used in the class definition of _CompoundModelMeta since
    # _CompoundModelMeta has not been defined yet.

    # Perform an arithmetic operation on two models.
    return lambda left, right: _CompoundModelMeta._from_operator(oper,
            left, right, **kwargs)


class _ModelMeta(OrderedDescriptorContainer, InheritDocstrings, abc.ABCMeta):
    """
    Metaclass for Model.

    Currently just handles auto-generating the param_names list based on
    Parameter descriptors declared at the class-level of Model subclasses.
    """

    registry = set()
    """
    A registry of all known concrete (non-abstract) Model subclasses.
    """

    _is_dynamic = False
    """
    This flag signifies whether this class was created in the "normal" way,
    with a class statement in the body of a module, as opposed to a call to
    `type` or some other metaclass constructor, such that the resulting class
    does not belong to a specific module.  This is important for pickling of
    dynamic classes.

    This flag is always forced to False for new classes, so code that creates
    dynamic classes should manually set it to True on those classes when
    creating them.
    """

    # Default empty dict for _parameters_, which will be empty on model
    # classes that don't have any Parameters
    _parameters_ = OrderedDict()

    def __new__(mcls, name, bases, members):
        # See the docstring for _is_dynamic above
        if '_is_dynamic' not in members:
            members['_is_dynamic'] = mcls._is_dynamic

        return super(_ModelMeta, mcls).__new__(mcls, name, bases, members)

    def __init__(cls, name, bases, members):
        # Make sure OrderedDescriptorContainer gets to run before doing
        # anything else
        super(_ModelMeta, cls).__init__(name, bases, members)

        if cls._parameters_:
            if hasattr(cls, '_param_names'):
                # Slight kludge to support compound models, where
                # cls.param_names is a property; could be improved with a
                # little refactoring but fine for now
                cls._param_names = tuple(cls._parameters_)
            else:
                cls.param_names = tuple(cls._parameters_)

        cls._create_inverse_property(members)
        cls._create_bounding_box_property(members)
        cls._handle_backwards_compat(name, members)
        cls._handle_special_methods(members)

        if not inspect.isabstract(cls) and not name.startswith('_'):
            cls.registry.add(cls)

    def __repr__(cls):
        """
        Custom repr for Model subclasses.
        """

        return cls._format_cls_repr()

    def _repr_pretty_(cls, p, cycle):
        """
        Repr for IPython's pretty printer.

        By default IPython "pretty prints" classes, so we need to implement
        this so that IPython displays the custom repr for Models.
        """

        p.text(repr(cls))

    def __reduce__(cls):
        if not cls._is_dynamic:
            # Just return a string specifying where the class can be imported
            # from
            return cls.__name__
        else:
            members = dict(cls.__dict__)
            # Delete any ABC-related attributes--these will be restored when
            # the class is reconstructed:
            for key in list(members):
                if key.startswith('_abc_'):
                    del members[key]

            # Delete custom __init__ and __call__ if they exist:
            for key in ('__init__', '__call__'):
                if key in members:
                    del members[key]

            return (type(cls), (cls.__name__, cls.__bases__, members))

    @property
    def name(cls):
        """
        The name of this model class--equivalent to ``cls.__name__``.

        This attribute is provided for symmetry with the `Model.name` attribute
        of model instances.
        """

        return cls.__name__

    @property
    def n_inputs(cls):
        return len(cls.inputs)

    @property
    def n_outputs(cls):
        return len(cls.outputs)

    @property
    def _is_concrete(cls):
        """
        A class-level property that determines whether the class is a concrete
        implementation of a Model--i.e. it is not some abstract base class or
        internal implementation detail (i.e. begins with '_').
        """
        return not (cls.__name__.startswith('_') or inspect.isabstract(cls))

    def rename(cls, name):
        """
        Creates a copy of this model class with a new name.

        The new class is technically a subclass of the original class, so that
        instance and type checks will still work.  For example::

            >>> from astropy.modeling.models import Rotation2D
            >>> SkyRotation = Rotation2D.rename('SkyRotation')
            >>> SkyRotation
            <class '__main__.SkyRotation'>
            Name: SkyRotation (Rotation2D)
            Inputs: ('x', 'y')
            Outputs: ('x', 'y')
            Fittable parameters: ('angle',)
            >>> issubclass(SkyRotation, Rotation2D)
            True
            >>> r = SkyRotation(90)
            >>> isinstance(r, Rotation2D)
            True
        """

        if six.PY2 and isinstance(name, six.text_type):
            # Unicode names are not allowed in Python 2, so just convert to
            # ASCII.  As such, for cross-compatibility all model names should
            # just be ASCII for now.
            name = name.encode('ascii')

        mod = find_current_module(2)
        if mod:
            modname = mod.__name__
        else:
            modname = '__main__'

        new_cls = type(name, (cls,), {})
        # On Python 2 __module__ must be a str, not unicode
        new_cls.__module__ = str(modname)

        if hasattr(cls, '__qualname__'):
            if new_cls.__module__ == '__main__':
                # __main__ is not added to a class's qualified name
                new_cls.__qualname__ = name
            else:
                new_cls.__qualname__ = '{0}.{1}'.format(modname, name)

        return new_cls

    def _create_inverse_property(cls, members):
        inverse = members.get('inverse')
        if inverse is None or cls.__bases__[0] is object:
            # The latter clause is the prevent the below code from running on
            # the Model base class, which implements the default getter and
            # setter for .inverse
            return

        if isinstance(inverse, property):
            # We allow the @property decorator to be omitted entirely from
            # the class definition, though its use should be encouraged for
            # clarity
            inverse = inverse.fget

        # Store the inverse getter internally, then delete the given .inverse
        # attribute so that cls.inverse resolves to Model.inverse instead
        cls._inverse = inverse
        del cls.inverse

    def _create_bounding_box_property(cls, members):
        """
        Takes any bounding_box defined on a concrete Model subclass (either
        as a fixed tuple or a property or method) and wraps it in the generic
        getter/setter interface for the bounding_box attribute.
        """

        # TODO: Much of this is verbatim from _create_inverse_property--I feel
        # like there could be a way to generify properties that work this way,
        # but for the time being that would probably only confuse things more.
        bounding_box = members.get('bounding_box')
        if bounding_box is None or cls.__bases__[0] is object:
            return

        if isinstance(bounding_box, property):
            bounding_box = bounding_box.fget

        if not callable(bounding_box):
            # See if it's a hard-coded bounding_box (as a sequence) and
            # normalize it
            try:
                bounding_box = _BoundingBox.validate(cls, bounding_box)
            except AssertionError as exc:
                raise ModelDefinitionError(exc.args[0])
        else:
            sig = signature(bounding_box)
            # May be a method that only takes 'self' as an argument (like a
            # property, but the @property decorator was forgotten)
            # TODO: Maybe warn in the above case?
            #
            # However, if the method takes additional arguments then this is a
            # parameterized bounding box and should be callable
            if len(sig.parameters) > 1:
                bounding_box = \
                        cls._create_bounding_box_subclass(bounding_box, sig)

        if six.PY2 and isinstance(bounding_box, types.MethodType):
            bounding_box = bounding_box.__func__

        # See the Model.bounding_box getter definition for how this attribute
        # is used
        cls._bounding_box = bounding_box
        del cls.bounding_box

    def _create_bounding_box_subclass(cls, func, sig):
        """
        For Models that take optional arguments for defining their bounding
        box, we create a subclass of _BoundingBox with a ``__call__`` method
        that supports those additional arguments.

        Takes the function's Signature as an argument since that is already
        computed in _create_bounding_box_property, so no need to duplicate that
        effort.
        """

        # TODO: Might be convenient if calling the bounding box also
        # automatically sets the _user_bounding_box.  So that
        #
        #    >>> model.bounding_box(arg=1)
        #
        # in addition to returning the computed bbox, also sets it, so that
        # it's a shortcut for
        #
        #    >>> model.bounding_box = model.bounding_box(arg=1)
        #
        # Not sure if that would be non-obvious / confusing though...

        def __call__(self, **kwargs):
            return func(self._model, **kwargs)

        kwargs = []
        for idx, param in enumerate(sig.parameters.values()):
            if idx == 0:
                # Presumed to be a 'self' argument
                continue

            if param.default is param.empty:
                raise ModelDefinitionError(
                    'The bounding_box method for {0} is not correctly '
                    'defined: If defined as a method all arguments to that '
                    'method (besides self) must be keyword arguments with '
                    'default values that can be used to compute a default '
                    'bounding box.'.format(cls.name))

            kwargs.append((param.name, param.default))

        __call__ = make_function_with_signature(__call__, ('self',), kwargs)

        return type(str('_{0}BoundingBox'.format(cls.name)), (_BoundingBox,),
                    {'__call__': __call__})

    def _handle_backwards_compat(cls, name, members):
        # Backwards compatibility check for 'eval' -> 'evaluate'
        # TODO: Remove sometime after Astropy 1.0 release.
        if 'eval' in members and 'evaluate' not in members:
            warnings.warn(
                "Use of an 'eval' method when defining subclasses of "
                "FittableModel is deprecated; please rename this method to "
                "'evaluate'.  Otherwise its semantics remain the same.",
                AstropyDeprecationWarning)
            cls.evaluate = members['eval']
        elif ('evaluate' in members and callable(members['evaluate']) and
                not getattr(members['evaluate'], '__isabstractmethod__',
                            False)):
            # Don't bother making a deprecated eval() except for concrete
            # implementations of evaluate, so that we don't end up with an eval
            # abstractmethod as well
            alt = '.'.join((name, 'evaluate'))
            deprecate = deprecated('1.0', alternative=alt, name='eval')
            cls.eval = deprecate(members['evaluate'])

    def _handle_special_methods(cls, members):
        # Handle init creation from inputs
        def update_wrapper(wrapper, cls):
            # Set up the new __call__'s metadata attributes as though it were
            # manually defined in the class definition
            # A bit like functools.update_wrapper but uses the class instead of
            # the wrapped function
            wrapper.__module__ = cls.__module__
            wrapper.__doc__ = getattr(cls, wrapper.__name__).__doc__
            if hasattr(cls, '__qualname__'):
                wrapper.__qualname__ = '{0}.{1}'.format(
                        cls.__qualname__, wrapper.__name__)

        if ('__call__' not in members and 'inputs' in members and
                isinstance(members['inputs'], tuple)):
            inputs = members['inputs']
            # Done create a custom __call__ for classes that already have one
            # explicitly defined (this includes the Model base class, and any
            # other classes that manually override __call__
            def __call__(self, *inputs, **kwargs):
                """Evaluate this model on the supplied inputs."""

                return super(cls, self).__call__(*inputs, **kwargs)

            args = ('self',) + inputs
            new_call = make_function_with_signature(
                    __call__, args, [('model_set_axis', None)])
            update_wrapper(new_call, cls)
            cls.__call__ = new_call

        if ('__init__' not in members and not inspect.isabstract(cls) and
                cls._parameters_):
            # If *all* the parameters have default values we can make them
            # keyword arguments; otherwise they must all be positional
            # arguments
            if all(p.default is not None
                   for p in six.itervalues(cls._parameters_)):
                args = ('self',)
                kwargs = [(name, cls._parameters_[name].default)
                          for name in cls.param_names]
            else:
                args = ('self',) + cls.param_names
                kwargs = {}

            def __init__(self, *params, **kwargs):
                return super(cls, self).__init__(*params, **kwargs)

            new_init = make_function_with_signature(
                    __init__, args, kwargs, varkwargs='kwargs')
            update_wrapper(new_init, cls)
            cls.__init__ = new_init

    # *** Arithmetic operators for creating compound models ***
    __add__ =     _model_oper('+')
    __sub__ =     _model_oper('-')
    __mul__ =     _model_oper('*')
    __truediv__ = _model_oper('/')
    __pow__ =     _model_oper('**')
    __or__ =      _model_oper('|')
    __and__ =     _model_oper('&')

    if not six.PY3:
        # The classic __div__ operator need only be implemented for Python 2
        # without from __future__ import division
        __div__ = _model_oper('/')

    # *** Other utilities ***

    def _format_cls_repr(cls, keywords=[]):
        """
        Internal implementation of ``__repr__``.

        This is separated out for ease of use by subclasses that wish to
        override the default ``__repr__`` while keeping the same basic
        formatting.
        """

        # For the sake of familiarity start the output with the standard class
        # __repr__
        parts = [super(_ModelMeta, cls).__repr__()]

        if not cls._is_concrete:
            return parts[0]

        def format_inheritance(cls):
            bases = []
            for base in cls.mro()[1:]:
                if not issubclass(base, Model):
                    continue
                elif (inspect.isabstract(base) or
                        base.__name__.startswith('_')):
                    break
                bases.append(base.name)
            if bases:
                return '{0} ({1})'.format(cls.name, ' -> '.join(bases))
            else:
                return cls.name

        try:
            default_keywords = [
                ('Name', format_inheritance(cls)),
                ('Inputs', cls.inputs),
                ('Outputs', cls.outputs),
            ]

            if cls.param_names:
                default_keywords.append(('Fittable parameters',
                                         cls.param_names))

            for keyword, value in default_keywords + keywords:
                if value is not None:
                    parts.append('{0}: {1}'.format(keyword, value))

            return '\n'.join(parts)
        except:
            # If any of the above formatting fails fall back on the basic repr
            # (this is particularly useful in debugging)
            return parts[0]


@six.add_metaclass(_ModelMeta)
class Model(object):
    """
    Base class for all models.

    This is an abstract class and should not be instantiated directly.

    This class sets the constraints and other properties for all individual
    parameters and performs parameter validation.

    The following initialization arguments apply to the majority of Model
    subclasses by default (exceptions include specialized utility models
    like `~astropy.modeling.mappings.Mapping`).  Parametric models take all
    their parameters as arguments, followed by any of the following optional
    keyword arguments:

    Parameters
    ----------
    name : str, optional
        A human-friendly name associated with this model instance
        (particularly useful for identifying the individual components of a
        compound model).

    meta : dict, optional
        An optional dict of user-defined metadata to attach to this model.
        How this is used and interpreted is up to the user or individual use
        case.

    n_models : int, optional
        If given an integer greater than 1, a *model set* is instantiated
        instead of a single model.  This affects how the parameter arguments
        are interpreted.  In this case each parameter must be given as a list
        or array--elements of this array are taken along the first axis (or
        ``model_set_axis`` if specified), such that the Nth element is the
        value of that parameter for the Nth model in the set.

        See the section on model sets in the documentation for more details.

    model_set_axis : int, optional
        This argument only applies when creating a model set (i.e. ``n_models >
        1``).  It changes how parameter values are interpreted.  Normally the
        first axis of each input parameter array (properly the 0th axis) is
        taken as the axis corresponding to the model sets.  However, any axis
        of an input array may be taken as this "model set axis".  This accepts
        negative integers as well--for example use ``model_set_axis=-1`` if the
        last (most rapidly changing) axis should be associated with the model
        sets.

    fixed : dict, optional
        Dictionary ``{parameter_name: bool}`` setting the fixed constraint
        for one or more parameters.  `True` means the parameter is held fixed
        during fitting and is prevented from updates once an instance of the
        model has been created.

        Alternatively the `~astropy.modeling.Parameter.fixed` property of a
        parameter may be used to lock or unlock individual parameters.

    tied : dict, optional
        Dictionary ``{parameter_name: callable}`` of parameters which are
        linked to some other parameter. The dictionary values are callables
        providing the linking relationship.

        Alternatively the `~astropy.modeling.Parameter.tied` property of a
        parameter may be used to set the ``tied`` constraint on individual
        parameters.

    bounds : dict, optional
        Dictionary ``{parameter_name: value}`` of lower and upper bounds of
        parameters. Keys are parameter names. Values are a list of length 2
        giving the desired range for the parameter.

        Alternatively the `~astropy.modeling.Parameter.min` and
        `~astropy.modeling.Parameter.max` or
        ~astropy.modeling.Parameter.bounds` properties of a parameter may be
        used to set bounds on individual parameters.

    eqcons : list, optional
        List of functions of length n such that ``eqcons[j](x0, *args) == 0.0``
        in a successfully optimized problem.

    ineqcons : list, optional
        List of functions of length n such that ``ieqcons[j](x0, *args) >=
        0.0`` is a successfully optimized problem.

    Examples
    --------
    >>> from astropy.modeling import models
    >>> def tie_center(model):
    ...         mean = 50 * model.stddev
    ...         return mean
    >>> tied_parameters = {'mean': tie_center}

    Specify that ``'mean'`` is a tied parameter in one of two ways:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                        tied=tied_parameters)

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.mean.tied
    False
    >>> g1.mean.tied = tie_center
    >>> g1.mean.tied
    <function tie_center at 0x...>

    Fixed parameters:

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3,
    ...                        fixed={'stddev': True})
    >>> g1.stddev.fixed
    True

    or

    >>> g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=.3)
    >>> g1.stddev.fixed
    False
    >>> g1.stddev.fixed = True
    >>> g1.stddev.fixed
    True
    """

    parameter_constraints = Parameter.constraints
    """
    Primarily for informational purposes, these are the types of constraints
    that can be set on a model's parameters.
    """
    model_constraints = ('eqcons', 'ineqcons')
    """
    Primarily for informational purposes, these are the types of constraints
    that constrain model evaluation.
    """

    param_names = ()
    """
    Names of the parameters that describe models of this type.

    The parameters in this tuple are in the same order they should be passed in
    when initializing a model of a specific type.  Some types of models, such
    as polynomial models, have a different number of parameters depending on
    some other property of the model, such as the degree.

    When defining a custom model class the value of this attribute is
    automatically set by the `~astropy.modeling.Parameter` attributes defined
    in the class body.
    """

    inputs = ()
    """The name(s) of the input variable(s) on which a model is evaluated."""
    outputs = ()
    """The name(s) of the output(s) of the model."""

    standard_broadcasting = True
    fittable = False
    linear = True

    meta = metadata.MetaData()
    """A dict-like object to store optional information."""

    # By default models either use their own inverse property or have no
    # inverse at all, but users may also assign a custom inverse to a model,
    # optionally; in that case it is of course up to the user to determine
    # whether their inverse is *actually* an inverse to the model they assign
    # it to.
    _inverse = None
    _user_inverse = None

    _bounding_box = None
    _user_bounding_box = None

    # Default n_models attribute, so that __len__ is still defined even when a
    # model hasn't completed initialization yet
    _n_models = 1

    def __init__(self, *args, **kwargs):
        super(Model, self).__init__()
        meta = kwargs.pop('meta', None)
        if meta is not None:
            self.meta = meta

        self._name = kwargs.pop('name', None)

        self._initialize_constraints(kwargs)
        # Remaining keyword args are either parameter values or invalid
        # Parameter values must be passed in as keyword arguments in order to
        # distinguish them
        self._initialize_parameters(args, kwargs)

    def __repr__(self):
        return self._format_repr()

    def __str__(self):
        return self._format_str()

    def __len__(self):
        return self._n_models

    def __call__(self, *inputs, **kwargs):
        """
        Evaluate this model using the given input(s) and the parameter values
        that were specified when the model was instantiated.
        """

        inputs, format_info = self.prepare_inputs(*inputs, **kwargs)
        parameters = self._param_sets(raw=True)

        outputs = self.evaluate(*chain(inputs, parameters))

        if self.n_outputs == 1:
            outputs = (outputs,)

        return self.prepare_outputs(format_info, *outputs, **kwargs)

    # *** Arithmetic operators for creating compound models ***
    __add__ =     _model_oper('+')
    __sub__ =     _model_oper('-')
    __mul__ =     _model_oper('*')
    __truediv__ = _model_oper('/')
    __pow__ =     _model_oper('**')
    __or__ =      _model_oper('|')
    __and__ =     _model_oper('&')

    if not six.PY3:
        __div__ = _model_oper('/')

    # *** Properties ***
    @property
    def name(self):
        """User-provided name for this model instance."""

        return self._name

    @property
    @deprecated('0.4', alternative='len(model)')
    def param_dim(self):
        return self._n_models

    @property
    def n_inputs(self):
        """
        The number of inputs to this model.

        Equivalent to ``len(model.inputs)``.
        """

        return len(self.inputs)

    @property
    def n_outputs(self):
        """
        The number of outputs from this model.

        Equivalent to ``len(model.outputs)``.
        """
        return len(self.outputs)

    @property
    def model_set_axis(self):
        """
        The index of the model set axis--that is the axis of a parameter array
        that pertains to which model a parameter value pertains to--as
        specified when the model was initialized.

        See the documentation on `Model Sets
        <http://docs.astropy.org/en/stable/modeling/models.html#model-sets>`_
        for more details.
        """

        return self._model_set_axis

    @property
    def param_sets(self):
        """
        Return parameters as a pset.

        This is a list with one item per parameter set, which is an array of
        that parameter's values across all parameter sets, with the last axis
        associated with the parameter set.
        """

        return self._param_sets()

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
        if not self.param_names:
            # Trivial, but not unheard of
            return self._parameters

        start = self._param_metrics[self.param_names[0]]['slice'].start
        stop = self._param_metrics[self.param_names[-1]]['slice'].stop

        return self._parameters[start:stop]

    @parameters.setter
    def parameters(self, value):
        """
        Assigning to this attribute updates the parameters array rather than
        replacing it.
        """

        if not self.param_names:
            return

        start = self._param_metrics[self.param_names[0]]['slice'].start
        stop = self._param_metrics[self.param_names[-1]]['slice'].stop

        try:
            value = np.array(value).flatten()
            self._parameters[start:stop] = value
        except ValueError as e:
            raise InputParameterError(
                "Input parameter values not compatible with the model "
                "parameters array: {0}".format(e))

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

    @property
    def eqcons(self):
        """List of parameter equality constraints."""

        return self._constraints['eqcons']

    @property
    def ineqcons(self):
        """List of parameter inequality constraints."""

        return self._constraints['ineqcons']

    @property
    def inverse(self):
        """
        Returns a new `~astropy.modeling.Model` instance which performs the
        inverse transform, if an analytic inverse is defined for this model.

        Even on models that don't have an inverse defined, this property can be
        set with a manually-defined inverse, such a pre-computed or
        experimentally determined inverse (often given as a
        `~astropy.modeling.polynomial.PolynomialModel`, but not by
        requirement).

        A custom inverse can be deleted with ``del model.inverse``.  In this
        case the model's inverse is reset to its default, if a default exists
        (otherwise the default is to raise `NotImplementedError`).

        Note to authors of `~astropy.modeling.Model` subclasses:  To define an
        inverse for a model simply override this property to return the
        appropriate model representing the inverse.  The machinery that will
        make the inverse manually-overridable is added automatically by the
        base class.
        """

        if self._user_inverse is not None:
            return self._user_inverse
        elif self._inverse is not None:
            return self._inverse()

        raise NotImplementedError("An analytical inverse transform has not "
                                  "been implemented for this model.")

    @inverse.setter
    def inverse(self, value):
        if not isinstance(value, (Model, type(None))):
            raise ValueError(
                "The ``inverse`` attribute may be assigned a `Model` "
                "instance or `None` (where `None` restores the default "
                "inverse for this model if one is defined.")

        if value is None:
            warnings.warn(
                "Currently setting `model.inverse = None` resets the inverse "
                "to the default inverse (if one exists).  However, starting "
                "in Astropy 1.2, setting `model.inverse = None` explicitly "
                "forces a model to have no inverse (such that accessing "
                "`model.inverse` raises a NotImplementedError) even if that "
                "model's class has a default inverse.\n\n"
                "Instead, call `del model.inverse` to reset the inverse to "
                "its default (if a default exists for that model's class--"
                "otherwise the model is reset to having no inverse.",
                AstropyDeprecationWarning)

        self._user_inverse = value

    @inverse.deleter
    def inverse(self):
        """
        Resets the model's inverse to its default (if one exists, otherwise
        the model will have no inverse).
        """

        del self._user_inverse

    @property
    def has_user_inverse(self):
        """
        A flag indicating whether or not a custom inverse model has been
        assigned to this model by a user, via assignment to ``model.inverse``.
        """

        return self._user_inverse is not None

    @property
    def _custom_inverse(self):
        # Deprecated alias for _user_inverse--included solely for temporary
        # backwards compatibility for pyasdf--remove after Astropy v1.1
        # release.
        return self._user_inverse

    @property
    def bounding_box(self):
        """
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
            except AssertionError as exc:
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

    # *** Public methods ***

    @abc.abstractmethod
    def evaluate(self, *args, **kwargs):
        """Evaluate the model on some input variables."""

    def render(self, out=None, coords=None):
        """
        Evaluates a model on an input array. Evaluation is limited to
        a bounding box if the `Model.bounding_box` attribute is set.

        Parameters
        ----------
        out : `numpy.ndarray`, optional
            The array on which the model is to be evaluated.
        coords : array-like, optional
            Coordinate arrays mapping to ``arr``, such that
            ``arr[coords] == arr``.

        Returns
        -------
        out : `numpy.ndarray`
            The model evaluated on the input array if given, or else a new array from
            ``coords``.
            If ``out`` and ``coords`` are both `None`, the returned array is
            limited to the `Model.bounding_box` limits. If
            `Model.bounding_box` is `None`, ``arr`` or ``coords`` must be passed.

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
            # Check dimensions match out and model
            assert len(coords) == ndim
            if out is not None:
                assert coords[0].shape == out.shape
            else:
                out = np.zeros(coords[0].shape)

        if out is not None:
            try:
                assert out.ndim == ndim
            except AssertionError:
                raise AssertionError(
                    'The array and model must have the same number '
                    'of dimensions.')

        if bbox is not None:

            # assures position is at center pixel, important when using add_array
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

    def prepare_inputs(self, *inputs, **kwargs):
        """
        This method is used in `~astropy.modeling.Model.__call__` to ensure
        that all the inputs to the model can be broadcast into compatible
        shapes (if one or both of them are input as arrays), particularly if
        there are more than one parameter sets.
        """

        model_set_axis = kwargs.pop('model_set_axis', None)

        if model_set_axis is None:
            # By default the model_set_axis for the input is assumed to be the
            # same as that for the parameters the model was defined with
            # TODO: Ensure that negative model_set_axis arguments are respected
            model_set_axis = self.model_set_axis

        n_models = len(self)

        params = [getattr(self, name) for name in self.param_names]
        inputs = [np.asanyarray(_input, dtype=float) for _input in inputs]

        _validate_input_shapes(inputs, self.inputs, n_models,
                               model_set_axis, self.standard_broadcasting)

        # The input formatting required for single models versus a multiple
        # model set are different enough that they've been split into separate
        # subroutines
        if n_models == 1:
            return _prepare_inputs_single_model(self, params, inputs,
                                                **kwargs)
        else:
            return _prepare_inputs_model_set(self, params, inputs, n_models,
                                             model_set_axis, **kwargs)

    def prepare_outputs(self, format_info, *outputs, **kwargs):
        if len(self) == 1:
            return _prepare_outputs_single_model(self, outputs, format_info)
        else:
            return _prepare_outputs_model_set(self, outputs, format_info)

    @deprecated('1.0',
                alternative='Use Model operators (TODO: link to compound '
                            'model docs')
    def add_model(self, model, mode):
        """
        Create a CompositeModel by chaining the current model with the new one
        using the specified mode.

        Parameters
        ----------
        model : an instance of a subclass of Model
        mode :  string
               'parallel', 'serial', 'p' or 's'
               a flag indicating whether to combine the models
               in series or in parallel

        Returns
        -------
        model : CompositeModel
            an instance of CompositeModel
        """

        from ._compound_deprecated import (SummedCompositeModel,
                                           SerialCompositeModel)

        if mode in ['parallel', 'p']:
            return SummedCompositeModel([self, model])
        elif mode in ['serial', 's']:
            return SerialCompositeModel([self, model])
        else:
            raise InputParameterError("Unrecognized mode {0}".format(mode))

    def copy(self):
        """
        Return a copy of this model.

        Uses a deep copy so that all model attributes, including parameter
        values, are copied as well.
        """

        return copy.deepcopy(self)

    @sharedmethod
    def rename(self, name):
        """
        Return a copy of this model with a new name.
        """

        new_model = self.copy()
        new_model._name = name
        return new_model

    # *** Internal methods ***
    @sharedmethod
    def _from_existing(self, existing, param_names):
        """
        Creates a new instance of ``cls`` that shares its underlying parameter
        values with an existing model instance given by ``existing``.

        This is used primarily by compound models to return a view of an
        individual component of a compound model.  ``param_names`` should be
        the names of the parameters in the *existing* model to use as the
        parameters in this new model.  Its length should equal the number of
        parameters this model takes, so that it can map parameters on the
        existing model to parameters on this model one-to-one.
        """

        # Basically this is an alternative __init__
        if isinstance(self, type):
            # self is a class, not an instance
            needs_initialization = True
            dummy_args = (0,) * len(param_names)
            self = self.__new__(self, *dummy_args)
        else:
            needs_initialization = False
            self = self.copy()

        aliases = dict(zip(self.param_names, param_names))
        # This is basically an alternative _initialize_constraints
        constraints = {}
        for cons_type in self.parameter_constraints:
            orig = existing._constraints[cons_type]
            constraints[cons_type] = AliasDict(orig, aliases)

        self._constraints = constraints

        self._n_models = existing._n_models
        self._model_set_axis = existing._model_set_axis
        self._parameters = existing._parameters

        self._param_metrics = defaultdict(dict)
        for param_a, param_b in six.iteritems(aliases):
            # Take the param metrics info for the giving parameters in the
            # existing model, and hand them to the appropriate parameters in
            # the new model
            self._param_metrics[param_a] = existing._param_metrics[param_b]

        if needs_initialization:
            self.__init__(*dummy_args)

        return self

    def _initialize_constraints(self, kwargs):
        """
        Pop parameter constraint values off the keyword arguments passed to
        `Model.__init__` and store them in private instance attributes.
        """

        if hasattr(self, '_constraints'):
            # Skip constraint initialization if it has already been handled via
            # an alternate initialization
            return

        self._constraints = {}
        # Pop any constraints off the keyword arguments
        for constraint in self.parameter_constraints:
            values = kwargs.pop(constraint, {})
            self._constraints[constraint] = values.copy()

            # Update with default parameter constraints
            for param_name in self.param_names:
                param = getattr(self, param_name)

                # Parameters don't have all constraint types
                value = getattr(param, constraint)
                if value is not None:
                    self._constraints[constraint][param_name] = value

        for constraint in self.model_constraints:
            values = kwargs.pop(constraint, [])
            self._constraints[constraint] = values

    def _initialize_parameters(self, args, kwargs):
        """
        Initialize the _parameters array that stores raw parameter values for
        all parameter sets for use with vectorized fitting algorithms; on
        FittableModels the _param_name attributes actually just reference
        slices of this array.
        """

        if hasattr(self, '_parameters'):
            # Skip parameter initialization if it has already been handled via
            # an alternate initialization
            return

        n_models = None
        # Pop off param_dim and handle backwards compatibility
        if 'param_dim' in kwargs:
            n_models = kwargs.pop('param_dim')
            warnings.warn(
                'The param_dim argument to {0}.__init__ is deprecated; '
                'use n_models instead.  See also the model_set_axis argument '
                'and related discussion in the docstring for Model.'.format(
                    self.__class__.__name__), AstropyDeprecationWarning)
            if 'n_models' in kwargs:
                raise TypeError(
                    "param_dim and n_models cannot both be specified; use "
                    "n_models, as param_dim is deprecated")
        else:
            n_models = kwargs.pop('n_models', None)

        if not (n_models is None or
                    (isinstance(n_models, int) and n_models >=1)):
            raise ValueError(
                "n_models must be either None (in which case it is "
                "determined from the model_set_axis of the parameter initial "
                "values) or it must be a positive integer "
                "(got {0!r})".format(n_models))

        model_set_axis = kwargs.pop('model_set_axis', None)
        if model_set_axis is None:
            if n_models is not None and n_models > 1:
                # Default to zero
                model_set_axis = 0
            else:
                # Otherwise disable
                model_set_axis = False
        else:
            if not (model_set_axis is False or
                    (isinstance(model_set_axis, int) and
                        not isinstance(model_set_axis, bool))):
                raise ValueError(
                    "model_set_axis must be either False or an integer "
                    "specifying the parameter array axis to map to each "
                    "model in a set of models (got {0!r}).".format(
                        model_set_axis))

        # Process positional arguments by matching them up with the
        # corresponding parameters in self.param_names--if any also appear as
        # keyword arguments this presents a conflict
        params = {}
        if len(args) > len(self.param_names):
            raise TypeError(
                "{0}.__init__() takes at most {1} positional arguments ({2} "
                "given)".format(self.__class__.__name__, len(self.param_names),
                                len(args)))

        for idx, arg in enumerate(args):
            if arg is None:
                # A value of None implies using the default value, if exists
                continue
            params[self.param_names[idx]] = np.asanyarray(arg, dtype=np.float)

        # At this point the only remaining keyword arguments should be
        # parameter names; any others are in error.
        for param_name in self.param_names:
            if param_name in kwargs:
                if param_name in params:
                    raise TypeError(
                        "{0}.__init__() got multiple values for parameter "
                        "{1!r}".format(self.__class__.__name__, param_name))
                value = kwargs.pop(param_name)
                if value is None:
                    continue
                params[param_name] = np.asanyarray(value, dtype=np.float)

        if kwargs:
            # If any keyword arguments were left over at this point they are
            # invalid--the base class should only be passed the parameter
            # values, constraints, and param_dim
            for kwarg in kwargs:
                # Just raise an error on the first unrecognized argument
                raise TypeError(
                    '{0}.__init__() got an unrecognized parameter '
                    '{1!r}'.format(self.__class__.__name__, kwarg))

        self._model_set_axis = model_set_axis
        self._param_metrics = defaultdict(dict)

        # Determine the number of model sets: If the model_set_axis is
        # None then there is just one parameter set; otherwise it is determined
        # by the size of that axis on the first parameter--if the other
        # parameters don't have the right number of axes or the sizes of their
        # model_set_axis don't match an error is raised
        if model_set_axis is not False and n_models != 1 and params:
            max_ndim = 0
            if model_set_axis < 0:
                min_ndim = abs(model_set_axis)
            else:
                min_ndim = model_set_axis + 1

            for name, value in six.iteritems(params):
                param_ndim = np.ndim(value)
                if param_ndim < min_ndim:
                    raise InputParameterError(
                        "All parameter values must be arrays of dimension "
                        "at least {0} for model_set_axis={1} (the value "
                        "given for {2!r} is only {3}-dimensional)".format(
                            min_ndim, model_set_axis, name, param_ndim))

                max_ndim = max(max_ndim, param_ndim)

                if n_models is None:
                    # Use the dimensions of the first parameter to determine
                    # the number of model sets
                    n_models = value.shape[model_set_axis]
                elif value.shape[model_set_axis] != n_models:
                    raise InputParameterError(
                        "Inconsistent dimensions for parameter {0!r} for "
                        "{1} model sets.  The length of axis {2} must be the "
                        "same for all input parameter values".format(
                        name, n_models, model_set_axis))

            self._check_param_broadcast(params, max_ndim)
        else:
            if n_models is None:
                n_models = 1

            self._check_param_broadcast(params, None)

        self._n_models = n_models
        self._initialize_parameter_values(params)

    def _initialize_parameter_values(self, params):
        # self._param_metrics should have been initialized in
        # self._initialize_parameters
        param_metrics = self._param_metrics
        total_size = 0

        for name in self.param_names:
            if params.get(name) is None:
                default = getattr(self, name).default

                if default is None:
                    # No value was supplied for the parameter, and the
                    # parameter does not have a default--therefor the model is
                    # underspecified
                    raise TypeError(
                        "{0}.__init__() requires a value for parameter "
                        "{1!r}".format(self.__class__.__name__, name))

                value = params[name] = default
            else:
                value = params[name]

            param_size = np.size(value)
            param_shape = np.shape(value)

            param_slice = slice(total_size, total_size + param_size)

            param_metrics[name]['slice'] = param_slice
            param_metrics[name]['shape'] = param_shape
            total_size += param_size

        self._param_metrics = param_metrics
        self._parameters = np.empty(total_size, dtype=np.float64)

        # Now set the parameter values (this will also fill
        # self._parameters)
        # TODO: This is a bit ugly, but easier to deal with than how this was
        # done previously.  There's still lots of opportunity for refactoring
        # though, in particular once we move the _get/set_model_value methods
        # out of Parameter and into Model (renaming them
        # _get/set_parameter_value)
        for name, value in params.items():
            param_descr = getattr(self, name)
            value = np.array(value)
            if param_descr._setter is not None:
                value = param_descr._setter(value)
            self._parameters[param_metrics[name]['slice']] = value.ravel()

        # Finally validate all the parameters; we do this last so that
        # validators that depend on one of the other parameters' values will
        # work
        for name in params:
            param_descr = getattr(self, name)
            param_descr.validator(param_descr.value)

    def _check_param_broadcast(self, params, max_ndim):
        """
        This subroutine checks that all parameter arrays can be broadcast
        against each other, and determines the shapes parameters must have in
        order to broadcast correctly.

        If model_set_axis is None this merely checks that the parameters
        broadcast and returns an empty dict if so.  This mode is only used for
        single model sets.
        """

        all_shapes = []
        param_names = []
        model_set_axis = self._model_set_axis

        for name in self.param_names:
            # Previously this just used iteritems(params), but we loop over all
            # param_names instead just to ensure some determinism in the
            # ordering behavior
            if name not in params:
                continue

            value = params[name]
            param_names.append(name)
            # We've already checked that each parameter array is compatible in
            # the model_set_axis dimension, but now we need to check the
            # dimensions excluding that axis
            # Split the array dimensions into the axes before model_set_axis
            # and after model_set_axis
            param_shape = np.shape(value)

            param_ndim = len(param_shape)
            if max_ndim is not None and param_ndim < max_ndim:
                # All arrays have the same number of dimensions up to the
                # model_set_axis dimension, but after that they may have a
                # different number of trailing axes.  The number of trailing
                # axes must be extended for mutual compatibility.  For example
                # if max_ndim = 3 and model_set_axis = 0, an array with the
                # shape (2, 2) must be extended to (2, 1, 2).  However, an
                # array with shape (2,) is extended to (2, 1).
                new_axes = (1,) * (max_ndim - param_ndim)

                if model_set_axis < 0:
                    # Just need to prepend axes to make up the difference
                    broadcast_shape = new_axes + param_shape
                else:
                    broadcast_shape = (param_shape[:model_set_axis + 1] +
                                       new_axes +
                                       param_shape[model_set_axis + 1:])
                self._param_metrics[name]['broadcast_shape'] = broadcast_shape
                all_shapes.append(broadcast_shape)
            else:
                all_shapes.append(param_shape)

        # Now check mutual broadcastability of all shapes
        try:
            check_broadcast(*all_shapes)
        except IncompatibleShapeError as exc:
            shape_a, shape_a_idx, shape_b, shape_b_idx = exc.args
            param_a = param_names[shape_a_idx]
            param_b = param_names[shape_b_idx]

            raise InputParameterError(
                "Parameter {0!r} of shape {1!r} cannot be broadcast with "
                "parameter {2!r} of shape {3!r}.  All parameter arrays "
                "must have shapes that are mutually compatible according "
                "to the broadcasting rules.".format(param_a, shape_a,
                                                    param_b, shape_b))

    def _param_sets(self, raw=False):
        """
        Implementation of the Model.param_sets property.

        This internal implementation has a ``raw`` argument which controls
        whether or not to return the raw parameter values (i.e. the values that
        are actually stored in the ._parameters array, as opposed to the values
        displayed to users.  In most cases these are one in the same but there
        are currently a few exceptions.

        Note: This is notably an overcomplicated device and may be removed
        entirely in the near future.
        """

        param_metrics = self._param_metrics
        values = []
        for name in self.param_names:
            if raw:
                value = getattr(self, name)._raw_value
            else:
                value = getattr(self, name).value

            broadcast_shape = param_metrics[name].get('broadcast_shape')
            if broadcast_shape is not None:
                value = value.reshape(broadcast_shape)

            values.append(value)

        shapes = [np.shape(value) for value in values]

        if len(self) == 1:
            # Add a single param set axis to the parameter's value (thus
            # converting scalars to shape (1,) array values) for consistency
            values = [np.array([value]) for value in values]

        if len(set(shapes)) != 1:
            # If the parameters are not all the same shape, converting to an
            # array is going to produce an object array
            # However the way Numpy creates object arrays is tricky in that it
            # will recurse into array objects in the list and break them up
            # into separate objects.  Doing things this way ensures a 1-D
            # object array the elements of which are the individual parameter
            # arrays.  There's not much reason to do this over returning a list
            # except for consistency
            psets = np.empty(len(values), dtype=object)
            psets[:] = values
            return psets

        return np.array(values)

    def _format_repr(self, args=[], kwargs={}, defaults={}):
        """
        Internal implementation of ``__repr__``.

        This is separated out for ease of use by subclasses that wish to
        override the default ``__repr__`` while keeping the same basic
        formatting.
        """

        # TODO: I think this could be reworked to preset model sets better

        parts = [repr(a) for a in args]

        parts.extend(
            "{0}={1}".format(name,
                             array_repr_oneline(getattr(self, name).value))
            for name in self.param_names)

        if self.name is not None:
            parts.append('name={0!r}'.format(self.name))

        for kwarg, value in kwargs.items():
            if kwarg in defaults and defaults[kwarg] != value:
                continue
            parts.append('{0}={1!r}'.format(kwarg, value))

        if len(self) > 1:
            parts.append("n_models={0}".format(len(self)))

        return '<{0}({1})>'.format(self.__class__.__name__, ', '.join(parts))

    def _format_str(self, keywords=[]):
        """
        Internal implementation of ``__str__``.

        This is separated out for ease of use by subclasses that wish to
        override the default ``__str__`` while keeping the same basic
        formatting.
        """

        default_keywords = [
            ('Model', self.__class__.__name__),
            ('Name', self.name),
            ('Inputs', self.inputs),
            ('Outputs', self.outputs),
            ('Model set size', len(self))
        ]

        parts = ['{0}: {1}'.format(keyword, value)
                 for keyword, value in default_keywords + keywords
                 if value is not None]

        parts.append('Parameters:')

        if len(self) == 1:
            columns = [[getattr(self, name).value]
                       for name in self.param_names]
        else:
            columns = [getattr(self, name).value
                       for name in self.param_names]

        if columns:
            param_table = Table(columns, names=self.param_names)
            parts.append(indent(str(param_table), width=4))

        return '\n'.join(parts)

class FittableModel(Model):
    """
    Base class for models that can be fitted using the built-in fitting
    algorithms.
    """

    linear = False
    # derivative with respect to parameters
    fit_deriv = None
    """
    Function (similar to the model's `~Model.evaluate`) to compute the
    derivatives of the model with respect to its parameters, for use by fitting
    algorithms.  In other words, this computes the Jacobian matrix with respect
    to the model's parameters.
    """
    # Flag that indicates if the model derivatives with respect to parameters
    # are given in columns or rows
    col_fit_deriv = True
    fittable = True


class Fittable1DModel(FittableModel):
    """
    Base class for one-dimensional fittable models.

    This class provides an easier interface to defining new models.
    Examples can be found in `astropy.modeling.functional_models`.
    """

    inputs = ('x',)
    outputs = ('y',)


class Fittable2DModel(FittableModel):
    """
    Base class for two-dimensional fittable models.

    This class provides an easier interface to defining new models.
    Examples can be found in `astropy.modeling.functional_models`.
    """

    inputs = ('x', 'y')
    outputs = ('z',)


def _make_arithmetic_operator(oper):
    # We don't bother with tuple unpacking here for efficiency's sake, but for
    # documentation purposes:
    #
    #     f_eval, f_n_inputs, f_n_outputs = f
    #
    # and similarly for g
    def op(f, g):
        return (make_binary_operator_eval(oper, f[0], g[0]), f[1], f[2])

    return op


def _composition_operator(f, g):
    # We don't bother with tuple unpacking here for efficiency's sake, but for
    # documentation purposes:
    #
    #     f_eval, f_n_inputs, f_n_outputs = f
    #
    # and similarly for g
    return (lambda inputs, params: g[0](f[0](inputs, params), params),
            f[1], g[2])


def _join_operator(f, g):
    # We don't bother with tuple unpacking here for efficiency's sake, but for
    # documentation purposes:
    #
    #     f_eval, f_n_inputs, f_n_outputs = f
    #
    # and similarly for g
    return (lambda inputs, params: (f[0](inputs[:f[1]], params) +
                                    g[0](inputs[f[1]:], params)),
            f[1] + g[1], f[2] + g[2])


# TODO: Support a couple unary operators--at least negation?
BINARY_OPERATORS = {
    '+': _make_arithmetic_operator(operator.add),
    '-': _make_arithmetic_operator(operator.sub),
    '*': _make_arithmetic_operator(operator.mul),
    '/': _make_arithmetic_operator(operator.truediv),
    '**': _make_arithmetic_operator(operator.pow),
    '|': _composition_operator,
    '&': _join_operator
}


_ORDER_OF_OPERATORS = [('|',), ('&',), ('+', '-'), ('*', '/'), ('**',)]
OPERATOR_PRECEDENCE = {}
for idx, ops in enumerate(_ORDER_OF_OPERATORS):
    for op in ops:
        OPERATOR_PRECEDENCE[op] = idx
del idx, op, ops


class _CompoundModelMeta(_ModelMeta):
    _tree = None
    _submodels = None
    _submodel_names = None
    _nextid = 0

    _param_names = None
    # _param_map is a mapping of the compound model's generated param names to
    # the parameters of submodels they are associated with.  The values in this
    # mapping are (idx, name) tuples were idx is the index of the submodel this
    # parameter is associated with, and name is the same parameter's name on
    # the submodel
    # In principle this will allow compound models to give entirely new names
    # to parameters that don't have to be the same as their original names on
    # the submodels, but right now that isn't taken advantage of
    _param_map = None

    _slice_offset = 0
    # When taking slices of a compound model, this keeps track of how offset
    # the first model in the slice is from the first model in the original
    # compound model it was taken from

    # This just inverts _param_map, swapping keys with values.  This is also
    # useful to have.
    _param_map_inverse = None
    _fittable = None

    _evaluate = None

    def __getitem__(cls, index):
        index = cls._normalize_index(index)

        if isinstance(index, int):
            return cls._get_submodels()[index]
        else:
            return cls._get_slice(index.start, index.stop)

    def __getattr__(cls, attr):
        # Make sure the _tree attribute is set; otherwise we are not looking up
        # an attribute on a concrete compound model class and should just raise
        # the AttributeError
        if cls._tree is not None and attr in cls.param_names:
            cls._init_param_descriptors()
            return getattr(cls, attr)

        raise AttributeError(attr)

    def __repr__(cls):
        if cls._tree is None:
            # This case is mostly for debugging purposes
            return cls._format_cls_repr()

        expression = cls._format_expression()
        components = '\n\n'.join('[{0}]: {1!r}'.format(idx, m)
                                 for idx, m in enumerate(cls._get_submodels()))
        keywords = [
            ('Expression', expression),
            ('Components', '\n' + indent(components))
        ]

        return cls._format_cls_repr(keywords=keywords)

    def __dir__(cls, *args):
        """
        Returns a list of attributes defined on a compound model, including
        all of its parameters.
        """

        # The *args is to address a bug (?) on Python 2.6 where the dir()
        # builtin calls __dir__ with an additional (unused) argument

        try:
            # Annoyingly, this will only work for Python 3.3+
            basedir = super(_CompoundModelMeta, cls).__dir__()
        except AttributeError:
            basedir = list(set((dir(type(cls)) + list(cls.__dict__))))

        if cls._tree is not None:
            for name in cls.param_names:
                basedir.append(name)

            basedir.sort()

        return basedir

    def __reduce__(cls):
        rv = super(_CompoundModelMeta, cls).__reduce__()

        if isinstance(rv, tuple):
            # Delete _evaluate from the members dict
            with ignored(KeyError):
                del rv[1][2]['_evaluate']

        return rv


    @property
    def submodel_names(cls):
        if cls._submodel_names is not None:
            return cls._submodel_names

        by_name = defaultdict(list)

        for idx, submodel in enumerate(cls._get_submodels()):
            # Keep track of the original sort order of the submodels
            by_name[submodel.name].append(idx)

        names = []
        for basename, indices in six.iteritems(by_name):
            if len(indices) == 1:
                # There is only one model with this name, so it doesn't need an
                # index appended to its name
                names.append((basename, indices[0]))
            else:
                for idx in indices:
                    names.append(('{0}_{1}'.format(basename, idx), idx))

        # Sort according to the models' original sort orders
        names.sort(key=lambda k: k[1])

        names = tuple(k[0] for k in names)

        cls._submodels_names = names
        return names

    @property
    def param_names(cls):
        if cls._param_names is None:
            cls._init_param_names()

        return cls._param_names

    @property
    def fittable(cls):
        if cls._fittable is None:
            cls._fittable = all(m.fittable for m in cls._get_submodels())

        return cls._fittable

    # TODO: Maybe we could use make_function_with_signature for evaluate, but
    # it's probably not worth it (and I'm not sure what the limit is on number
    # of function arguments/local variables but we could break that limit for
    # complicated compound models...
    def evaluate(cls, *args):
        if cls._evaluate is None:
            func = cls._tree.evaluate(BINARY_OPERATORS,
                                      getter=cls._model_evaluate_getter)[0]
            # Making this a staticmethod isn't strictly necessary for Python 3,
            # but it is necessary on Python 2 since looking up cls._evaluate
            # will return an unbound method otherwise
            cls._evaluate = staticmethod(func)

        inputs = args[:cls.n_inputs]
        params = iter(args[cls.n_inputs:])
        result = cls._evaluate(inputs, params)

        if cls.n_outputs == 1:
            return result[0]
        else:
            return result

    # TODO: This supports creating a new compound model from two existing
    # compound models (or normal models) and a single operator.  However, it
    # ought also to be possible to create a new model from an *entire*
    # expression, represented as a sequence of operators and their operands (or
    # an exiting ExpressionTree) and build that into a compound model without
    # creating an intermediate _CompoundModel class for every single operator
    # in the expression.  This will prove to be a useful optimization in many
    # cases
    @classmethod
    def _from_operator(mcls, operator, left, right, additional_members={}):
        """
        Given a Python operator (represented by a string, such as ``'+'``
        or ``'*'``, and two model classes or instances, return a new compound
        model that evaluates the given operator on the outputs of the left and
        right input models.

        If either of the input models are a model *class* (i.e. a subclass of
        `~astropy.modeling.Model`) then the returned model is a new subclass of
        `~astropy.modeling.Model` that may be instantiated with any parameter
        values.  If both input models are *instances* of a model, a new class
        is still created, but this method returns an *instance* of that class,
        taking the parameter values from the parameters of the input model
        instances.

        If given, the ``additional_members`` `dict` may provide additional
        class members that should be added to the generated
        `~astropy.modeling.Model` subclass.  Some members that are generated by
        this method should not be provided by ``additional_members``.  These
        include ``_tree``, ``inputs``, ``outputs``, ``linear``,
        ``standard_broadcasting``, and ``__module__`.  This is currently for
        internal use only.
        """

        # Note, currently this only supports binary operators, but could be
        # easily extended to support unary operators (namely '-') if/when
        # needed
        children = []
        for child in (left, right):
            if isinstance(child, (_CompoundModelMeta, _CompoundModel)):
                children.append(child._tree)
            else:
                children.append(ExpressionTree(child))

        tree = ExpressionTree(operator, left=children[0], right=children[1])

        name = str('CompoundModel{0}'.format(_CompoundModelMeta._nextid))
        _CompoundModelMeta._nextid += 1

        mod = find_current_module(3)
        if mod:
            modname = mod.__name__
        else:
            modname = '__main__'

        # TODO: These aren't the full rules for handling inputs and outputs, but
        # this will handle most basic cases correctly
        if operator == '|':
            inputs = left.inputs
            outputs = right.outputs

            if left.n_outputs != right.n_inputs:
                raise ModelDefinitionError(
                    "Unsupported operands for |: {0} (n_inputs={1}, "
                    "n_outputs={2}) and {3} (n_inputs={4}, n_outputs={5}); "
                    "n_outputs for the left-hand model must match n_inputs "
                    "for the right-hand model.".format(
                        left.name, left.n_inputs, left.n_outputs, right.name,
                        right.n_inputs, right.n_outputs))
        elif operator == '&':
            inputs = combine_labels(left.inputs, right.inputs)
            outputs = combine_labels(left.outputs, right.outputs)
        else:
            # Without loss of generality
            inputs = left.inputs
            outputs = left.outputs

            if (left.n_inputs != right.n_inputs or
                    left.n_outputs != right.n_outputs):
                raise ModelDefinitionError(
                    "Unsupported operands for {0}: {1} (n_inputs={2}, "
                    "n_outputs={3}) and {4} (n_inputs={5}, n_outputs={6}); "
                    "models must have the same n_inputs and the same "
                    "n_outputs for this operator".format(
                        operator, left.name, left.n_inputs, left.n_outputs,
                        right.name, right.n_inputs, right.n_outputs))

        if operator in ('|', '+', '-'):
            linear = left.linear and right.linear
        else:
            # Which is not to say it is *definitely* not linear but it would be
            # trickier to determine
            linear = False

        standard_broadcasting = \
                left.standard_broadcasting and right.standard_broadcasting

        # Note: If any other members are added here, make sure to mention them
        # in the docstring of this method.
        members = additional_members
        members.update({
            '_tree': tree,
            '_is_dynamic': True,  # See docs for _ModelMeta._is_dynamic
            'inputs': inputs,
            'outputs': outputs,
            'linear': linear,
            'standard_broadcasting': standard_broadcasting,
            '__module__': str(modname)})

        new_cls = mcls(name, (_CompoundModel,), members)

        if isinstance(left, Model) and isinstance(right, Model):
            # Both models used in the operator were already instantiated models,
            # not model *classes*.  As such it's not particularly useful to return
            # the class itself, but to instead produce a new instance:
            instance = new_cls()

            # Workaround for https://github.com/astropy/astropy/issues/3542
            # TODO: Any effort to restructure the tree-like data structure for
            # compound models should try to obviate this workaround--if
            # intermediate compound models are stored in the tree as well then
            # we can immediately check for custom inverses on sub-models when
            # computing the inverse
            instance._user_inverse = mcls._make_user_inverse(
                    operator, left, right)

            return instance

        # Otherwise return the new uninstantiated class itself
        return new_cls

    @classmethod
    def _handle_backwards_compat(mcls, name, members):
        # Override _handle_backwards_compat from _ModelMeta to be a no-op; it
        # is not needed since compound models did not exist before version 1.0
        # anyways.

        # TODO: Remove this at the same time as removing
        # _ModelMeta._handle_backwards_compat
        return

    @classmethod
    def _make_user_inverse(mcls, operator, left, right):
        """
        Generates an inverse `Model` for this `_CompoundModel` when either
        model in the operation has a *custom inverse* that was manually
        assigned by the user.

        If either model has a custom inverse, and in particular if another
        `_CompoundModel` has a custom inverse, then none of that model's
        sub-models should be considered at all when computing the inverse.
        So in that case we just compute the inverse ahead of time and set
        it as the new compound model's custom inverse.

        Note, this use case only applies when combining model instances,
        since model classes don't currently have a notion of a "custom
        inverse" (though it could probably be supported by overriding the
        class's inverse property).

        TODO: Consider fixing things so the aforementioned class-based case
        works as well.  However, for the present purposes this is good enough.
        """

        if not (operator in ('&', '|') and
                (left._user_inverse or right._user_inverse)):
            # These are the only operators that support an inverse right now
            return None

        try:
            left_inv = left.inverse
            right_inv = right.inverse
        except NotImplementedError:
            # If either inverse is undefined then just return False; this
            # means the normal _CompoundModel.inverse routine will fail
            # naturally anyways, since it requires all sub-models to have
            # an inverse defined
            return None

        if operator == '&':
            return left_inv & right_inv
        else:
            return right_inv | left_inv

    # TODO: Perhaps, just perhaps, the post-order (or ???-order) ordering of
    # leaf nodes is something the ExpressionTree class itself could just know
    def _get_submodels(cls):
        # Would make this a lazyproperty but those don't currently work with
        # type objects
        if cls._submodels is not None:
            return cls._submodels

        submodels = [c.value for c in cls._tree.traverse_postorder()
                     if c.isleaf]
        cls._submodels = submodels
        return submodels

    def _init_param_descriptors(cls):
        """
        This routine sets up the names for all the parameters on a compound
        model, including figuring out unique names for those parameters and
        also mapping them back to their associated parameters of the underlying
        submodels.

        Setting this all up is costly, and only necessary for compound models
        that a user will directly interact with.  For example when building an
        expression like::

            >>> M = (Model1 + Model2) * Model3  # doctest: +SKIP

        the user will generally never interact directly with the temporary
        result of the subexpression ``(Model1 + Model2)``.  So there's no need
        to setup all the parameters for that temporary throwaway.  Only once
        the full expression is built and the user initializes or introspects
        ``M`` is it necessary to determine its full parameterization.
        """

        # Accessing cls.param_names will implicitly call _init_param_names if
        # needed and thus also set up the _param_map; I'm not crazy about that
        # design but it stands for now
        for param_name in cls.param_names:
            submodel_idx, submodel_param = cls._param_map[param_name]
            submodel = cls[submodel_idx]

            orig_param = getattr(submodel, submodel_param, None)
            if not isinstance(orig_param, Parameter):
                # This is just a pathological case that is only really needed
                # to support the deprecated _CompositeModel--composite models
                # claim to have some parameters, but don't actually implement
                # the parameter descriptors, so we just make one up basically,
                # with a default value of zero.  This value will just be thrown
                # away, basically.
                # TODO: Remove this special case once the legacy interfaces
                # have been removed (basically this entire if statement--keep
                # only the parts in the else: clause.
                new_param = Parameter(name=param_name, default=0)
            else:
                if isinstance(submodel, Model):
                    # Take the parameter's default from the model's value for that
                    # parameter
                    default = orig_param.value
                else:
                    default = orig_param.default

                # Copy constraints
                constraints = dict((key, getattr(orig_param, key))
                                   for key in Model.parameter_constraints)

                # Note: Parameter.copy() returns a new unbound Parameter, never
                # a bound Parameter even if submodel is a Model instance (as
                # opposed to a Model subclass)
                new_param = orig_param.copy(name=param_name, default=default,
                                            **constraints)

            setattr(cls, param_name, new_param)

    def _init_param_names(cls):
        """
        This subroutine is solely for setting up the ``param_names`` attribute
        itself.

        See ``_init_param_descriptors`` for the full parameter setup.
        """

        # Currently this skips over Model *instances* in the expression tree;
        # basically these are treated as constants and do not add
        # fittable/tunable parameters to the compound model.
        # TODO: I'm not 100% happy with this design, and maybe we need some
        # interface for distinguishing fittable/settable parameters with
        # *constant* parameters (which would be distinct from parameters with
        # fixed constraints since they're permanently locked in place). But I'm
        # not sure if this is really the best way to treat the issue.

        names = []
        param_map = {}

        # Start counting the suffix indices to put on parameter names from the
        # slice_offset.  Usually this will just be zero, but for compound
        # models that were sliced from another compound model this may be > 0
        param_suffix = cls._slice_offset

        for idx, model in enumerate(cls._get_submodels()):
            if not model.param_names:
                # Skip models that don't have parameters in the numbering
                # TODO: Reevaluate this if it turns out to be confusing, though
                # parameter-less models are not very common in practice (there
                # are a few projections that don't take parameters)
                continue

            for param_name in model.param_names:
                # This is sort of heuristic, but we want to check that
                # model.param_name *actually* returns a Parameter descriptor,
                # and that the model isn't some inconsistent type that happens
                # to have a param_names attribute but does not actually
                # implement settable parameters.
                # In the future we can probably remove this check, but this is
                # here specifically to support the legacy compat
                # _CompositeModel which can be considered a pathological case
                # in the context of the new framework
                #if not isinstance(getattr(model, param_name, None),
                #                  Parameter):
                #    break
                name = '{0}_{1}'.format(param_name, param_suffix + idx)
                names.append(name)
                param_map[name] = (idx, param_name)

        cls._param_names = tuple(names)
        cls._param_map = param_map
        cls._param_map_inverse = dict((v, k) for k, v in param_map.items())

    def _format_expression(cls):
        # TODO: At some point might be useful to make a public version of this,
        # albeit with more formatting options
        return cls._tree.format_expression(OPERATOR_PRECEDENCE)

    def _normalize_index(cls, index):
        """
        Converts an index given to __getitem__ to either an integer, or
        a slice with integer start and stop values.

        If the length of the slice is exactly 1 this converts the index to a
        simple integer lookup.

        Negative integers are converted to positive integers.
        """

        def get_index_from_name(name):
            try:
                return cls.submodel_names.index(name)
            except ValueError:
                raise IndexError(
                    'Compound model {0} does not have a component named '
                    '{1}'.format(cls.name, name))

        def check_for_negative_index(index):
            if index < 0:
                new_index = len(cls.submodel_names) + index
                if new_index < 0:
                    # If still < 0 then this is an invalid index
                    raise IndexError(
                            "Model index {0} out of range.".format(index))
                else:
                    index = new_index

            return index

        if isinstance(index, six.string_types):
            return get_index_from_name(index)
        elif isinstance(index, slice):
            if index.step not in (1, None):
                # In principle it could be but I can scarcely imagine a case
                # where it would be useful.  If someone can think of one then
                # we can enable it.
                raise ValueError(
                    "Step not supported for compound model slicing.")
            start = index.start if index.start is not None else 0
            stop = (index.stop
                    if index.stop is not None else len(cls.submodel_names))
            if isinstance(start, int):
                start = check_for_negative_index(start)
            if isinstance(stop, int):
                stop = check_for_negative_index(stop)
            if isinstance(start, six.string_types):
                start = get_index_from_name(start)
            if isinstance(stop, six.string_types):
                stop = get_index_from_name(stop) + 1
            length = stop - start

            if length == 1:
                return start
            elif length <= 0:
                raise ValueError("Empty slice of a compound model.")

            return slice(start, stop)
        elif isinstance(index, int):
            if index >= len(cls.submodel_names):
                raise IndexError(
                        "Model index {0} out of range.".format(index))

            return check_for_negative_index(index)

        raise TypeError(
            'Submodels can be indexed either by their integer order or '
            'their name (got {0!r}).'.format(index))

    def _get_slice(cls, start, stop):
        """
        Return a new model build from a sub-expression of the expression
        represented by this model.

        Right now this is highly inefficient, as it creates a new temporary
        model for each operator that appears in the sub-expression.  It would
        be better if this just built a new expression tree, and the new model
        instantiated directly from that tree.

        Once tree -> model instantiation is possible this should be fixed to
        use that instead.
        """

        members = {'_slice_offset': cls._slice_offset + start}
        operators = dict((oper, _model_oper(oper, additional_members=members))
                         for oper in BINARY_OPERATORS)

        return cls._tree.evaluate(operators, start=start, stop=stop)

    @staticmethod
    def _model_evaluate_getter(idx, model):
        n_params = len(model.param_names)
        n_inputs = model.n_inputs
        n_outputs = model.n_outputs

        # There is currently an unfortunate inconsistency in some models, which
        # requires them to be instantiated for their evaluate to work.  I think
        # that needs to be reconsidered and fixed somehow, but in the meantime
        # we need to check for that case
        if (not isinstance(model, Model) and
                isinstancemethod(model, model.evaluate)):
            if n_outputs == 1:
                # Where previously model was a class, now make an instance
                def f(inputs, params):
                    param_values = tuple(islice(params, n_params))
                    return (model(*param_values).evaluate(
                        *chain(inputs, param_values)),)
            else:
                def f(inputs, params):
                    param_values = tuple(islice(params, n_params))
                    return model(*param_values).evaluate(
                        *chain(inputs, param_values))
        else:
            evaluate = model.evaluate
            if n_outputs == 1:
                f = lambda inputs, params: \
                    (evaluate(*chain(inputs, islice(params, n_params))),)
            else:
                f = lambda inputs, params: \
                    evaluate(*chain(inputs, islice(params, n_params)))

        return (f, n_inputs, n_outputs)


@six.add_metaclass(_CompoundModelMeta)
class _CompoundModel(Model):
    fit_deriv = None
    col_fit_deriv = False

    _submodels = None

    def __getattr__(self, attr):
        value = getattr(self.__class__, attr)
        if hasattr(value, '__get__'):
            # Object is a descriptor, so we should really return the result of
            # its __get__
            value = value.__get__(self, self.__class__)
        return value

    def __getitem__(self, index):
        index = self.__class__._normalize_index(index)
        model = self.__class__[index]

        if isinstance(index, slice):
            param_names = model.param_names
        else:
            param_map = self.__class__._param_map_inverse
            param_names = tuple(param_map[index, name]
                                for name in model.param_names)

        return model._from_existing(self, param_names)

    if sys.version_info[:3] < (2, 7, 3):
        def __reduce__(self):
            # _CompoundModel classes have a generated evaluate() that is cached
            # off in the _evaluate attribute.  This can't be pickled, and so
            # should be regenerated after unpickling (alas)
            if find_current_module(2) is not copy:
                # The copy module also uses __reduce__, but there's no problem
                # there.
                raise RuntimeError(
                    "Pickling of compound models is not possible using Python "
                    "versions less than 2.7.3 due to a bug in Python.  See "
                    "http://docs.astropy.org/en/v1.0.4/known_issues.html#"
                    "pickling-error-on-compound-models for more information ("
                    "tried to pickle {0!r}).".format(self))
            else:
                return super(_CompoundModel, self).__reduce__()

    @property
    def submodel_names(self):
        return self.__class__.submodel_names

    @property
    def param_names(self):
        return self.__class__.param_names

    @property
    def fittable(self):
        return self.__class__.fittable

    @sharedmethod
    def evaluate(self, *args):
        return self.__class__.evaluate(*args)

    # TODO: The way this works is highly inefficient--the inverse is created by
    # making a new model for each operator in the compound model, which could
    # potentially mean creating a large number of temporary throwaway model
    # classes.  This can definitely be optimized in the future by implementing
    # a way to construct a single model class from an existing tree
    @property
    def inverse(self):
        def _not_implemented(oper):
            def _raise(x, y):
                raise NotImplementedError(
                    "The inverse is not currently defined for compound "
                    "models created using the {0} operator.".format(oper))
            return _raise

        operators = dict((oper, _not_implemented(oper))
                         for oper in ('+', '-', '*', '/', '**'))
        operators['&'] = operator.and_
        # Reverse the order of compositions
        operators['|'] = lambda x, y: operator.or_(y, x)

        leaf_idx = -1

        def getter(idx, model):
            try:
                # By indexing on self[] this will return an instance of the
                # model, with all the appropriate parameters set, which is
                # currently required to return an inverse
                return self[idx].inverse
            except NotImplementedError:
                raise NotImplementedError(
                    "All models in a composite model must have an inverse "
                    "defined in order for the composite model to have an "
                    "inverse.  {0!r} does not have an inverse.".format(model))


        return self._tree.evaluate(operators, getter=getter)

    @sharedmethod
    def _get_submodels(self):
        return self.__class__._get_submodels()


def custom_model(*args, **kwargs):
    """
    Create a model from a user defined function. The inputs and parameters of
    the model will be inferred from the arguments of the function.

    This can be used either as a function or as a decorator.  See below for
    examples of both usages.

    .. note::

        All model parameters have to be defined as keyword arguments with
        default values in the model function.  Use `None` as a default argument
        value if you do not want to have a default value for that parameter.

    Parameters
    ----------
    func : function
        Function which defines the model.  It should take N positional
        arguments where ``N`` is dimensions of the model (the number of
        independent variable in the model), and any number of keyword arguments
        (the parameters).  It must return the value of the model (typically as
        an array, but can also be a scalar for scalar inputs).  This
        corresponds to the `~astropy.modeling.Model.evaluate` method.
    fit_deriv : function, optional
        Function which defines the Jacobian derivative of the model. I.e., the
        derivative with respect to the *parameters* of the model.  It should
        have the same argument signature as ``func``, but should return a
        sequence where each element of the sequence is the derivative
        with respect to the corresponding argument. This corresponds to the
        :meth:`~astropy.modeling.FittableModel.fit_deriv` method.

    Examples
    --------
    Define a sinusoidal model function as a custom 1D model::

        >>> from astropy.modeling.models import custom_model
        >>> import numpy as np
        >>> def sine_model(x, amplitude=1., frequency=1.):
        ...     return amplitude * np.sin(2 * np.pi * frequency * x)
        >>> def sine_deriv(x, amplitude=1., frequency=1.):
        ...     return 2 * np.pi * amplitude * np.cos(2 * np.pi * frequency * x)
        >>> SineModel = custom_model(sine_model, fit_deriv=sine_deriv)

    Create an instance of the custom model and evaluate it::

        >>> model = SineModel()
        >>> model(0.25)
        1.0

    This model instance can now be used like a usual astropy model.

    The next example demonstrates a 2D Moffat function model, and also
    demonstrates the support for docstrings (this example could also include
    a derivative, but it has been omitted for simplicity)::

        >>> @custom_model
        ... def Moffat2D(x, y, amplitude=1.0, x_0=0.0, y_0=0.0, gamma=1.0,
        ...            alpha=1.0):
        ...     \"\"\"Two dimensional Moffat function.\"\"\"
        ...     rr_gg = ((x - x_0) ** 2 + (y - y_0) ** 2) / gamma ** 2
        ...     return amplitude * (1 + rr_gg) ** (-alpha)
        ...
        >>> print(Moffat2D.__doc__)
        Two dimensional Moffat function.
        >>> model = Moffat2D()
        >>> model(1, 1)  # doctest: +FLOAT_CMP
        0.3333333333333333
    """

    fit_deriv = kwargs.get('fit_deriv', None)

    if len(args) == 1 and six.callable(args[0]):
        return _custom_model_wrapper(args[0], fit_deriv=fit_deriv)
    elif not args:
        return functools.partial(_custom_model_wrapper, fit_deriv=fit_deriv)
    else:
        raise TypeError(
            "{0} takes at most one positional argument (the callable/"
            "function to be turned into a model.  When used as a decorator "
            "it should be passed keyword arguments only (if "
            "any).".format(__name__))


def _custom_model_wrapper(func, fit_deriv=None):
    """
    Internal implementation `custom_model`.

    When `custom_model` is called as a function its arguments are passed to
    this function, and the result of this function is returned.

    When `custom_model` is used as a decorator a partial evaluation of this
    function is returned by `custom_model`.
    """

    if not six.callable(func):
        raise ModelDefinitionError(
            "func is not callable; it must be a function or other callable "
            "object")

    if fit_deriv is not None and not six.callable(fit_deriv):
        raise ModelDefinitionError(
            "fit_deriv not callable; it must be a function or other "
            "callable object")

    model_name = func.__name__

    inputs, params = get_inputs_and_params(func)

    if (fit_deriv is not None and
            len(six.get_function_defaults(fit_deriv)) != len(params)):
        raise ModelDefinitionError("derivative function should accept "
                                   "same number of parameters as func.")

    # TODO: Maybe have a clever scheme for default output name?
    if inputs:
        output_names = (inputs[0].name,)
    else:
        output_names = ('x',)

    params = dict((param.name, Parameter(param.name, default=param.default))
                  for param in params)

    mod = find_current_module(2)
    if mod:
        modname = mod.__name__
    else:
        modname = '__main__'

    members = {
        '__module__': str(modname),
        '__doc__': func.__doc__,
        'inputs': tuple(x.name for x in inputs),
        'outputs': output_names,
        'evaluate': staticmethod(func),
    }

    if fit_deriv is not None:
        members['fit_deriv'] = staticmethod(fit_deriv)

    members.update(params)

    return type(model_name, (FittableModel,), members)


def render_model(model, arr=None, coords=None):
    """
    Evaluates a model on an input array. Evaluation is limited to
    a bounding box if the `Model.bounding_box` attribute is set.

    Parameters
    ----------
    model : `Model`
        Model to be evaluated.
    arr : `numpy.ndarray`, optional
        Array on which the model is evaluated.
    coords : array-like, optional
        Coordinate arrays mapping to ``arr``, such that
        ``arr[coords] == arr``.

    Returns
    -------
    array : `numpy.ndarray`
        The model evaluated on the input ``arr`` or a new array from ``coords``.
        If ``arr`` and ``coords`` are both `None`, the returned array is
        limited to the `Model.bounding_box` limits. If
        `Model.bounding_box` is `None`, ``arr`` or ``coords`` must be passed.

    Examples
    --------
    :ref:`bounding-boxes`
    """

    bbox = model.bounding_box

    if (coords is None) & (arr is None) & (bbox is None):
        raise AssertionError('If no bounding_box is set, coords or arr must be input.')

    # for consistent indexing
    if model.n_inputs == 1:
        if coords is not None:
            coords = [coords]
        if bbox is not None:
            bbox = [bbox]

    if arr is not None:
        arr = arr.copy()
        # Check dimensions match model
        assert arr.ndim == model.n_inputs

    if coords is not None:
        # Check dimensions match arr and model
        coords = np.array(coords)
        assert len(coords) == model.n_inputs
        if arr is not None:
            assert coords[0].shape == arr.shape
        else:
            arr = np.zeros(coords[0].shape)

    if bbox is not None:
        # assures position is at center pixel, important when using add_array
        pd = pos, delta = np.array([(np.mean(bb), np.ceil((bb[1] - bb[0]) / 2))
                                    for bb in bbox]).astype(int).T

        if coords is not None:
            sub_shape = tuple(delta * 2 + 1)
            sub_coords = np.array([extract_array(c, sub_shape, pos) for c in coords])
        else:
            limits = [slice(p - d, p + d + 1, 1) for p, d in pd.T]
            sub_coords = np.mgrid[limits]

        sub_coords = sub_coords[::-1]

        if arr is None:
            arr = model(*sub_coords)
        else:
            try:
                arr = add_array(arr, model(*sub_coords), pos)
            except ValueError:
                raise ValueError('The `bounding_box` is larger than the input'
                                ' arr in one or more dimensions. Set '
                                '`model.bounding_box = None`.')
    else:

        if coords is None:
            im_shape = arr.shape
            limits = [slice(i) for i in im_shape]
            coords = np.mgrid[limits]

        arr += model(*coords[::-1])

    return arr


def _prepare_inputs_single_model(model, params, inputs, **kwargs):
    broadcasts = []

    for idx, _input in enumerate(inputs):
        input_shape = _input.shape

        # Ensure that array scalars are always upgrade to 1-D arrays for the
        # sake of consistency with how parameters work.  They will be cast back
        # to scalars at the end
        if not input_shape:
            inputs[idx] = _input.reshape((1,))

        if not params:
            max_broadcast = input_shape
        else:
            max_broadcast = ()

        for param in params:
            try:
                if model.standard_broadcasting:
                    broadcast = check_broadcast(input_shape, param.shape)
                else:
                    broadcast = input_shape
            except IncompatibleShapeError:
                raise ValueError(
                    "Model input argument {0!r} of shape {1!r} cannot be "
                    "broadcast with parameter {2!r} of shape "
                    "{3!r}.".format(model.inputs[idx], input_shape,
                                    param.name, param.shape))

            if len(broadcast) > len(max_broadcast):
                max_broadcast = broadcast
            elif len(broadcast) == len(max_broadcast):
                max_broadcast = max(max_broadcast, broadcast)

        broadcasts.append(max_broadcast)

    if model.n_outputs > model.n_inputs:
        if len(set(broadcasts)) > 1:
            raise ValueError(
                "For models with n_outputs > n_inputs, the combination of "
                "all inputs and parameters must broadcast to the same shape, "
                "which will be used as the shape of all outputs.  In this "
                "case some of the inputs had different shapes, so it is "
                "ambiguous how to format outputs for this model.  Try using "
                "inputs that are all the same size and shape.")
        else:
            # Extend the broadcasts list to include shapes for all outputs
            extra_outputs = model.n_outputs - model.n_inputs
            if not broadcasts:
                # If there were no inputs then the broadcasts list is empty
                # just add a None since there is no broadcasting of outputs and
                # inputs necessary (see _prepare_outputs_single_model)
                broadcasts.append(None)
            broadcasts.extend([broadcasts[0]] * extra_outputs)

    return inputs, (broadcasts,)


def _prepare_outputs_single_model(model, outputs, format_info):
    broadcasts = format_info[0]

    outputs = list(outputs)

    for idx, output in enumerate(outputs):
        broadcast_shape = broadcasts[idx]
        if broadcast_shape is not None:
            if not broadcast_shape:
                # Shape is (), i.e. a scalar should be returned
                outputs[idx] = np.asscalar(output)
            else:
                outputs[idx] = output.reshape(broadcast_shape)

    if model.n_outputs == 1:
        return outputs[0]
    else:
        return tuple(outputs)


def _prepare_inputs_model_set(model, params, inputs, n_models, model_set_axis,
                              **kwargs):
    reshaped = []
    pivots = []

    for idx, _input in enumerate(inputs):
        max_param_shape = ()

        if n_models > 1 and model_set_axis is not False:
            # Use the shape of the input *excluding* the model axis
            input_shape = (_input.shape[:model_set_axis] +
                           _input.shape[model_set_axis + 1:])
        else:
            input_shape = _input.shape

        for param in params:
            try:
                check_broadcast(input_shape, param.shape)
            except IncompatibleShapeError:
                raise ValueError(
                    "Model input argument {0!r} of shape {1!r} cannot be "
                    "broadcast with parameter {2!r} of shape "
                    "{3!r}.".format(model.inputs[idx], input_shape,
                                    param.name, param.shape))

            if len(param.shape) > len(max_param_shape):
                max_param_shape = param.shape

        # We've now determined that, excluding the model_set_axis, the
        # input can broadcast with all the parameters
        input_ndim = len(input_shape)
        if model_set_axis is False:
            if len(max_param_shape) > input_ndim:
                # Just needs to prepend new axes to the input
                n_new_axes = 1 + len(max_param_shape) - input_ndim
                new_axes = (1,) * n_new_axes
                new_shape = new_axes + _input.shape
                pivot = model.model_set_axis
            else:
                pivot = input_ndim - len(max_param_shape)
                new_shape = (_input.shape[:pivot] + (1,) +
                             _input.shape[pivot:])
            new_input = _input.reshape(new_shape)
        else:
            if len(max_param_shape) >= input_ndim:
                n_new_axes = len(max_param_shape) - input_ndim
                pivot = model.model_set_axis
                new_axes = (1,) * n_new_axes
                new_shape = (_input.shape[:pivot + 1] + new_axes +
                             _input.shape[pivot + 1:])
                new_input = _input.reshape(new_shape)
            else:
                pivot = _input.ndim - len(max_param_shape) - 1
                new_input = np.rollaxis(_input, model_set_axis,
                                        pivot + 1)

        pivots.append(pivot)
        reshaped.append(new_input)

    if model.n_inputs < model.n_outputs:
        pivots.extend([model_set_axis] * (model.n_outputs - model.n_inputs))

    return reshaped, (pivots,)


def _prepare_outputs_model_set(model, outputs, format_info):
    pivots = format_info[0]

    outputs = list(outputs)

    for idx, output in enumerate(outputs):
        pivot = pivots[idx]
        if pivot < output.ndim and pivot != model.model_set_axis:
            outputs[idx] = np.rollaxis(output, pivot,
                                       model.model_set_axis)

    if model.n_outputs == 1:
        return outputs[0]
    else:
        return tuple(outputs)


def _validate_input_shapes(inputs, argnames, n_models, model_set_axis,
                           validate_broadcasting):
    """
    Perform basic validation of model inputs--that they are mutually
    broadcastable and that they have the minimum dimensions for the given
    model_set_axis.

    If validation succeeds, returns the total shape that will result from
    broadcasting the input arrays with each other.
    """

    check_model_set_axis = n_models > 1 and model_set_axis is not False

    if not (validate_broadcasting or check_model_set_axis):
        # Nothing else needed here
        return

    all_shapes = []

    for idx, _input in enumerate(inputs):
        input_shape = np.shape(_input)
        # Ensure that the input's model_set_axis matches the model's
        # n_models
        if input_shape and check_model_set_axis:
            # Note: Scalar inputs *only* get a pass on this
            if len(input_shape) < model_set_axis + 1:
                raise ValueError(
                    "For model_set_axis={0}, all inputs must be at "
                    "least {1}-dimensional.".format(
                        model_set_axis, model_set_axis + 1))
            elif input_shape[model_set_axis] != n_models:
                raise ValueError(
                    "Input argument {0!r} does not have the correct "
                    "dimensions in model_set_axis={1} for a model set with "
                    "n_models={2}.".format(argnames[idx], model_set_axis,
                                           n_models))
        all_shapes.append(input_shape)

    if not validate_broadcasting:
        return

    try:
        input_broadcast = check_broadcast(*all_shapes)
    except IncompatibleShapeError as exc:
        shape_a, shape_a_idx, shape_b, shape_b_idx = exc.args
        arg_a = argnames[shape_a_idx]
        arg_b = argnames[shape_b_idx]

        raise ValueError(
            "Model input argument {0!r} of shape {1!r} cannot "
            "be broadcast with input {2!r} of shape {3!r}".format(
                arg_a, shape_a, arg_b, shape_b))

    return input_broadcast


copyreg.pickle(_ModelMeta, _ModelMeta.__reduce__)
copyreg.pickle(_CompoundModelMeta, _CompoundModelMeta.__reduce__)

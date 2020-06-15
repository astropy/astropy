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
# pylint: disable=invalid-name, protected-access, redefined-outer-name
import abc
import copy
import inspect
import functools
import operator
import types
import warnings

from collections import defaultdict, OrderedDict, deque
from inspect import signature
from itertools import chain

import numpy as np

from astropy.utils import indent, metadata
from astropy.table import Table
from astropy.units import Quantity, UnitsError, dimensionless_unscaled
from astropy.units.utils import quantity_asanyarray
from astropy.utils import (sharedmethod, find_current_module,
                           check_broadcast, IncompatibleShapeError, isiterable)
from astropy.utils.codegen import make_function_with_signature
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.misc import get_parameters
from astropy.nddata.utils import add_array, extract_array
from .utils import (combine_labels, make_binary_operator_eval,
                    get_inputs_and_params, _BoundingBox, _combine_equivalency_dict,
                    _ConstraintsDict)
from .parameters import (Parameter, InputParameterError,
                         param_repr_oneline, _tofloat)


__all__ = ['Model', 'FittableModel', 'Fittable1DModel', 'Fittable2DModel',
           'CompoundModel', 'fix_inputs', 'custom_model', 'ModelDefinitionError']


def _model_oper(oper, **kwargs):
    """
    Returns a function that evaluates a given Python arithmetic operator
    between two models.  The operator should be given as a string, like ``'+'``
    or ``'**'``.
    """
    return lambda left, right: CompoundModel(oper, left, right, **kwargs)


class ModelDefinitionError(TypeError):
    """Used for incorrect models definitions."""


class _ModelMeta(abc.ABCMeta):
    """
    Metaclass for Model.

    Currently just handles auto-generating the param_names list based on
    Parameter descriptors declared at the class-level of Model subclasses.
    """
    @classmethod
    def __prepare__(mcls, name, bases):
        return OrderedDict()

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

    def __new__(mcls, name, bases, members):
        # See the docstring for _is_dynamic above
        if '_is_dynamic' not in members:
            members['_is_dynamic'] = mcls._is_dynamic
        get_parameters(members)
        opermethods = [
            ('__add__', _model_oper('+')),
            ('__sub__', _model_oper('-')),
            ('__mul__', _model_oper('*')),
            ('__truediv__', _model_oper('/')),
            ('__pow__', _model_oper('**')),
            ('__or__', _model_oper('|')),
            ('__and__', _model_oper('&')),
            ('_fix_inputs', _model_oper('fix_inputs'))
        ]

        for opermethod, opercall in opermethods:
            members[opermethod] = opercall
        cls = super().__new__(mcls, name, bases, members)

        param_names = list(members['_parameters_'])

        # Need to walk each base MRO to collect all parameter names
        for base in bases:
            for tbase in base.__mro__:
                if issubclass(tbase, Model):
                    # Preserve order of definitions
                    param_names = list(tbase._parameters_) + param_names
        # Remove duplicates (arising from redefintion in subclass).
        param_names = list(dict.fromkeys(param_names))
        if cls._parameters_:
            if hasattr(cls, '_param_names'):
                # Slight kludge to support compound models, where
                # cls.param_names is a property; could be improved with a
                # little refactoring but fine for now
                cls._param_names = tuple(param_names)
            else:
                cls.param_names = tuple(param_names)

        return cls

    def __init__(cls, name, bases, members):
        super(_ModelMeta, cls).__init__(name, bases, members)
        if cls.__name__ != "CompoundModel":
            cls._create_inverse_property(members)
        cls._create_bounding_box_property(members)
        pdict = OrderedDict()
        for base in bases:
            for tbase in base.__mro__:
                if issubclass(tbase, Model):
                    for parname, val in cls._parameters_.items():
                        pdict[parname] = val
        cls._handle_special_methods(members, pdict)

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
    def _is_concrete(cls):
        """
        A class-level property that determines whether the class is a concrete
        implementation of a Model--i.e. it is not some abstract base class or
        internal implementation detail (i.e. begins with '_').
        """
        return not (cls.__name__.startswith('_') or inspect.isabstract(cls))

    def rename(cls, name=None, inputs=None, outputs=None):
        """
        Creates a copy of this model class with a new name, inputs or outputs.

        The new class is technically a subclass of the original class, so that
        instance and type checks will still work.  For example::

            >>> from astropy.modeling.models import Rotation2D
            >>> SkyRotation = Rotation2D.rename('SkyRotation')
            >>> SkyRotation
            <class 'astropy.modeling.core.SkyRotation'>
            Name: SkyRotation (Rotation2D)
            N_inputs: 2
            N_outputs: 2
            Fittable parameters: ('angle',)
            >>> issubclass(SkyRotation, Rotation2D)
            True
            >>> r = SkyRotation(90)
            >>> isinstance(r, Rotation2D)
            True
        """

        mod = find_current_module(2)
        if mod:
            modname = mod.__name__
        else:
            modname = '__main__'

        if name is None:
            name = cls.name
        if inputs is None:
            inputs = cls.inputs
        else:
            if not isinstance(inputs, tuple):
                raise TypeError("Expected 'inputs' to be a tuple of strings.")
            elif len(inputs) != len(cls.inputs):
                raise ValueError(f'{cls.name} expects {len(cls.inputs)} inputs')
        if outputs is None:
            outputs = cls.outputs
        else:
            if not isinstance(outputs, tuple):
                raise TypeError("Expected 'outputs' to be a tuple of strings.")
            elif len(outputs) != len(cls.outputs):
                raise ValueError(f'{cls.name} expects {len(cls.outputs)} outputs')
        new_cls = type(name, (cls,), {"inputs": inputs, "outputs": outputs})
        new_cls.__module__ = modname
        new_cls.__qualname__ = name

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
            except ValueError as exc:
                raise ModelDefinitionError(exc.args[0])
        else:
            sig = signature(bounding_box)
            # May be a method that only takes 'self' as an argument (like a
            # property, but the @property decorator was forgotten)
            #
            # However, if the method takes additional arguments then this is a
            # parameterized bounding box and should be callable
            if len(sig.parameters) > 1:
                bounding_box = \
                        cls._create_bounding_box_subclass(bounding_box, sig)

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

        __call__.__signature__ = sig

        return type('_{0}BoundingBox'.format(cls.name), (_BoundingBox,),
                    {'__call__': __call__})

    def _handle_special_methods(cls, members, pdict):

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

        if ('__call__' not in members and 'n_inputs' in members and
                isinstance(members['n_inputs'], int) and members['n_inputs'] > 0):

            # Don't create a custom __call__ for classes that already have one
            # explicitly defined (this includes the Model base class, and any
            # other classes that manually override __call__

            def __call__(self, *inputs, **kwargs):
                """Evaluate this model on the supplied inputs."""
                return super(cls, self).__call__(*inputs, **kwargs)

            # When called, models can take two optional keyword arguments:
            #
            # * model_set_axis, which indicates (for multi-dimensional input)
            #   which axis is used to indicate different models
            #
            # * equivalencies, a dictionary of equivalencies to be applied to
            #   the input values, where each key should correspond to one of
            #   the inputs.
            #
            # The following code creates the __call__ function with these
            # two keyword arguments.

            args = ('self',)
            kwargs = dict([('model_set_axis', None),
                           ('with_bounding_box', False),
                           ('fill_value', np.nan),
                           ('equivalencies', None),
                           ('inputs_map', None)])

            new_call = make_function_with_signature(
                __call__, args, kwargs, varargs='inputs', varkwargs='new_inputs')

            # The following makes it look like __call__
            # was defined in the class
            update_wrapper(new_call, cls)

            cls.__call__ = new_call

        if ('__init__' not in members and not inspect.isabstract(cls) and
                cls._parameters_):
            # Build list of all parameters including inherited ones

            # If *all* the parameters have default values we can make them
            # keyword arguments; otherwise they must all be positional
            # arguments
            if all(p.default is not None for p in pdict.values()):
                args = ('self',)
                kwargs = []
                for param_name, param_val in pdict.items():
                    default = param_val.default
                    unit = param_val.unit
                    # If the unit was specified in the parameter but the
                    # default is not a Quantity, attach the unit to the
                    # default.
                    if unit is not None:
                        default = Quantity(default, unit, copy=False)
                    kwargs.append((param_name, default))
            else:
                args = ('self',) + tuple(pdict.keys())
                kwargs = {}

            def __init__(self, *params, **kwargs):
                return super(cls, self).__init__(*params, **kwargs)

            new_init = make_function_with_signature(
                __init__, args, kwargs, varkwargs='kwargs')
            update_wrapper(new_init, cls)
            cls.__init__ = new_init

    # *** Arithmetic operators for creating compound models ***
    __add__ = _model_oper('+')
    __sub__ = _model_oper('-')
    __mul__ = _model_oper('*')
    __truediv__ = _model_oper('/')
    __pow__ = _model_oper('**')
    __or__ = _model_oper('|')
    __and__ = _model_oper('&')
    _fix_inputs = _model_oper('fix_inputs')

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
        parts = [super().__repr__()]

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
            return cls.name

        try:
            default_keywords = [
                ('Name', format_inheritance(cls)),
                ('N_inputs', cls.n_inputs),
                ('N_outputs', cls.n_outputs),
            ]

            if cls.param_names:
                default_keywords.append(('Fittable parameters',
                                         cls.param_names))

            for keyword, value in default_keywords + keywords:
                if value is not None:
                    parts.append('{0}: {1}'.format(keyword, value))

            return '\n'.join(parts)
        except Exception:
            # If any of the above formatting fails fall back on the basic repr
            # (this is particularly useful in debugging)
            return parts[0]


class Model(metaclass=_ModelMeta):
    """
    Base class for all models.

    This is an abstract class and should not be instantiated directly.

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
        sets. Also, ``model_set_axis=False`` can be used to tell that a given
        input should be used to evaluate all the models in the model set.

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
        A dictionary ``{parameter_name: value}`` of lower and upper bounds of
        parameters. Keys are parameter names. Values are a list or a tuple
        of length 2 giving the desired range for the parameter.

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

    n_inputs = 0
    """The number of inputs."""
    n_outputs = 0
    """ The number of outputs."""

    standard_broadcasting = True
    fittable = False
    linear = True
    _separable = None
    """ A boolean flag to indicate whether a model is separable."""
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

    # New classes can set this as a boolean value.
    # It is converted to a dictionary mapping input name to a boolean value.
    _input_units_strict = False

    # Allow dimensionless input (and corresponding output). If this is True,
    # input values to evaluate will gain the units specified in input_units. If
    # this is a dictionary then it should map input name to a bool to allow
    # dimensionless numbers for that input.
    # Only has an effect if input_units is defined.
    _input_units_allow_dimensionless = False

    # Default equivalencies to apply to input values. If set, this should be a
    # dictionary where each key is a string that corresponds to one of the
    # model inputs. Only has an effect if input_units is defined.
    input_units_equivalencies = None

    def __init__(self, *args, meta=None, name=None, **kwargs):
        super().__init__()
        self._default_inputs_outputs()
        if meta is not None:
            self.meta = meta
        self._name = name
        # add parameters to instance level by walking MRO list
        mro = self.__class__.__mro__
        for cls in mro:
            if issubclass(cls, Model):
                for parname, val in cls._parameters_.items():
                    newpar = copy.deepcopy(val)
                    newpar.model = self
                    if parname not in self.__dict__:
                        self.__dict__[parname] = newpar

        self._initialize_constraints(kwargs)
        # Remaining keyword args are either parameter values or invalid
        # Parameter values must be passed in as keyword arguments in order to
        # distinguish them
        self._initialize_parameters(args, kwargs)
        self._initialize_slices()
        self._initialize_unit_support()

        # Raise DeprecationWarning on classes with class attributes
        # ``inputs`` and ``outputs``.
        self._inputs_deprecation()

    def _inputs_deprecation(self):
        if hasattr(self.__class__, 'inputs') and isinstance(self.__class__.inputs, tuple):
            warnings.warn(
                f"""Class {self.__class__.__name__} defines class attributes ``inputs``.
                This has been deprecated in v4.0 and support will be removed in v4.1.
                Starting with v4.0 classes must define a class attribute ``n_inputs``.
                Please consult the documentation for details.
                """, AstropyDeprecationWarning)

    def _default_inputs_outputs(self):
        if self.n_inputs == 1 and self.n_outputs == 1:
            self._inputs = ("x",)
            self._outputs = ("y",)
        elif self.n_inputs == 2 and self.n_outputs == 1:
            self._inputs = ("x", "y")
            self._outputs = ("z",)
        else:
            try:
                self._inputs = tuple("x" + str(idx) for idx in range(self.n_inputs))
                self._outputs = tuple("x" + str(idx) for idx in range(self.n_outputs))
            except TypeError:
                # self.n_inputs and self.n_outputs are properties
                # This is the case when subclasses of Model do not define
                # ``n_inputs``, ``n_outputs``, ``inputs`` or ``outputs``.
                self._inputs = ()
                self._outputs = ()

    @property
    def inputs(self):
        return self._inputs

    @inputs.setter
    def inputs(self, val):
        if len(val) != self.n_inputs:
            raise ValueError(f"Expected {self.n_inputs} number of inputs, got {len(val)}.")
        self._inputs = val
        self._initialize_unit_support()

    @property
    def outputs(self):
        return self._outputs

    @outputs.setter
    def outputs(self, val):
        if len(val) != self.n_outputs:
            raise ValueError(f"Expected {self.n_outputs} number of outputs, got {len(val)}.")
        self._outputs = val

    @property
    def n_inputs(self):
        # TODO: remove the code in the ``if`` block when support
        # for models with ``inputs`` as class variables is removed.
        if hasattr(self.__class__, 'n_inputs') and isinstance(self.__class__.n_inputs, property):
            try:
                return len(self.__class__.inputs)
            except TypeError:
                try:
                    return len(self.inputs)
                except AttributeError:
                    return 0

        return self.__class__.n_inputs

    @property
    def n_outputs(self):
        # TODO: remove the code in the ``if`` block when support
        # for models with ``outputs`` as class variables is removed.
        if hasattr(self.__class__, 'n_outputs') and isinstance(self.__class__.n_outputs, property):
            try:
                return len(self.__class__.outputs)
            except TypeError:
                try:
                    return len(self.outputs)
                except AttributeError:
                    return 0

        return self.__class__.n_outputs

    def _initialize_unit_support(self):
        """
        Convert self._input_units_strict and
        self.input_units_allow_dimensionless to dictionaries
        mapping input name to a boolean value.
        """
        if isinstance(self._input_units_strict, bool):
            self._input_units_strict = {key: self._input_units_strict for
                                        key in self.inputs}

        if isinstance(self._input_units_allow_dimensionless, bool):
            self._input_units_allow_dimensionless = {key: self._input_units_allow_dimensionless
                                                     for key in self.inputs}

    @property
    def input_units_strict(self):
        """
        Enforce strict units on inputs to evaluate. If this is set to True,
        input values to evaluate will be in the exact units specified by
        input_units. If the input quantities are convertible to input_units,
        they are converted. If this is a dictionary then it should map input
        name to a bool to set strict input units for that parameter.
        """
        val = self._input_units_strict
        if isinstance(val, bool):
            return {key: val for key in self.inputs}
        return dict(zip(self.inputs, val.values()))

    @property
    def input_units_allow_dimensionless(self):
        """
        Allow dimensionless input (and corresponding output). If this is True,
        input values to evaluate will gain the units specified in input_units. If
        this is a dictionary then it should map input name to a bool to allow
        dimensionless numbers for that input.
        Only has an effect if input_units is defined.
        """

        val = self._input_units_allow_dimensionless
        if isinstance(val, bool):
            return {key: val for key in self.inputs}
        return dict(zip(self.inputs, val.values()))

    @property
    def uses_quantity(self):
        """
        True if this model has been created with `~astropy.units.Quantity`
        objects or if there are no parameters.

        This can be used to determine if this model should be evaluated with
        `~astropy.units.Quantity` or regular floats.
        """
        pisq = [isinstance(p, Quantity) for p in self._param_sets(units=True)]
        return (len(pisq) == 0) or any(pisq)

    def __repr__(self):
        return self._format_repr()

    def __str__(self):
        return self._format_str()

    def __len__(self):
        return self._n_models

    def __setattr__(self, attr, value):
        if isinstance(self, CompoundModel):
            param_names = self._param_names
        param_names = self.param_names

        if param_names is not None and attr in self.param_names:
            param = self.__dict__[attr]
            value = _tofloat(value)
            if param._validator is not None:
                param._validator(self, value)
            # check consistency with previous shape and size
            eshape = self._param_metrics[attr]['shape']
            if eshape == ():
                eshape = (1,)
            vshape = np.array(value).shape
            if vshape == ():
                vshape = (1,)
            esize = self._param_metrics[attr]['size']
            if (np.size(value) != esize or
                    _strip_ones(vshape) != _strip_ones(eshape)):
                raise InputParameterError(
                    "Value for parameter {0} does not match shape or size\n"
                    "expected by model ({1}, {2}) vs ({3}, {4})".format(
                        attr, vshape, np.size(value), eshape, esize))
            if param.unit is None:
                if isinstance(value, Quantity):
                    param._unit = value.unit
                    param.value = value.value
                else:
                    param.value = value
            else:
                if not isinstance(value, Quantity):
                    raise UnitsError(f"The '{param.name}' parameter should be given as a"
                                     " Quantity because it was originally "
                                     "initialized as a Quantity")
                param._unit = value.unit
                param.value = value.value
        else:
            if attr in ['fittable', 'linear']:
                self.__dict__[attr] = value
            else:
                super().__setattr__(attr, value)

    def __call__(self, *args, **kwargs):
        """
        Evaluate this model using the given input(s) and the parameter values
        that were specified when the model was instantiated.
        """
        new_args, kwargs = self._get_renamed_inputs_as_positional(*args, **kwargs)

        return generic_call(self, *new_args, **kwargs)

    def _get_renamed_inputs_as_positional(self, *args, **kwargs):
        def _keyword2positional(kwargs):
            # Inputs were passed as keyword (not positional) arguments.
            # Because the signature of the ``__call__`` is defined at
            # the class level, the name of the inputs cannot be changed at
            # the instance level and the old names are always present in the
            # signature of the method. In order to use the new names of the
            # inputs, the old names are taken out of ``kwargs``, the input
            # values are sorted in the order of self.inputs and passed as
            # positional arguments to ``__call__``.

            # These are the keys that are always present as keyword arguments.
            keys = ['model_set_axis', 'with_bounding_box', 'fill_value',
                    'equivalencies', 'inputs_map']

            new_inputs = {}
            # kwargs contain the names of the new inputs + ``keys``
            allkeys = list(kwargs.keys())
            # Remove the names of the new inputs from kwargs and save them
            # to a dict ``new_inputs``.
            for key in allkeys:
                if key not in keys:
                    new_inputs[key] = kwargs[key]
                    del kwargs[key]
            return new_inputs, kwargs
        n_args = len(args)

        new_inputs, kwargs = _keyword2positional(kwargs)
        n_all_args = n_args + len(new_inputs)

        if n_all_args < self.n_inputs:
            raise ValueError(f"Missing input arguments - expected {self.n_inputs}, got {n_all_args}")
        elif n_all_args > self.n_inputs:
            raise ValueError(f"Too many input arguments - expected {self.n_inputs}, got {n_all_args}")
        if n_args == 0:
            # Create positional arguments from the keyword arguments in ``new_inputs``.
            new_args = []
            for k in self.inputs:
                new_args.append(new_inputs[k])
        elif n_args != self.n_inputs:
            # Some inputs are passed as positional, others as keyword arguments.
            args = list(args)

            # Create positional arguments from the keyword arguments in ``new_inputs``.
            new_args = []
            for k in self.inputs:
                if k in new_inputs:
                    new_args.append(new_inputs[k])
                else:
                    new_args.append(args[0])
                    del args[0]
        else:
            new_args = args
        return new_args, kwargs

    # *** Properties ***
    @property
    def name(self):
        """User-provided name for this model instance."""

        return self._name

    @name.setter
    def name(self, val):
        """Assign a (new) name to this model."""

        self._name = val

    @property
    def model_set_axis(self):
        """
        The index of the model set axis--that is the axis of a parameter array
        that pertains to which model a parameter value pertains to--as
        specified when the model was initialized.

        See the documentation on :ref:`modeling-model-sets`
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
        self._array_to_parameters()

    @property
    def fixed(self):
        """
        A ``dict`` mapping parameter names to their fixed constraint.
        """

        return _ConstraintsDict(self, 'fixed')

    @property
    def bounds(self):
        """
        A ``dict`` mapping parameter names to their upper and lower bounds as
        ``(min, max)`` tuples or ``[min, max]`` lists.
        """
        return _ConstraintsDict(self, 'bounds')

    @property
    def tied(self):
        """
        A ``dict`` mapping parameter names to their tied constraint.
        """

        return _ConstraintsDict(self, 'tied')

    @property
    def eqcons(self):
        """List of parameter equality constraints."""

        return self._mconstraints['eqcons']

    @property
    def ineqcons(self):
        """List of parameter inequality constraints."""

        return self._mconstraints['ineqcons']

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
                "instance or `None` (where `None` explicitly forces the "
                "model to have no inverse.")

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
    def separable(self):
        """ A flag indicating whether a model is separable."""

        if self._separable is not None:
            return self._separable
        raise NotImplementedError(
            'The "separable" property is not defined for '
            'model {}'.format(self.__class__.__name__))

    # *** Public methods ***

    def without_units_for_data(self, **kwargs):
        """
        Return an instance of the model for which the parameter values have
        been converted to the right units for the data, then the units have
        been stripped away.

        The input and output Quantity objects should be given as keyword
        arguments.

        Notes
        -----

        This method is needed in order to be able to fit models with units in
        the parameters, since we need to temporarily strip away the units from
        the model during the fitting (which might be done by e.g. scipy
        functions).

        The units that the parameters should be converted to are not
        necessarily the units of the input data, but are derived from them.
        Model subclasses that want fitting to work in the presence of
        quantities need to define a ``_parameter_units_for_data_units`` method
        that takes the input and output units (as two dictionaries) and
        returns a dictionary giving the target units for each parameter.

        """
        model = self.copy()

        inputs_unit = {inp: getattr(kwargs[inp], 'unit', dimensionless_unscaled)
                       for inp in self.inputs if kwargs[inp] is not None}

        outputs_unit = {out: getattr(kwargs[out], 'unit', dimensionless_unscaled)
                        for out in self.outputs if kwargs[out] is not None}
        parameter_units = self._parameter_units_for_data_units(inputs_unit,
                                                               outputs_unit)
        for name, unit in parameter_units.items():
            parameter = getattr(model, name)
            if parameter.unit is not None:
                parameter.value = parameter.quantity.to(unit).value
                parameter._set_unit(None, force=True)

        if isinstance(model, CompoundModel):
            model.strip_units_from_tree()

        return model

    def strip_units_from_tree(self):
        for item in self._leaflist:
            for parname in item.param_names:
                par = getattr(item, parname)
                par._set_unit(None, force=True)

    def with_units_from_data(self, **kwargs):
        """
        Return an instance of the model which has units for which the parameter
        values are compatible with the data units specified.

        The input and output Quantity objects should be given as keyword
        arguments.

        Notes
        -----

        This method is needed in order to be able to fit models with units in
        the parameters, since we need to temporarily strip away the units from
        the model during the fitting (which might be done by e.g. scipy
        functions).

        The units that the parameters will gain are not necessarily the units
        of the input data, but are derived from them. Model subclasses that
        want fitting to work in the presence of quantities need to define a
        ``_parameter_units_for_data_units`` method that takes the input and output
        units (as two dictionaries) and returns a dictionary giving the target
        units for each parameter.
        """

        model = self.copy()
        inputs_unit = {inp: getattr(kwargs[inp], 'unit', dimensionless_unscaled)
                       for inp in self.inputs if kwargs[inp] is not None}

        outputs_unit = {out: getattr(kwargs[out], 'unit', dimensionless_unscaled)
                        for out in self.outputs if kwargs[out] is not None}

        parameter_units = self._parameter_units_for_data_units(inputs_unit,
                                                               outputs_unit)

        # We are adding units to parameters that already have a value, but we
        # don't want to convert the parameter, just add the unit directly,
        # hence the call to ``_set_unit``.
        for name, unit in parameter_units.items():
            parameter = getattr(model, name)
            parameter._set_unit(unit, force=True)

        return model

    @property
    def _has_units(self):
        # Returns True if any of the parameters have units
        for param in self.param_names:
            if getattr(self, param).unit is not None:
                return True
        else:
            return False

    @property
    def _supports_unit_fitting(self):
        # If the model has a ``_parameter_units_for_data_units`` method, this
        # indicates that we have enough information to strip the units away
        # and add them back after fitting, when fitting quantities
        return hasattr(self, '_parameter_units_for_data_units')

    @abc.abstractmethod
    def evaluate(self, *args, **kwargs):
        """Evaluate the model on some input variables."""

    def sum_of_implicit_terms(self, *args, **kwargs):
        """
        Evaluate the sum of any implicit model terms on some input variables.
        This includes any fixed terms used in evaluating a linear model that
        do not have corresponding parameters exposed to the user. The
        prototypical case is `astropy.modeling.functional_models.Shift`, which
        corresponds to a function y = a + bx, where b=1 is intrinsically fixed
        by the type of model, such that sum_of_implicit_terms(x) == x. This
        method is needed by linear fitters to correct the dependent variable
        for the implicit term(s) when solving for the remaining terms
        (ie. a = y - bx).
        """

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
        coords : array_like, optional
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
            # Assures position is at center pixel,
            # important when using add_array.
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

    @property
    def input_units(self):
        """
        This property is used to indicate what units or sets of units the
        evaluate method expects, and returns a dictionary mapping inputs to
        units (or `None` if any units are accepted).

        Model sub-classes can also use function annotations in evaluate to
        indicate valid input units, in which case this property should
        not be overridden since it will return the input units based on the
        annotations.
        """
        if hasattr(self, '_input_units'):
            return self._input_units
        elif hasattr(self.evaluate, '__annotations__'):
            annotations = self.evaluate.__annotations__.copy()
            annotations.pop('return', None)
            if annotations:
                # If there are not annotations for all inputs this will error.
                return dict((name, annotations[name]) for name in self.inputs)
        else:
            # None means any unit is accepted
            return None

    @property
    def return_units(self):
        """
        This property is used to indicate what units or sets of units the
        output of evaluate should be in, and returns a dictionary mapping
        outputs to units (or `None` if any units are accepted).

        Model sub-classes can also use function annotations in evaluate to
        indicate valid output units, in which case this property should not be
        overridden since it will return the return units based on the
        annotations.
        """
        if hasattr(self, '_return_units'):
            return self._return_units
        elif hasattr(self.evaluate, '__annotations__'):
            return self.evaluate.__annotations__.get('return', None)
        else:
            # None means any unit is accepted
            return None

    def prepare_inputs(self, *inputs, model_set_axis=None, equivalencies=None,
                       **kwargs):
        """
        This method is used in `~astropy.modeling.Model.__call__` to ensure
        that all the inputs to the model can be broadcast into compatible
        shapes (if one or both of them are input as arrays), particularly if
        there are more than one parameter sets. This also makes sure that (if
        applicable) the units of the input will be compatible with the evaluate
        method.
        """
        # When we instantiate the model class, we make sure that __call__ can
        # take the following two keyword arguments: model_set_axis and
        # equivalencies.

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

        inputs_map = kwargs.get('inputs_map', None)

        inputs = self._validate_input_units(inputs, equivalencies, inputs_map)

        # The input formatting required for single models versus a multiple
        # model set are different enough that they've been split into separate
        # subroutines
        if n_models == 1:
            return _prepare_inputs_single_model(self, params, inputs,
                                                **kwargs)
        else:
            return _prepare_inputs_model_set(self, params, inputs, n_models,
                                             model_set_axis, **kwargs)

    def _validate_input_units(self, inputs, equivalencies=None, inputs_map=None):
        inputs = list(inputs)
        name = self.name or self.__class__.__name__
        # Check that the units are correct, if applicable

        if self.input_units is not None:
            # If a leaflist is provided that means this is in the context of
            # a compound model and it is necessary to create the appropriate
            # alias for the input coordinate name for the equivalencies dict
            if inputs_map:
                edict = {}
                for mod, mapping in inputs_map:
                    if self is mod:
                        edict[mapping[0]] = equivalencies[mapping[1]]
            else:
                edict = equivalencies
            # We combine any instance-level input equivalencies with user
            # specified ones at call-time.
            input_units_equivalencies = _combine_equivalency_dict(self.inputs,
                                                                  edict,
                                                                  self.input_units_equivalencies)

            # We now iterate over the different inputs and make sure that their
            # units are consistent with those specified in input_units.
            for i in range(len(inputs)):

                input_name = self.inputs[i]
                input_unit = self.input_units.get(input_name, None)

                if input_unit is None:
                    continue

                if isinstance(inputs[i], Quantity):

                    # We check for consistency of the units with input_units,
                    # taking into account any equivalencies

                    if inputs[i].unit.is_equivalent(
                            input_unit,
                            equivalencies=input_units_equivalencies[input_name]):

                        # If equivalencies have been specified, we need to
                        # convert the input to the input units - this is
                        # because some equivalencies are non-linear, and
                        # we need to be sure that we evaluate the model in
                        # its own frame of reference. If input_units_strict
                        # is set, we also need to convert to the input units.
                        if len(input_units_equivalencies) > 0 or self.input_units_strict[input_name]:
                            inputs[i] = inputs[i].to(input_unit,
                                                     equivalencies=input_units_equivalencies[input_name])

                    else:

                        # We consider the following two cases separately so as
                        # to be able to raise more appropriate/nicer exceptions

                        if input_unit is dimensionless_unscaled:
                            raise UnitsError("{0}: Units of input '{1}', {2} ({3}),"
                                             "could not be converted to "
                                             "required dimensionless "
                                             "input".format(name,
                                                            self.inputs[i],
                                                            inputs[i].unit,
                                                            inputs[i].unit.physical_type))
                        else:
                            raise UnitsError("{0}: Units of input '{1}', {2} ({3}),"
                                             " could not be "
                                             "converted to required input"
                                             " units of {4} ({5})".format(
                                                 name,
                                                 self.inputs[i],
                                                 inputs[i].unit,
                                                 inputs[i].unit.physical_type,
                                                 input_unit,
                                                 input_unit.physical_type))
                else:

                    # If we allow dimensionless input, we add the units to the
                    # input values without conversion, otherwise we raise an
                    # exception.

                    if (not self.input_units_allow_dimensionless[input_name] and
                        input_unit is not dimensionless_unscaled and
                        input_unit is not None):
                        if np.any(inputs[i] != 0):
                            raise UnitsError("{0}: Units of input '{1}', (dimensionless), could not be "
                                             "converted to required input units of "
                                             "{2} ({3})".format(name, self.inputs[i], input_unit,
                                                                input_unit.physical_type))
        return inputs

    def _process_output_units(self, inputs, outputs):
        inputs_are_quantity = any([isinstance(i, Quantity) for i in inputs])
        if self.return_units and inputs_are_quantity:
            # We allow a non-iterable unit only if there is one output
            if self.n_outputs == 1 and not isiterable(self.return_units):
                return_units = {self.outputs[0]: self.return_units}
            else:
                return_units = self.return_units

            outputs = tuple([Quantity(out, return_units.get(out_name, None), subok=True)
                             for out, out_name in zip(outputs, self.outputs)])
        return outputs

    def prepare_outputs(self, format_info, *outputs, **kwargs):
        model_set_axis = kwargs.get('model_set_axis', None)

        if len(self) == 1:
            return _prepare_outputs_single_model(outputs, format_info)
        else:
            return _prepare_outputs_model_set(self, outputs, format_info, model_set_axis)

    def copy(self):
        """
        Return a copy of this model.

        Uses a deep copy so that all model attributes, including parameter
        values, are copied as well.
        """

        return copy.deepcopy(self)

    def deepcopy(self):
        """
        Return a deep copy of this model.

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

    @property
    def n_submodels(self):
        """
        Return the number of components in a single model, which is
        obviously 1.
        """
        return 1

    def _initialize_constraints(self, kwargs):
        """
        Pop parameter constraint values off the keyword arguments passed to
        `Model.__init__` and store them in private instance attributes.
        """

        # Pop any constraints off the keyword arguments
        for constraint in self.parameter_constraints:
            values = kwargs.pop(constraint, {})
            for ckey, cvalue in values.items():
                param = getattr(self, ckey)
                setattr(param, constraint, cvalue)
        self._mconstraints = {}
        for constraint in self.model_constraints:
            values = kwargs.pop(constraint, [])
            self._mconstraints[constraint] = values

    def _initialize_parameters(self, args, kwargs):
        """
        Initialize the _parameters array that stores raw parameter values for
        all parameter sets for use with vectorized fitting algorithms; on
        FittableModels the _param_name attributes actually just reference
        slices of this array.
        """
        n_models = kwargs.pop('n_models', None)

        if not (n_models is None or
                (isinstance(n_models, (int, np.integer)) and n_models >= 1)):
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
        params = set()
        if len(args) > len(self.param_names):
            raise TypeError(
                "{0}.__init__() takes at most {1} positional arguments ({2} "
                "given)".format(self.__class__.__name__, len(self.param_names),
                                len(args)))

        self._model_set_axis = model_set_axis
        self._param_metrics = defaultdict(dict)

        for idx, arg in enumerate(args):
            if arg is None:
                # A value of None implies using the default value, if exists
                continue
            # We use quantity_asanyarray here instead of np.asanyarray because
            # if any of the arguments are quantities, we need to return a
            # Quantity object not a plain Numpy array.
            param_name = self.param_names[idx]
            params.add(param_name)
            if not isinstance(arg, Parameter):
                value = quantity_asanyarray(arg, dtype=float)
            else:
                value = arg
            self._initialize_parameter_value(param_name, value)

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
                # We use quantity_asanyarray here instead of np.asanyarray
                # because if any of the arguments are quantities, we need
                # to return a Quantity object not a plain Numpy array.
                value = quantity_asanyarray(value, dtype=float)
                params.add(param_name)
                self._initialize_parameter_value(param_name, value)
        # Now deal with case where param_name is not supplied by args or kwargs
        for param_name in self.param_names:
            if param_name not in params:
                self._initialize_parameter_value(param_name, None)

        if kwargs:
            # If any keyword arguments were left over at this point they are
            # invalid--the base class should only be passed the parameter
            # values, constraints, and param_dim
            for kwarg in kwargs:
                # Just raise an error on the first unrecognized argument
                raise TypeError(
                    '{0}.__init__() got an unrecognized parameter '
                    '{1!r}'.format(self.__class__.__name__, kwarg))

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

            for name in self.param_names:
                value = getattr(self, name)
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

            self._check_param_broadcast(max_ndim)
        else:
            if n_models is None:
                n_models = 1

            self._check_param_broadcast(None)

        self._n_models = n_models
        # now validate parameters
        for name in params:
            param = getattr(self, name)
            if param._validator is not None:
                param._validator(self, param.value)

    def _initialize_parameter_value(self, param_name, value):
        """Mostly deals with consistency checks and determining unit issues."""
        if isinstance(value, Parameter):
            self.__dict__[param_name] = value
            return
        param = getattr(self, param_name)
        # Use default if value is not provided
        if value is None:
            default = param.default
            if default is None:
                # No value was supplied for the parameter and the
                # parameter does not have a default, therefore the model
                # is underspecified
                raise TypeError("{0}.__init__() requires a value for parameter "
                                "{1!r}".format(self.__class__.__name__, param_name))
            value = default
            unit = param.unit
        else:
            if isinstance(value, Quantity):
                unit = value.unit
                value = value.value
            else:
                unit = None
        if unit is None and param.unit is not None:
            raise InputParameterError(
                "{0}.__init__() requires a Quantity for parameter "
                "{1!r}".format(self.__class__.__name__, param_name))
        param._unit = unit
        param.internal_unit = None
        if param._setter is not None:
            if unit is not None:
                _val = param._setter(value * unit)
            else:
                _val = param._setter(value)
            if isinstance(_val, Quantity):
                param.internal_unit = _val.unit
                param._internal_value = np.array(_val.value)
            else:
                param.internal_unit = None
                param._internal_value = np.array(_val)
        else:
            param._value = np.array(value)

    def _initialize_slices(self):

        param_metrics = self._param_metrics
        total_size = 0

        for name in self.param_names:
            param = getattr(self, name)
            value = param.value
            param_size = np.size(value)
            param_shape = np.shape(value)
            param_slice = slice(total_size, total_size + param_size)
            param_metrics[name]['slice'] = param_slice
            param_metrics[name]['shape'] = param_shape
            param_metrics[name]['size'] = param_size
            total_size += param_size
        self._parameters = np.empty(total_size, dtype=np.float64)

    def _parameters_to_array(self):
        # Now set the parameter values (this will also fill
        # self._parameters)
        param_metrics = self._param_metrics
        for name in self.param_names:
            param = getattr(self, name)
            value = param.value
            if not isinstance(value, np.ndarray):
                value = np.array([value])
            self._parameters[param_metrics[name]['slice']] = value.ravel()

        # Finally validate all the parameters; we do this last so that
        # validators that depend on one of the other parameters' values will
        # work

    def _array_to_parameters(self):
        param_metrics = self._param_metrics
        for name in self.param_names:
            param = getattr(self, name)
            value = self._parameters[param_metrics[name]['slice']]
            value.shape = param_metrics[name]['shape']
            param.value = value

    def _check_param_broadcast(self, max_ndim):
        """
        This subroutine checks that all parameter arrays can be broadcast
        against each other, and determines the shapes parameters must have in
        order to broadcast correctly.

        If model_set_axis is None this merely checks that the parameters
        broadcast and returns an empty dict if so.  This mode is only used for
        single model sets.
        """
        all_shapes = []
        model_set_axis = self._model_set_axis

        for name in self.param_names:
            param = getattr(self, name)
            value = param.value
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
            param_a = self.param_names[shape_a_idx]
            param_b = self.param_names[shape_b_idx]

            raise InputParameterError(
                "Parameter {0!r} of shape {1!r} cannot be broadcast with "
                "parameter {2!r} of shape {3!r}.  All parameter arrays "
                "must have shapes that are mutually compatible according "
                "to the broadcasting rules.".format(param_a, shape_a,
                                                    param_b, shape_b))

    def _param_sets(self, raw=False, units=False):
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

        values = []
        shapes = []
        for name in self.param_names:
            param = getattr(self, name)

            if raw and param._setter:
                value = param._internal_value
            else:
                value = param.value

            broadcast_shape = self._param_metrics[name].get('broadcast_shape')
            if broadcast_shape is not None:
                value = value.reshape(broadcast_shape)

            shapes.append(np.shape(value))

            if len(self) == 1:
                # Add a single param set axis to the parameter's value (thus
                # converting scalars to shape (1,) array values) for
                # consistency
                value = np.array([value])

            if units:
                if raw and param.internal_unit is not None:
                    unit = param.internal_unit
                else:
                    unit = param.unit
                if unit is not None:
                    value = Quantity(value, unit)

            values.append(value)

        if len(set(shapes)) != 1 or units:
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

        parts = [repr(a) for a in args]

        parts.extend(
            "{0}={1}".format(name,
                             param_repr_oneline(getattr(self, name)))
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
            # Set units on the columns
            for name in self.param_names:
                param_table[name].unit = getattr(self, name).unit
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
    n_inputs = 1
    n_outputs = 1
    _separable = True


class Fittable2DModel(FittableModel):
    """
    Base class for two-dimensional fittable models.

    This class provides an easier interface to defining new models.
    Examples can be found in `astropy.modeling.functional_models`.
    """

    n_inputs = 2
    n_outputs = 1


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


BINARY_OPERATORS = {
    '+': _make_arithmetic_operator(operator.add),
    '-': _make_arithmetic_operator(operator.sub),
    '*': _make_arithmetic_operator(operator.mul),
    '/': _make_arithmetic_operator(operator.truediv),
    '**': _make_arithmetic_operator(operator.pow),
    '|': _composition_operator,
    '&': _join_operator
}

SPECIAL_OPERATORS = {}


def _add_special_operator(sop_name, sop):
    SPECIAL_OPERATORS[sop_name] = sop


class CompoundModel(Model):
    '''
    Base class for compound models.

    While it can be used directly, the recommended way
    to combine models is through the model operators.
    '''

    def __init__(self, op, left, right, name=None, inverse=None):
        self.__dict__['_param_names'] = None
        self._n_submodels = None
        self.op = op
        self.left = left
        self.right = right
        self._bounding_box = None
        self._user_bounding_box = None
        self._leaflist = None
        self._tdict = None
        self._parameters = None
        self._parameters_ = None
        self._param_metrics = None
        self._has_inverse = False  # may be set to True in following code
        if inverse:
            self.inverse = inverse

        if op != 'fix_inputs' and len(left) != len(right):
            raise ValueError(
                'Both operands must have equal values for n_models')
        self._n_models = len(left)

        if op != 'fix_inputs' and ((left.model_set_axis != right.model_set_axis)
                                   or left.model_set_axis):  # not False and not 0
            raise ValueError("model_set_axis must be False or 0 and consistent for operands")
        self._model_set_axis = left.model_set_axis

        if op in ['+', '-', '*', '/', '**'] or op in SPECIAL_OPERATORS:
            if (left.n_inputs != right.n_inputs) or \
               (left.n_outputs != right.n_outputs):
                raise ModelDefinitionError(
                    'Both operands must match numbers of inputs and outputs')
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
                inv = CompoundModel('&',
                                    self.left.inverse,
                                    self.right.inverse,
                                    inverse=self)
                if left._user_inverse is not None or right._user_inverse is not None:
                    self._user_inverse = inv
                    if self.inverse._has_inverse:
                        del self._user_inverse._user_inverse
                        self.inverse._has_inverse = False
                else:
                    self._inverse = inv
        elif op == '|':
            if left.n_outputs != right.n_inputs:
                raise ModelDefinitionError(
                    "Unsupported operands for |: {0} (n_inputs={1}, "
                    "n_outputs={2}) and {3} (n_inputs={4}, n_outputs={5}); "
                    "n_outputs for the left-hand model must match n_inputs "
                    "for the right-hand model.".format(
                        left.name, left.n_inputs, left.n_outputs, right.name,
                        right.n_inputs, right.n_outputs))

            self.n_inputs = left.n_inputs
            self.n_outputs = right.n_outputs
            self.inputs = left.inputs
            self.outputs = right.outputs
            if inverse is None and self.both_inverses_exist():
                self._has_inverse = True
                inv = CompoundModel('|',
                                    self.right.inverse,
                                    self.left.inverse,
                                    inverse=self)
                if left._user_inverse is not None or right._user_inverse is not None:
                    self._user_inverse = inv
                    if self.inverse._has_inverse:
                        del self._user_inverse._user_inverse
                    self.inverse._has_inverse = False
                else:
                    self._inverse = inv
        elif op == 'fix_inputs':
            if not isinstance(left, Model):
                raise ValueError('First argument to "fix_inputs" must be an instance of an astropy Model.')
            if not isinstance(right, dict):
                raise ValueError('Expected a dictionary for second argument of "fix_inputs".')

            # Dict keys must match either possible indices
            # for model on left side, or names for inputs.
            self.n_inputs = left.n_inputs - len(right)
            # Assign directly to the private attribute (instead of using the setter)
            # to avoid asserting the new number of outputs matches the old one.
            self._outputs = left.outputs
            self.n_outputs = left.n_outputs
            newinputs = list(left.inputs)
            keys = right.keys()
            input_ind = []
            for key in keys:
                if isinstance(key, int):
                    if key >= left.n_inputs or key < 0:
                        raise ValueError(
                            'Substitution key integer value '
                            'not among possible input choices.')
                    if key in input_ind:
                        raise ValueError("Duplicate specification of "
                                         "same input (index/name).")
                    input_ind.append(key)
                elif isinstance(key, str):
                    if key not in left.inputs:
                        raise ValueError(
                            'Substitution key string not among possible '
                            'input choices.')
                    # Check to see it doesn't match positional
                    # specification.
                    ind = left.inputs.index(key)
                    if ind in input_ind:
                        raise ValueError("Duplicate specification of "
                                         "same input (index/name).")
                    input_ind.append(ind)
            # Remove substituted inputs
            input_ind.sort()
            input_ind.reverse()
            for ind in input_ind:
                del newinputs[ind]
            self.inputs = tuple(newinputs)
            # Now check to see if the input model has bounding_box defined.
            # If so, remove the appropriate dimensions and set it for this
            # instance.
            try:
                bounding_box = self.left.bounding_box
                self._fix_input_bounding_box(input_ind)
            except NotImplementedError:
                pass

        else:
            raise ModelDefinitionError('Illegal operator: ', self.op)
        self.name = name
        self._fittable = None
        self.fit_deriv = None
        self.col_fit_deriv = None
        if op in ('|', '+', '-'):
            self.linear = left.linear and right.linear
        else:
            self.linear = False
        self.eqcons = []
        self.ineqcons = []
        self._map_parameters()

    def evaluate(self, *args, **kwargs):
        pass

    @property
    def n_submodels(self):
        if self._leaflist is None:
            self._make_leaflist()
        return len(self._leaflist)

    @property
    def submodel_names(self):
        """ Return the names of submodels in a ``CompoundModel``."""
        if self._leaflist is None:
            self._make_leaflist()
        names = [item.name for item in self._leaflist]
        nonecount = 0
        newnames = []
        for item in names:
            if item is None:
                newnames.append('None_{}'.format(nonecount))
                nonecount += 1
            else:
                newnames.append(item)
        return tuple(newnames)

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
        # Turn any keyword arguments into positional arguments.
        args, kw = self._get_renamed_inputs_as_positional(*args, **kw)

        # If equivalencies are provided, necessary to map parameters and pass
        # the leaflist as a keyword input for use by model evaluation so that
        # the compound model input names can be matched to the model input
        # names.

        if 'equivalencies' in kw:
            # Restructure to be useful for the individual model lookup
            kw['inputs_map'] = [(value[0], (value[1], key)) for
                                key, value in self.inputs_map().items()]
        with_bbox = kw.pop('with_bounding_box', False)
        fill_value = kw.pop('fill_value', np.nan)
        # Use of bounding box for compound models requires special treatment
        # in selecting only valid inputs to pass along to constituent models.
        bbox = get_bounding_box(self)
        if with_bbox and bbox is not None:
            # first check inputs are consistent in shape
            input_shape = _validate_input_shapes(args, (), self._n_models,
                                                 self.model_set_axis, self.standard_broadcasting)
            vinputs, valid_ind, allout = prepare_bounding_box_inputs(self, input_shape, args, bbox)
            if not allout:
                valid_result = self._evaluate(*vinputs, **kw)
                if self.n_outputs == 1:
                    valid_result = [valid_result]
                outputs = prepare_bounding_box_outputs(valid_result, valid_ind,
                                                       input_shape, fill_value)
            else:
                outputs = [np.zeros(input_shape) + fill_value for i in range(self.n_outputs)]
            if self.n_outputs == 1:
                return outputs[0]
            return outputs
        else:
            return self._evaluate(*args, **kw)

    def _evaluate(self, *args, **kw):
        op = self.op
        if op != 'fix_inputs':
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
            elif op in SPECIAL_OPERATORS:
                return binary_operation(SPECIAL_OPERATORS[op], leftval, rightval)
            else:

                raise ModelDefinitionError('Unrecognized operator {op}')
        else:
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
                                         "positional value supplied.")
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
                # subindsorted, subvalsorted = list(zip(*subargs))
            # The substitutions must be inserted in order
            for ind, val in subargs:
                newargs.insert(ind, val)
            return self.left(*newargs, **kw)

    @property
    def param_names(self):
        """ An ordered list of parameter names."""
        return self._param_names

    def _make_leaflist(self):
        tdict = {}
        leaflist = []
        make_subtree_dict(self, '', tdict, leaflist)
        self._leaflist = leaflist
        self._tdict = tdict

    def __getattr__(self, name):
        """
        If someone accesses an attribute not already defined, map the
        parameters, and then see if the requested attribute is one of
        the parameters
        """
        # The following test is needed to avoid infinite recursion
        # caused by deepcopy. There may be other such cases discovered.
        if name == '__setstate__':
            raise AttributeError
        if name in self._param_names:
            return self.__dict__[name]
        else:
            raise AttributeError('Attribute "{}" not found'.format(name))

    def __getitem__(self, index):
        if self._leaflist is None:
            self._make_leaflist()
        leaflist = self._leaflist
        tdict = self._tdict
        if isinstance(index, slice):
            if index.step:
                raise ValueError('Steps in slices not supported '
                                 'for compound models')
            if index.start is not None:
                if isinstance(index.start, str):
                    start = self._str_index_to_int(index.start)
                else:
                    start = index.start
            else:
                start = 0
            if index.stop is not None:
                if isinstance(index.stop, str):
                    stop = self._str_index_to_int(index.stop)
                else:
                    stop = index.stop - 1
            else:
                stop = len(leaflist) - 1
            if index.stop == 0:
                raise ValueError("Slice endpoint cannot be 0")
            if start < 0:
                start = len(leaflist) + start
            if stop < 0:
                stop = len(leaflist) + stop
            # now search for matching node:
            if stop == start:  # only single value, get leaf instead in code below
                index = start
            else:
                for key in tdict:
                    node, leftind, rightind = tdict[key]
                    if leftind == start and rightind == stop:
                        return node
                raise IndexError("No appropriate subtree matches slice")
        if isinstance(index, type(0)):
            return leaflist[index]
        elif isinstance(index, type('')):
            return leaflist[self._str_index_to_int(index)]
        else:
            raise TypeError('index must be integer, slice, or model name string')

    def _str_index_to_int(self, str_index):
        # Search through leaflist for item with that name
        found = []
        for nleaf, leaf in enumerate(self._leaflist):
            if getattr(leaf, 'name', None) == str_index:
                found.append(nleaf)
        if len(found) == 0:
            raise IndexError("No component with name '{}' found".format(str_index))
        if len(found) > 1:
            raise IndexError("Multiple components found using '{}' as name\n"
                             "at indices {}".format(str_index, found))
        return found[0]

    @property
    def n_inputs(self):
        """ The number of inputs of a model."""
        return self._n_inputs

    @n_inputs.setter
    def n_inputs(self, value):
        self._n_inputs = value

    @property
    def n_outputs(self):
        """ The number of outputs of a model."""
        return self._n_outputs

    @n_outputs.setter
    def n_outputs(self, value):
        self._n_outputs = value

    @property
    def eqcons(self):
        return self._eqcons

    @eqcons.setter
    def eqcons(self, value):
        self._eqcons = value

    @property
    def ineqcons(self):
        return self._eqcons

    @ineqcons.setter
    def ineqcons(self, value):
        self._eqcons = value

    def traverse_postorder(self, include_operator=False):
        """ Postorder traversal of the CompoundModel tree."""
        res = []
        if isinstance(self.left, CompoundModel):
            res = res + self.left.traverse_postorder(include_operator)
        else:
            res = res + [self.left]
        if isinstance(self.right, CompoundModel):
            res = res + self.right.traverse_postorder(include_operator)
        else:
            res = res + [self.right]
        if include_operator:
            res.append(self.op)
        else:
            res.append(self)
        return res

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

    def _format_components(self):
        if self._parameters_ is None:
            self._map_parameters()
        return '\n\n'.join('[{0}]: {1!r}'.format(idx, m)
                           for idx, m in enumerate(self._leaflist))

    def __str__(self):
        expression = self._format_expression()
        components = self._format_components()
        keywords = [
            ('Expression', expression),
            ('Components', '\n' + indent(components))
        ]
        return super()._format_str(keywords=keywords)

    def rename(self, name):
        self.name = name
        return self

    @property
    def isleaf(self):
        return False

    def has_inverse(self):
        return self._has_inverse

    @property
    def inverse(self):
        if self.has_inverse():
            if self._user_inverse is not None:
                return self._user_inverse
            return self._inverse
        else:
            raise NotImplementedError("Inverse function not provided")

    @inverse.setter
    def inverse(self, value):
        if not isinstance(value, Model):
            raise ValueError("Attempt to assign non model to inverse")
        self._user_inverse = value
        self._has_inverse = True

    @inverse.deleter
    def inverse(self):
        self._has_inverse = False
        self._user_inverse = None

    @property
    def fittable(self):
        """ Set the fittable attribute on a compound model."""
        if self._fittable is None:
            if self._leaflist is None:
                self._map_parameters()
            self._fittable = all(m.fittable for m in self._leaflist)
        return self._fittable

    __add__ = _model_oper('+')
    __sub__ = _model_oper('-')
    __mul__ = _model_oper('*')
    __truediv__ = _model_oper('/')
    __pow__ = _model_oper('**')
    __or__ = _model_oper('|')
    __and__ = _model_oper('&')

    def _map_parameters(self):
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
        param_map = {}
        self._param_names = []
        for lindex, leaf in enumerate(self._leaflist):
            if not isinstance(leaf, dict):
                for param_name in leaf.param_names:
                    param = getattr(leaf, param_name)
                    new_param_name = "{}_{}".format(param_name, lindex)
                    self.__dict__[new_param_name] = param
                    self._parameters_[new_param_name] = param
                    self._param_names.append(new_param_name)
                    param_map[new_param_name] = (lindex, param_name)
        self._param_metrics = {}
        self._param_map = param_map
        self._param_map_inverse = dict((v, k) for k, v in param_map.items())
        self._initialize_slices()
        self._param_names = tuple(self._param_names)

    def _initialize_slices(self):
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
            param_metrics[name]['size'] = param_size
            total_size += param_size
        self._parameters = np.empty(total_size, dtype=np.float64)

    @staticmethod
    def _recursive_lookup(branch, adict, key):
        if isinstance(branch, CompoundModel):
            return adict[key]
        return branch, key

    def inputs_map(self):
        """
        Map the names of the inputs to this ExpressionTree to the inputs to the leaf models.
        """
        inputs_map = {}
        if not isinstance(self.op, str):  # If we don't have an operator the mapping is trivial
            return {inp: (self, inp) for inp in self.inputs}

        elif self.op == '|':
            if isinstance(self.left, CompoundModel):
                l_inputs_map = self.left.inputs_map()
            for inp in self.inputs:
                if isinstance(self.left, CompoundModel):
                    inputs_map[inp] = l_inputs_map[inp]
                else:
                    inputs_map[inp] = self.left, inp
        elif self.op == '&':
            if isinstance(self.left, CompoundModel):
                l_inputs_map = self.left.inputs_map()
            if isinstance(self.right, CompoundModel):
                r_inputs_map = self.right.inputs_map()
            for i, inp in enumerate(self.inputs):
                if i < len(self.left.inputs):  # Get from left
                    if isinstance(self.left, CompoundModel):
                        inputs_map[inp] = l_inputs_map[self.left.inputs[i]]
                    else:
                        inputs_map[inp] = self.left, self.left.inputs[i]
                else:  # Get from right
                    if isinstance(self.right, CompoundModel):
                        inputs_map[inp] = r_inputs_map[self.right.inputs[i - len(self.left.inputs)]]
                    else:
                        inputs_map[inp] = self.right, self.right.inputs[i - len(self.left.inputs)]
        elif self.op == 'fix_inputs':
            fixed_ind = list(self.right.keys())
            ind = [list(self.left.inputs).index(i) if isinstance(i, str) else i for i in fixed_ind]
            inp_ind = list(range(self.left.n_inputs))
            for i in ind:
                inp_ind.remove(i)
            for i in inp_ind:
                inputs_map[self.left.inputs[i]] = self.left, self.left.inputs[i]
        else:
            if isinstance(self.left, CompoundModel):
                l_inputs_map = self.left.inputs_map()
            for inp in self.left.inputs:
                if isinstance(self.left, CompoundModel):
                    inputs_map[inp] = l_inputs_map[inp]
                else:
                    inputs_map[inp] = self.left, inp
        return inputs_map

    def _parameter_units_for_data_units(self, input_units, output_units):
        if self._leaflist is None:
            self._map_parameters()
        units_for_data = {}
        for imodel, model in enumerate(self._leaflist):
            units_for_data_leaf = model._parameter_units_for_data_units(input_units, output_units)
            for param_leaf in units_for_data_leaf:
                param = self._param_map_inverse[(imodel, param_leaf)]
                units_for_data[param] = units_for_data_leaf[param_leaf]
        return units_for_data

    @property
    def input_units(self):
        inputs_map = self.inputs_map()
        input_units_dict = {key: inputs_map[key][0].input_units[orig_key]
                            for key, (mod, orig_key) in inputs_map.items()
                            if inputs_map[key][0].input_units is not None}
        if input_units_dict:
            return input_units_dict
        return None

    @property
    def input_units_equivalencies(self):
        inputs_map = self.inputs_map()
        return {key: inputs_map[key][0].input_units_equivalencies[orig_key]
                for key, (mod, orig_key) in inputs_map.items()
                if inputs_map[key][0].input_units_equivalencies is not None}

    @property
    def input_units_allow_dimensionless(self):
        inputs_map = self.inputs_map()
        return {key: inputs_map[key][0].input_units_allow_dimensionless[orig_key]
                for key, (mod, orig_key) in inputs_map.items()}

    @property
    def input_units_strict(self):
        inputs_map = self.inputs_map()
        return {key: inputs_map[key][0].input_units_strict[orig_key]
                for key, (mod, orig_key) in inputs_map.items()}

    @property
    def return_units(self):
        outputs_map = self.outputs_map()
        return {key: outputs_map[key][0].return_units[orig_key]
                for key, (mod, orig_key) in outputs_map.items()
                if outputs_map[key][0].return_units is not None}

    def outputs_map(self):
        """
        Map the names of the outputs to this ExpressionTree to the outputs to the leaf models.
        """
        outputs_map = {}
        if not isinstance(self.op, str):  # If we don't have an operator the mapping is trivial
            return {out: (self, out) for out in self.outputs}

        elif self.op == '|':
            if isinstance(self.right, CompoundModel):
                r_outputs_map = self.right.outputs_map()
            for out in self.outputs:
                if isinstance(self.right, CompoundModel):
                    outputs_map[out] = r_outputs_map[out]
                else:
                    outputs_map[out] = self.right, out

        elif self.op == '&':
            if isinstance(self.left, CompoundModel):
                l_outputs_map = self.left.outputs_map()
            if isinstance(self.right, CompoundModel):
                r_outputs_map = self.right.outputs_map()
            for i, out in enumerate(self.outputs):
                if i < len(self.left.outputs):  # Get from left
                    if isinstance(self.left, CompoundModel):
                        outputs_map[out] = l_outputs_map[self.left.outputs[i]]
                    else:
                        outputs_map[out] = self.left, self.left.outputs[i]
                else:  # Get from right
                    if isinstance(self.right, CompoundModel):
                        outputs_map[out] = r_outputs_map[self.right.outputs[i - len(self.left.outputs)]]
                    else:
                        outputs_map[out] = self.right, self.right.outputs[i - len(self.left.outputs)]
        elif self.op == 'fix_inputs':
            return self.left.outputs_map()
        else:
            if isinstance(self.left, CompoundModel):
                l_outputs_map = self.left.outputs_map()
            for out in self.left.outputs:
                if isinstance(self.left, CompoundModel):
                    outputs_map[out] = l_outputs_map()[out]
                else:
                    outputs_map[out] = self.left, out
        return outputs_map

    def _fix_input_bounding_box(self, input_ind):
        """
        If the ``fix_inputs`` operator is used and the model it is applied to
        has a bounding box definition, delete the corresponding inputs from
        that bounding box. This method presumes the bounding_box is not None.
        This also presumes that the list of input indices to remove (i.e.,
        input_ind has already been put in reverse sorted order).
        """
        bounding_box = list(self.left.bounding_box)
        for ind in input_ind:
            del bounding_box[ind]
        if self.n_inputs == 1:
            bounding_box = bounding_box[0]
        self.bounding_box = bounding_box

    @property
    def has_user_bounding_box(self):
        """
        A flag indicating whether or not a custom bounding_box has been
        assigned to this model by a user, via assignment to
        ``model.bounding_box``.
        """

        return self._user_bounding_box is not None

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
        coords : array_like, optional
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
    return binoperator(left, right)


def get_ops(tree, opset):
    """
    Recursive function to collect operators used.
    """
    if isinstance(tree, CompoundModel):
        opset.add(tree.op)
        get_ops(tree.left, opset)
        get_ops(tree.right, opset)
    else:
        return


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


_ORDER_OF_OPERATORS = [('fix_inputs',), ('|',), ('&',), ('+', '-'), ('*', '/'), ('**',)]
OPERATOR_PRECEDENCE = {}
for idx, ops in enumerate(_ORDER_OF_OPERATORS):
    for op in ops:
        OPERATOR_PRECEDENCE[op] = idx
del idx, op, ops


def fix_inputs(modelinstance, values):
    """
    This function creates a compound model with one or more of the input
    values of the input model assigned fixed values (scalar or array).

    Parameters
    ----------
    modelinstance : Model instance. This is the model that one or more of the
        model input values will be fixed to some constant value.
    values : A dictionary where the key identifies which input to fix
        and its value is the value to fix it at. The key may either be the
        name of the input or a number reflecting its order in the inputs.

    Examples
    --------

    >>> from astropy.modeling.models import Gaussian2D
    >>> g = Gaussian2D(1, 2, 3, 4, 5)
    >>> gv = fix_inputs(g, {0: 2.5})

    Results in a 1D function equivalent to Gaussian2D(1, 2, 3, 4, 5)(x=2.5, y)
    """
    return CompoundModel('fix_inputs', modelinstance, values)


def custom_model(*args, fit_deriv=None, **kwargs):
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

    if kwargs:
        warnings.warn(
            "Function received unexpected arguments ({}) these "
            "are ignored but will raise an Exception in the "
            "future.".format(list(kwargs)),
            AstropyDeprecationWarning)

    if len(args) == 1 and callable(args[0]):
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

    if not callable(func):
        raise ModelDefinitionError(
            "func is not callable; it must be a function or other callable "
            "object")

    if fit_deriv is not None and not callable(fit_deriv):
        raise ModelDefinitionError(
            "fit_deriv not callable; it must be a function or other "
            "callable object")

    model_name = func.__name__

    inputs, params = get_inputs_and_params(func)

    if (fit_deriv is not None and
            len(fit_deriv.__defaults__) != len(params)):
        raise ModelDefinitionError("derivative function should accept "
                                   "same number of parameters as func.")

    # TODO: Maybe have a clever scheme for default output name?
    if inputs:
        output_names = (inputs[0].name,)
    else:
        output_names = ('x',)

    params = OrderedDict((param.name, Parameter(param.name,
                                                default=param.default)) for param in params)

    mod = find_current_module(2)
    if mod:
        modname = mod.__name__
    else:
        modname = '__main__'

    members = OrderedDict([
        ('__module__', str(modname)),
        ('__doc__', func.__doc__),
        ('n_inputs', len(inputs)),
        # tuple(x.name for x in inputs)),
        ('n_outputs', len(output_names)),
        ('evaluate', staticmethod(func))])

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
    coords : array_like, optional
        Coordinate arrays mapping to ``arr``, such that
        ``arr[coords] == arr``.

    Returns
    -------
    array : `numpy.ndarray`
        The model evaluated on the input ``arr`` or a new array from
        ``coords``.
        If ``arr`` and ``coords`` are both `None`, the returned array is
        limited to the `Model.bounding_box` limits. If
        `Model.bounding_box` is `None`, ``arr`` or ``coords`` must be passed.

    Examples
    --------
    :ref:`bounding-boxes`
    """

    bbox = model.bounding_box

    if (coords is None) & (arr is None) & (bbox is None):
        raise ValueError('If no bounding_box is set,'
                         'coords or arr must be input.')

    # for consistent indexing
    if model.n_inputs == 1:
        if coords is not None:
            coords = [coords]
        if bbox is not None:
            bbox = [bbox]

    if arr is not None:
        arr = arr.copy()
        # Check dimensions match model
        if arr.ndim != model.n_inputs:
            raise ValueError('number of array dimensions inconsistent with '
                             'number of model inputs.')
    if coords is not None:
        # Check dimensions match arr and model
        coords = np.array(coords)
        if len(coords) != model.n_inputs:
            raise ValueError('coordinate length inconsistent with the number '
                             'of model inputs.')
        if arr is not None:
            if coords[0].shape != arr.shape:
                raise ValueError('coordinate shape inconsistent with the '
                                 'array shape.')
        else:
            arr = np.zeros(coords[0].shape)

    if bbox is not None:
        # assures position is at center pixel, important when using add_array
        pd = pos, delta = np.array([(np.mean(bb), np.ceil((bb[1] - bb[0]) / 2))
                                    for bb in bbox]).astype(int).T

        if coords is not None:
            sub_shape = tuple(delta * 2 + 1)
            sub_coords = np.array([extract_array(c, sub_shape, pos)
                                   for c in coords])
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
        # Extend the broadcasts list to include shapes for all outputs
        extra_outputs = model.n_outputs - model.n_inputs
        if not broadcasts:
            # If there were no inputs then the broadcasts list is empty
            # just add a None since there is no broadcasting of outputs and
            # inputs necessary (see _prepare_outputs_single_model)
            broadcasts.append(None)
        broadcasts.extend([broadcasts[0]] * extra_outputs)

    return inputs, (broadcasts,)


def _prepare_outputs_single_model(outputs, format_info):
    broadcasts = format_info[0]
    outputs = list(outputs)
    for idx, output in enumerate(outputs):
        broadcast_shape = broadcasts[idx]
        if broadcast_shape is not None:
            if not broadcast_shape:
                # Shape is (), i.e. a scalar should be returned
                outputs[idx] = output.item()
            else:
                outputs[idx] = output.reshape(broadcast_shape)

    return tuple(outputs)


def _prepare_inputs_model_set(model, params, inputs, n_models, model_set_axis_input,
                              **kwargs):
    reshaped = []
    pivots = []

    model_set_axis_param = model.model_set_axis  # needed to reshape param
    for idx, _input in enumerate(inputs):
        max_param_shape = ()
        if n_models > 1 and model_set_axis_input is not False:
            # Use the shape of the input *excluding* the model axis
            input_shape = (_input.shape[:model_set_axis_input] +
                           _input.shape[model_set_axis_input + 1:])
        else:
            input_shape = _input.shape

        for param in params:
            try:
                check_broadcast(input_shape,
                                remove_axes_from_shape(param.shape,
                                                       model_set_axis_param))
            except IncompatibleShapeError:
                raise ValueError(
                    "Model input argument {0!r} of shape {1!r} cannot be "
                    "broadcast with parameter {2!r} of shape "
                    "{3!r}.".format(model.inputs[idx], input_shape,
                                    param.name,
                                    remove_axes_from_shape(param.shape,
                                                           model_set_axis_param)))

            if len(param.shape) - 1 > len(max_param_shape):
                max_param_shape = remove_axes_from_shape(param.shape,
                                                         model_set_axis_param)

        # We've now determined that, excluding the model_set_axis, the
        # input can broadcast with all the parameters
        input_ndim = len(input_shape)
        if model_set_axis_input is False:
            if len(max_param_shape) > input_ndim:
                # Just needs to prepend new axes to the input
                n_new_axes = 1 + len(max_param_shape) - input_ndim
                new_axes = (1,) * n_new_axes
                new_shape = new_axes + _input.shape
                pivot = model_set_axis_param
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
                new_input = np.rollaxis(_input, model_set_axis_input,
                                        pivot + 1)
        pivots.append(pivot)
        reshaped.append(new_input)

    if model.n_inputs < model.n_outputs:
        pivots.extend([model_set_axis_input] * (model.n_outputs - model.n_inputs))

    return reshaped, (pivots,)


def _prepare_outputs_model_set(model, outputs, format_info, model_set_axis):
    pivots = format_info[0]
    # If model_set_axis = False was passed then use
    # model._model_set_axis to format the output.
    if model_set_axis is None or model_set_axis is False:
        model_set_axis = model.model_set_axis
    outputs = list(outputs)
    for idx, output in enumerate(outputs):
        pivot = pivots[idx]
        if pivot < output.ndim and pivot != model_set_axis:
            outputs[idx] = np.rollaxis(output, pivot,
                                       model_set_axis)
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
            if input_shape[model_set_axis] != n_models:
                try:
                    argname = argnames[idx]
                except IndexError:
                    # the case of model.inputs = ()
                    argname = str(idx)
                raise ValueError(
                    "Input argument {0!r} does not have the correct "
                    "dimensions in model_set_axis={1} for a model set with "
                    "n_models={2}.".format(argname, model_set_axis,
                                           n_models))
        all_shapes.append(input_shape)

    input_shape = check_consistent_shapes(*all_shapes)
    if input_shape is None:
        raise ValueError(
            "All inputs must have identical shapes or must be scalars.")

    return input_shape


def remove_axes_from_shape(shape, axis):
    """
    Given a shape tuple as the first input, construct a new one by  removing
    that particular axis from the shape and all preceeding axes. Negative axis
    numbers are permittted, where the axis is relative to the last axis.
    """
    if len(shape) == 0:
        return shape
    if axis < 0:
        axis = len(shape) + axis
        return shape[:axis] + shape[axis+1:]
    if axis >= len(shape):
        axis = len(shape)-1
    shape = shape[axis+1:]
    return shape


def check_consistent_shapes(*shapes):
    """
    Given shapes as arguments, check to see if all are the same (excluding
    scalars, i.e., shape==(); if all the same, return the common shape; if
    not, return None)
    """
    # remove scalars from the list
    ashapes = [shape for shape in shapes if shape != ()]
    if len(ashapes) == 0:
        return ()
    if len(ashapes) == 1:
        return ashapes[0]
    rshape = ashapes[0]
    for shape in ashapes[1:]:
        if shape != rshape:
            return None
    return rshape


def get_bounding_box(self):
    """
    Return the ``bounding_box`` of a model.

    Raises
    ------
    NotImplementedError
        If ``bounding_box`` is not defined.
    """
    try:
        bbox = self.bounding_box
    except NotImplementedError:
        bbox = None
    return bbox


def generic_call(self, *inputs, **kwargs):
    """ The base ``Model. __call__`` method."""
    inputs, format_info = self.prepare_inputs(*inputs, **kwargs)
    if isinstance(self, CompoundModel):
        # CompoundModels do not normally hold parameters at that level
        parameters = ()
    else:
        parameters = self._param_sets(raw=True, units=True)
    with_bbox = kwargs.pop('with_bounding_box', False)
    fill_value = kwargs.pop('fill_value', np.nan)
    bbox = get_bounding_box(self)
    if with_bbox and bbox is not None:
        input_shape = _validate_input_shapes(
            inputs, self.inputs, self._n_models, self.model_set_axis,
            self.standard_broadcasting)
        vinputs, valid_ind, allout = prepare_bounding_box_inputs(
            self, input_shape, inputs, bbox)
        valid_result_unit = None
        if not allout:
            valid_result = self.evaluate(*chain(vinputs, parameters))
            valid_result_unit = getattr(valid_result, 'unit', None)
            if self.n_outputs == 1:
                valid_result = [valid_result]
            outputs = prepare_bounding_box_outputs(valid_result, valid_ind,
                                                   input_shape, fill_value)
        else:
            outputs = [np.zeros(input_shape) + fill_value for i in range(self.n_outputs)]
        if valid_result_unit is not None:
            outputs = Quantity(outputs, valid_result_unit, copy=False)
    else:
        outputs = self.evaluate(*chain(inputs, parameters))
        if self.n_outputs == 1:
            outputs = (outputs,)
    outputs = self.prepare_outputs(format_info, *outputs, **kwargs)

    outputs = self._process_output_units(inputs, outputs)

    if self.n_outputs == 1:
        return outputs[0]
    return outputs


def prepare_bounding_box_inputs(self, input_shape, inputs, bbox):
    """
    Assign a value of ``np.nan`` to indices outside the bounding box.
    """
    allout = False
    if self.n_inputs > 1:
        # bounding_box is in python order -
        # convert it to the order of the inputs
        bbox = bbox[::-1]
    if self.n_inputs == 1:
        bbox = [bbox]
    # indices where input is outside the bbox
    # have a value of 1 in ``nan_ind``
    nan_ind = np.zeros(input_shape, dtype=bool)
    for ind, inp in enumerate(inputs):
        inp = np.asanyarray(inp)
        outside = np.logical_or(inp < bbox[ind][0], inp > bbox[ind][1])
        if inp.shape:
            nan_ind[outside] = True
        else:
            nan_ind |= outside
            if nan_ind:
                allout = True
    # get an array with indices of valid inputs
    valid_ind = np.atleast_1d(np.logical_not(nan_ind)).nonzero()
    if len(valid_ind[0]) == 0:
        allout = True
    # inputs holds only inputs within the bbox
    args = []
    if not allout:
        for inp in inputs:
            if input_shape:
                args.append(np.array(inp)[valid_ind])
            else:
                args.append(inp)
    return args, valid_ind, allout


def prepare_bounding_box_outputs(valid_result, valid_ind, input_shape, fill_value):
    """
    Populate the output arrays with ``fill_value``.
    """
    result = [np.zeros(input_shape) + fill_value
              for vr in valid_result]
    for ind, r in enumerate(valid_result):
        if not result[ind].shape:
            result[ind] = np.array(r)
        else:
            result[ind][valid_ind] = r
    return result


def _strip_ones(intup):
    return tuple(item for item in intup if item != 1)


def hide_inverse(model):
    """
    This is a convenience function intended to disable automatic generation
    of the inverse in compound models by disabling one of the constituent
    model's inverse. This is to handle cases where user provided inverse
    functions are not compatible within an expression.

    Example:
        compound_model.inverse = hide_inverse(m1) + m2 + m3

    This will insure that the defined inverse itself won't attempt to
    build its own inverse, which would otherwise fail in this example
    (e.g., m = m1 + m2 + m3 happens to raises an exception for this
    reason.)

    Note that this permanently disables it. To prevent that either copy
    the model or restore the inverse later.
    """
    del model.inverse
    return model

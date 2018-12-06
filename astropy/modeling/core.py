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

import abc
import copy
import copyreg
import inspect
import functools
import operator
import types
import warnings

from collections import defaultdict, OrderedDict
from inspect import signature
from itertools import chain, islice

import numpy as np

from astropy.utils import indent, metadata
from astropy.table import Table
from astropy.units import Quantity, UnitsError, dimensionless_unscaled
from astropy.units.utils import quantity_asanyarray
from astropy.utils import (sharedmethod, find_current_module,
                           OrderedDescriptorContainer,
                           check_broadcast, IncompatibleShapeError, isiterable)
from astropy.utils.codegen import make_function_with_signature
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.misc import get_parameters
from .utils import (combine_labels, make_binary_operator_eval,
                    ExpressionTree, AliasDict, get_inputs_and_params,
                    _BoundingBox, _combine_equivalency_dict)
from astropy.nddata.utils import add_array, extract_array
#from . import compound
from .parameters import (Parameter, InputParameterError,
                         param_repr_oneline, _tofloat)

from collections import deque


from ..utils import indent
from .utils import combine_labels, _BoundingBox

__all__ = ['Model', 'FittableModel', 'Fittable1DModel', 'Fittable2DModel',
           'CompoundModel', 'custom_model', 'ModelDefinitionError']

def _model_oper(oper, **kwargs):
    """
    This is an alternate version of compound models intended to use
    much less memory than the default.
    """
    return lambda left, right: CompoundModel(oper, left, right, **kwargs)


class _CompoundModel:
    pass


class _CompoundModelMeta:
    pass


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

    def _opfunc(left, right):
        # Deprecation is for https://github.com/astropy/astropy/issues/8234
        if not (isinstance(left, Model) and isinstance(right, Model)):
            warnings.warn(
                'Composition of model classes will be removed in 4.0 '
                '(but composition of model instances is not affected)',
                AstropyDeprecationWarning)

        # Perform an arithmetic operation on two models.
        return _CompoundModelMeta._from_operator(oper, left, right, **kwargs)

    return _opfunc


class _ModelMeta(InheritDocstrings, abc.ABCMeta):
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
            ('__and__',_model_oper('&')),
            ###('__mod__', _model_oper('%'))
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
            <class 'astropy.modeling.core.SkyRotation'>
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

        mod = find_current_module(2)
        if mod:
            modname = mod.__name__
        else:
            modname = '__main__'

        new_cls = type(name, (cls,), {})
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
            # TODO: Maybe warn in the above case?
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

        if ('__call__' not in members and 'inputs' in members and
                isinstance(members['inputs'], tuple)):

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
            inputs = members['inputs']
            args = ('self',) + inputs
            new_call = make_function_with_signature(
                    __call__, args, [('model_set_axis', None),
                                     ('with_bounding_box', False),
                                     ('fill_value', np.nan),
                                     ('equivalencies', None)])

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
            if all(p.default is not None
                   for p in pdict.values()):
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
    __pow__ =     _model_oper('**')
    __or__ =      _model_oper('|')
    __and__ =     _model_oper('&')
    ###__mod__ =     _model_oper('%')

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
        except Exception:
            # If any of the above formatting fails fall back on the basic repr
            # (this is particularly useful in debugging)
            return parts[0]


class Model(metaclass=_ModelMeta):
    """
    Warning: DOCSTRING IS OUT OF DATE!
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
                    self.__dict__[parname] = newpar
                    #newpar._validator = val._validator
        self._initialize_constraints(kwargs)
        # Remaining keyword args are either parameter values or invalid
        # Parameter values must be passed in as keyword arguments in order to
        # distinguish them
        self._initialize_parameters(args, kwargs)
        self._initialize_slices()
        self._initialize_unit_support()

    def _initialize_unit_support(self):
        """
        Convert self._input_units_strict and
        self.input_units_allow_dimensionless to dictionaries
        mapping input name to a boolena value.
        """
        if isinstance(self._input_units_strict, bool):
            self._input_units_strict = {key: self._input_units_strict for
                                        key in self.__class__.inputs}

        if isinstance(self._input_units_allow_dimensionless, bool):
            self._input_units_allow_dimensionless = {key: self._input_units_allow_dimensionless
                                                     for key in self.__class__.inputs}

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
            return {key: val for key in self.__class__.inputs}
        else:
            return val

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
            return {key: val for key in self.__class__.inputs}
        else:
            return val

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
        else:
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
                    raise UnitsError("The '{0}' parameter should be given as a"
                                     " Quantity because it was originally "
                                     "initialized as a Quantity".format(
                                                            param.name))
                else:
                    param._unit = value.unit
                    param.value = value.value
        else:
            super().__setattr__(attr, value)

    def __call__(self, *inputs, **kwargs):
        """
        Evaluate this model using the given input(s) and the parameter values
        that were specified when the model was instantiated.
        """

        return generic_call(self, *inputs, **kwargs)

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
        ``(min, max)`` tuples or ``[min, max]`` lists.
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

        return self._user_inverse

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
        else:
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
        quantities need to define a _parameter_units_for_data_units method
        that takes the input and output units (as two dictionaries) and
        returns a dictionary giving the target units for each parameter.
        """

        model = self.copy()

        inputs_unit = {inp: getattr(kwargs[inp], 'unit',
                       dimensionless_unscaled)
                       for inp in self.inputs if kwargs[inp] is not None}

        outputs_unit = {out: getattr(kwargs[out], 'unit',
                        dimensionless_unscaled)
                        for out in self.outputs if kwargs[out] is not None}
        parameter_units = self._parameter_units_for_data_units(inputs_unit,
                                                               outputs_unit)

        for name, unit in parameter_units.items():
            parameter = getattr(model, name)
            if parameter.unit is not None:
                parameter.value = parameter.quantity.to(unit).value
                parameter._set_unit(None, force=True)

        if isinstance(model, _CompoundModel):
            model.strip_units_from_tree()

        return model

    def strip_units_from_tree(self):
        for item in self._tree.traverse_inorder():
            if isinstance(item.value, Model):
                for parname in item.value.param_names:
                    par = getattr(item.value, parname)
                    par._set_unit(None, force=True)
                    setattr(item.value, parname, par)

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
        _parameter_units_for_data_units method that takes the input and output
        units (as two dictionaries) and returns a dictionary giving the target
        units for each parameter.
        """

        model = self.copy()

        inputs_unit = {inp: getattr(kwargs[inp], 'unit',
                       dimensionless_unscaled)
                       for inp in self.inputs if kwargs[inp] is not None}

        outputs_unit = {out: getattr(kwargs[out], 'unit',
                        dimensionless_unscaled)
                        for out in self.outputs if kwargs[out] is not None}

        parameter_units = self._parameter_units_for_data_units(inputs_unit,
                                                               outputs_unit)

        # We are adding units to parameters that already have a value, but we
        # don't want to convert the parameter, just add the unit directly,
        # hence the call to _set_unit.
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
        # If the model has a '_parameter_units_for_data_units' method, this
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

        inputs = self._validate_input_units(inputs, equivalencies)

        # The input formatting required for single models versus a multiple
        # model set are different enough that they've been split into separate
        # subroutines
        if n_models == 1:
            return _prepare_inputs_single_model(self, params, inputs,
                                                **kwargs)
        else:
            return _prepare_inputs_model_set(self, params, inputs, n_models,
                                             model_set_axis, **kwargs)

    def _validate_input_units(self, inputs, equivalencies=None):

        inputs = list(inputs)
        name = self.name or self.__class__.__name__
        # Check that the units are correct, if applicable

        if self.input_units is not None:

            # We combine any instance-level input equivalencies with user
            # specified ones at call-time.
            input_units_equivalencies = \
                _combine_equivalency_dict(self.inputs,
                                          equivalencies,
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
                            equivalencies=
                            input_units_equivalencies[input_name]):

                        # If equivalencies have been specified, we need to
                        # convert the input to the input units - this is
                        # because some equivalencies are non-linear, and
                        # we need to be sure that we evaluate the model in
                        # its own frame of reference. If input_units_strict
                        # is set, we also need to convert to the input units.
                        if len(input_units_equivalencies) > 0 or \
                               self.input_units_strict[input_name]:
                            inputs[i] = \
                                inputs[i].to(
                                    input_unit,
                                    equivalencies=input_units_equivalencies[
                                        input_name])

                    else:

                        # We consider the following two cases separately so as
                        # to be able to raise more appropriate/nicer exceptions

                        if input_unit is dimensionless_unscaled:
                            raise UnitsError("Units of input '{0}', {1} ({2}),"
                                             "could not be converted to "
                                             "required dimensionless "
                                             "input".format(self.inputs[i],
                                                            inputs[i].unit,
                                                            inputs[i].unit.
                                                            physical_type))
                        else:
                            raise UnitsError("Units of input '{0}', {1} ({2}),"
                                             " could not be "
                                             "converted to required input"
                                             " units of {3} ({4})".format(
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
                            raise UnitsError("Units of input '{0}', (dimensionless), could not be "
                                             "converted to required input units of "
                                             "{1} ({2})".format(self.inputs[i], input_unit,
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
            return _prepare_outputs_single_model(self, outputs, format_info)
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

    @sharedmethod
    def n_submodels(self):
        """
        Return the number of components in a single model, which is
        obviously 1.
        """
        return 1

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
            self.__init__()
        if isinstance(self, CompoundModel):
            # Need to set parameter attributes
            self._parameters_ = \
                [getattr(existing, param_name) for param_name in param_names]
            for param_name in param_names:
                self.__dict__[param_name] = getattr(existing, param_name)
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
        # self._parameters = existing._parameters

        # self._param_metrics = defaultdict(dict)
        # for param_a, param_b in aliases.items():
        #     # Take the param metrics info for the giving parameters in the
        #     # existing model, and hand them to the appropriate parameters in
        #     # the new model
        #     self._param_metrics[param_a] = existing._param_metrics[param_b]

        for param_a, param_b in aliases.items():
            setattr(self, param_a, getattr(existing, param_b))
        if needs_initialization:
            self.__init__(*dummy_args)

        return self

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

    @property
    def _constraints(self):
        """
        Extract parameter constraints into dictionary of dictionaries
        """
        constraints = {}
        for constraint in self.parameter_constraints:
            tdict = {}
            for param_name in self.param_names:
                param = getattr(self, param_name)
                tdict[param_name] = getattr(param, constraint)
            constraints[constraint] = tdict
        for constraint in self.model_constraints:
            constraints[constraint] = self._mconstraints[constraint]

        return constraints

    @_constraints.setter
    def _constraints(self, constraints):
        self._initialize_constraints(constraints)

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

        supplied_parvalues = []
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
        ## now validate parameters
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
                    raise TypeError(
                        "{0}.__init__() requires a value for parameter "
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
            unit = None
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
        param_names = []
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

        #param_metrics = self._param_metrics
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

        # TODO: Returning an array from this method may be entirely pointless
        # for internal use--perhaps only the external param_sets method should
        # return an array (and just for backwards compat--I would prefer to
        # maybe deprecate that method)

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

    inputs = ('x',)
    outputs = ('y',)
    _separable = True


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

SPECIAL_OPERATORS = {}


def _add_special_operator(sop_name, sop):
    SPECIAL_OPERATORS[sop_name] = sop

"""
This module provides an alternate implementation of compound models that
is lighter weight than the default implementation.

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

- picklingt
  compound models (this implementation walks the tree every time)
- and other things I've overlooked at this moment...

Things that will never be supported:

- Compound models of model classes (as opposed to instances)
"""


class CompoundModel(Model):
    '''
    Lightweight compound model implementation
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
        self._parameters = None
        self._parameters_ = None
        self._param_metrics = None
        self._has_inverse = False  # may be set to True in following code
        if inverse:
            self._user_inverse = inverse
        else:
            self._user_inverse = None
        if op != '%' and len(left) != len(right):
            raise ValueError(
                'Both operands must have equal values for n_models')
        else:
            self._n_models = len(left)
        if op in ['+', '-', '*', '/', '**'] or op in SPECIAL_OPERATORS:
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
                self._inverse = CompoundModel('&',
                                              self.left.inverse,
                                              self.right.inverse,
                                              inverse=self)
        elif op == '|':
            if left.n_outputs != right.n_inputs:
                raise ModelDefinitionError(
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
        else:
            raise ModelDefinitionError('Illegal operator: ', self.op)
        if inverse is not None:
            self._inverse = inverse
            self._has_inverse = True
        self.name = name
        self._fittable = None
        self.fit_deriv = None
        self.col_fit_deriv = None
        if op in ('|', '+', '-'):
            self.linear = left.linear and right.linear
        else:
            self.linear = False
        self.eqcons = False
        self.ineqcons = False

    def __len__(self):
        return self._n_models

    def evaluate(self, *args, **kwargs):
        pass

    @property
    def n_submodels(self):
        if self._leaflist is None:
            self._make_leaflist()
        return len(self._leaflist)

    @property
    def submodel_names(self):
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
            elif op in SPECIAL_OPERATORS:
                return binary_operation(SPECIAL_OPERATORS[op], leftval, rightval)
        else:
            raise ModelDefinitionError('unrecognized operator {op}')

    @property
    def param_names(self):
        if self._param_names is None:
            self.map_parameters()
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
        self.map_parameters()
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
    @property
    def n_inputs(self):
        return self._n_inputs

    @n_inputs.setter
    def n_inputs(self, value):
        self._n_inputs = value

    @property
    def n_outputs(self):
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
            self.map_parameters()
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
        if self.has_inverse():
            if self._user_inverse is not None:
                return self._user_inverse
            else:
                return self._inverse
        else:
            raise NotImplementedError("Inverse function not provided")

    @inverse.setter
    def inverse(self, value):
        self._inverse = value

    @property
    def fittable(self):
        if self._fittable is None:
            if self._leaflist is None:
                self.map_parameters()
            self._fittable = all(m.fittable for m in self._leaflist)
        return self._fittable

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
        if not isinstance(invmodel, Model):
            raise ValueError("Attempt to assign non model to inverse")
        self._has_inverse = True
        self._user_inverse = invmodel

    __add__ =     _model_oper('+')
    __sub__ =     _model_oper('-')
    __mul__ =     _model_oper('*')
    __truediv__ = _model_oper('/')
    __pow__ =     _model_oper('**')
    __or__ =      _model_oper('|')
    __and__ =     _model_oper('&')
    ###__mod__ =     _model_oper('%')

    def map_parameters(self):
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
        self._param_names = []
        for lindex, leaf in enumerate(self._leaflist):
            for param_name in leaf.param_names:
                param = getattr(leaf, param_name)
                new_param_name = "{}_{}".format(param_name, lindex)
                self.__dict__[new_param_name] = param
                self._parameters_[new_param_name] = param
                self._param_names.append(new_param_name)
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
            param_metrics[name]['size'] = param_size
            total_size += param_size
        self._parameters = np.empty(total_size, dtype=np.float64)

    @property
    def _constraints(self):
        return self._constraints_compound

    @_constraints.setter
    def _constraints(self, value):
        self._constraints_compound = value

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
        ('inputs', tuple(x.name for x in inputs)),
        ('outputs', output_names),
        ('evaluate', staticmethod(func))]
    )

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
                outputs[idx] = output.item()
            else:
                outputs[idx] = output.reshape(broadcast_shape)

    return tuple(outputs)


def _prepare_inputs_model_set(model, params, inputs, n_models, model_set_axis_input,
                              **kwargs):
    reshaped = []
    pivots = []

    model_set_axis_param = model.model_set_axis # needed to reshape param
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


def generic_call(self, *inputs, **kwargs):
    inputs, format_info = self.prepare_inputs(*inputs, **kwargs)
    # Check whether any of the inputs are quantities
    inputs_are_quantity = any([isinstance(i, Quantity) for i in inputs])
    if isinstance(self, CompoundModel):
        # CompoundModels do not normally hold parameters at that level
        parameters = ()
    else:
        parameters = self._param_sets(raw=True, units=True)
    with_bbox = kwargs.pop('with_bounding_box', False)
    fill_value = kwargs.pop('fill_value', np.nan)
    bbox = None
    if with_bbox:
        try:
            bbox = self.bounding_box
        except NotImplementedError:
            bbox = None
        if self.n_inputs > 1 and bbox is not None:
            # bounding_box is in python order -
            # convert it to the order of the inputs
            bbox = bbox[::-1]
        if bbox is None:
            outputs = self.evaluate(*chain(inputs, parameters))
        else:
            if self.n_inputs == 1:
                bbox = [bbox]
            # indices where input is outside the bbox
            # have a value of 1 in ``nan_ind``
            nan_ind = np.zeros(inputs[0].shape, dtype=bool)
            for ind, inp in enumerate(inputs):
                # Pass an ``out`` array so that ``axis_ind``
                # is array for scalars as well.
                axis_ind = np.zeros(inp.shape, dtype=bool)
                axis_ind = np.logical_or(inp < bbox[ind][0],
                                         inp > bbox[ind][1], out=axis_ind)
                nan_ind[axis_ind] = 1
            # get an array with indices of valid inputs
            valid_ind = np.logical_not(nan_ind).nonzero()
            # inputs holds only inputs within the bbox
            args = []
            for input in inputs:
                if not input.shape:
                    # shape is ()
                    if nan_ind:
                        outputs = [fill_value for a in args]
                    else:
                        args.append(input)
                else:
                    args.append(input[valid_ind])
            valid_result = self.evaluate(*chain(args, parameters))
            if self.n_outputs == 1:
                valid_result = [valid_result]
            # combine the valid results with the ``fill_value`` values
            # outside the bbox
            result = [np.zeros(inputs[0].shape) + fill_value
                      for i in range(len(valid_result))]
            for ind, r in enumerate(valid_result):
                if not result[ind].shape:
                    # shape is ()
                    result[ind] = r
                else:
                    result[ind][valid_ind] = r
            # format output
            if self.n_outputs == 1:
                outputs = np.asarray(result[0])
            else:
                outputs = [np.asarray(r) for r in result]
    else:
        outputs = self.evaluate(*chain(inputs, parameters))
    if self.n_outputs == 1:
        outputs = (outputs,)
    outputs = self.prepare_outputs(format_info, *outputs, **kwargs)

    outputs = self._process_output_units(inputs, outputs)

    if self.n_outputs == 1:
        return outputs[0]
    else:
        return outputs

def _strip_ones(intup):
    return tuple(item for item in intup if item !=1)

# def ismodel(obj):
#     """
#     Returns True if the object is an instance of Model or CompoundModel.

#     This should be used instead of isinstance since CompoundModel does
#     not inherit from Model for efficiency reasons.
#     """
#     return isinstance(obj, Model) or isinstance(obj, CompoundModel)


copyreg.pickle(_ModelMeta, _ModelMeta.__reduce__)

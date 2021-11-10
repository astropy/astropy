# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['quantity_input', 'dequantity_input']

import inspect
from numbers import Number
from collections.abc import Sequence
from functools import wraps

import numpy as np

from . import _typing as T
from .core import (Unit, UnitBase, UnitsError,
                   add_enabled_equivalencies, dimensionless_unscaled)
from .function.core import FunctionUnitBase
from .physical import PhysicalType, get_physical_type
from .quantity import Quantity
from .structured import StructuredUnit
from .unitspec import UnitSpec, NullUnitSpec, UnitSpecBase


def _is_unitlike(target, allow_structured=True):
    """Check if target is `~astropy.units.Unit`-like.

    Parameters
    ----------
    target : Any
    allow_structured : bool
        Whether to count a `~astropy.units.StructuredUnit` as a
        `~astropy.units.Unit`, allowing for a distinction between a list of
        units and a structured unit.

    Returns
    -------
    bool
        True if target is `~astropy.units.Unit`-like, False otherwise.
    """
    # check if unit-like
    try:
        unit = Unit(target)
    except (TypeError, ValueError):
        return False

    # distinguish between list of units and a structured unit
    ulike = True
    if (not allow_structured
        and not isinstance(target, StructuredUnit)
        and isinstance(unit, StructuredUnit)):
        ulike = False

    return ulike


def _is_physicaltypelike(target):
    """Check if target is `~astropy.units.PhysicalType`-like.

    Parameters
    ----------
    target : Any
        recognized types are:
        - PhysicalType
        - str
        - `~astropy.units.Unit`
        - `~astropy.units.Quantity`

    Returns
    -------
    bool
        True if target is `~astropy.units.PhysicalType`-like, False otherwise.
    """
    try:
        get_physical_type(target)
    except TypeError:
        return False
    return True


def _parse_annotation(target):

    if target in (None, inspect._empty):
        return None

    # check if unit-like
    try:
        unit = Unit(target)
    except (TypeError, ValueError):
        try:
            ptype = get_physical_type(target)
        except (TypeError, ValueError, KeyError):  # KeyError for Enum
            if isinstance(target, str):
                raise ValueError(f"invalid unit or physical type {target!r}.") from None
        else:
            return ptype
    else:
        return unit

    # could be a type hint
    origin = T.get_origin(target)
    if origin is T.Union:
        return [_parse_annotation(t) for t in T.get_args(target)]
    elif origin is not T.Annotated:  # can't be Quantity[]
        return False

    # parse type hint
    cls, *annotations = T.get_args(target)
    if not issubclass(cls, Quantity) or not annotations:
        return False

    # get unit from type hint
    unit, *rest = annotations
    if not isinstance(unit, (UnitBase, PhysicalType)):
        return None

    return unit


def _parse_target(target, error=False):
    """Parse unit / physical type target.

    Parameters
    ----------
    target : Any
        recognized types are:
        - `~astropy.units.UnitSpecBase`
        - list (not tuple) of recognized types
        - unit-like
        - PhysicalType-like

    error : bool
        Whether to raise a `TypeError` if 'target' is not one of its recognized
        types.

    Returns
    -------
    `~astropy.units.UnitSpecBase`, |Unit|, |PhysicalType|, None, or list thereof

    Raises
    ------
    TypeError
        If 'error' is True and 'target' is not one of its recognized types.
    """
    if target is None:
        spec = NullUnitSpec()
    elif isinstance(target, UnitSpecBase):
        spec = target
    # determine it's a unit

    elif _is_unitlike(target, allow_structured=False):
        spec = Unit(target)
    # determine if PhysicalType-like
    elif _is_physicaltypelike(target):
        spec = get_physical_type(target)
    elif isinstance(target, Sequence):
        spec = [_parse_target(t, error=error) for t in target]
    # none of the recognized types. either error or pass through.
    elif not error:
        spec = None
    else:
        raise TypeError(f"{target} is neither a Unit or PhysicalType"
                        " and cannot be parsed into either.") from None

    return spec


class QuantityInput:
    """Quantity-aware functions, with behavior set by unit specifications.

    Parameters
    ----------
    equivalencies : list
    strict_dimensionless : bool, optional
    action : str, optional
        The `~astropy.units.UnitSpec` 'action' for each `~inspect.Parameter`
        in ``sig``.
    return_action : str or None, optional

    """

    @classmethod
    def as_decorator(cls, func=None, **kwargs):
        r"""
        A decorator for validating the units of arguments to functions.

        Unit specifications can be provided as keyword arguments to the
        decorator, or by using function annotation syntax. Arguments to the
        decorator take precedence over any function annotations present.

        A `~astropy.units.UnitsError` will be raised if the unit attribute of
        the argument is not equivalent to the unit specified to the decorator or
        in the annotation. If the argument has no unit attribute, i.e. it is not
        a Quantity object, a `ValueError` will be raised unless the argument is
        an annotation. This is to allow non Quantity annotations to pass
        through.

        Where an equivalency is specified in the decorator, the function will be
        executed with that equivalency in force.

        Notes
        -----
        The original function is accessible by the attributed ``__wrapped__``.
        See :func:`functools.wraps` for details.

        Examples
        --------

        .. code-block:: python

            import astropy.units as u
            @u.quantity_input(myangle=u.arcsec)
            def myfunction(myangle):
                return myangle**2


        .. code-block:: python

            import astropy.units as u
            @u.quantity_input
            def myfunction(myangle: u.arcsec):
                return myangle**2

        Also you can specify a return value annotation, which will
        cause the function to always return a `~astropy.units.Quantity` in that
        unit.

        .. code-block:: python

            import astropy.units as u
            @u.quantity_input
            def myfunction(myangle: u.arcsec) -> u.deg**2:
                return myangle**2

        Using equivalencies::

            import astropy.units as u
            @u.quantity_input(myenergy=u.eV, equivalencies=u.mass_energy())
            def myfunction(myenergy):
                return myenergy**2

        """
        self = cls(**kwargs)
        # allow for pie-syntax without parentheses -- ``@as_decorator``.
        if func is not None and not kwargs:
            return self(func)
        return self

    def __init__(self, bound=True, strict_dimensionless=False,
                 equivalencies=(), action="validate", return_action="to unit",
                 parameters={}, **kwargs):
        self.bound = bound
        self.strict_dimensionless = strict_dimensionless

        self.equivalencies = equivalencies
        self.decorator_kwargs = {**parameters, ** kwargs}

        # UnitSpec action
        self.action = action
        self.return_action = return_action

    def __call__(self, wrapped_function):

        wsig = inspect.signature(wrapped_function)
        wrapped_function.__signature__ = wsig
        uspecs = self._make_uspecs_from_signature(wsig)

        # Define the wrapper function
        # inside this function is evaluated every function call
        @wraps(wrapped_function)
        def wrapper(*args, **kwargs):

            # Bind the arguments to our new function to the signature of the original.
            ba = wrapper.__wrapped__.__signature__.bind(*args, **kwargs)
            ba.apply_defaults()  # ensures (kw)args exist, even if blank
            arguments = ba.arguments  # store reference. updates in-place.

            with add_enabled_equivalencies(wrapper.__units_specs__["equivalencies"]):

                # this updates with external changes to __units_specs__
                for name, uspec in wrapper.__units_specs__["uspecs"].items():
                    # normal argument
                    if name in arguments:
                        arguments[name] = uspec(arguments[name])
                    # variable keyword arguments (apply uspec over all)
                    elif name.startswith("**"):
                        kwargs = arguments[name[2:]]  # updates in-place.
                        for k, v in kwargs.items():
                            kwargs[k] = uspec(v)
                    # variable arguments (apply uspec over all)
                    elif name.startswith("*"):
                        args = [uspec(v) for v in arguments[name[1:]]]
                        arguments[name[1:]] = tuple(args)
                    # callable's return annotation (guaranteed to exist)
                    # TODO? speed up by making optional
                    # TODO! allow to handle multiple outputs
                    elif name == "return":
                        return_spec = uspec  # store for later
                    else:
                        raise Exception(f"no arg for uspec {name}")

                return_ = wrapped_function(*ba.args, **ba.kwargs)
                return return_spec(return_)

        # ---------

        # store all necessary information on the wrapper
        wrapper.__units_specs__ = dict(uspecs=uspecs,
                                       equivalencies=self.equivalencies)

        # method to update wrapper if e.g. parameter defaults are changed
        def update_wrapper(equivalencies=None, strict_dimensionless=None):
            """Update (in-place) {0!r}'s QuantityInput wrapper.

            Parameters
            ----------
            equivalencies : list or None, optional
                New equivalencies for wrapper. Only set if not None.
            strict_dimensionless : bool or None, optional
                New flag for wrapper. Only set if not None.
            """
            self._update_wrapper(wrapper, equivalencies=equivalencies,
                                 strict_dimensionless=strict_dimensionless)

        update_wrapper.__doc__ = update_wrapper.__doc__.format(wrapped_function)
        wrapper.update_wrapper = update_wrapper

        return wrapper

    def _make_uspecs_from_signature(self, sig):
        """Get `~astropy.units.Unit` specifications from `~inspect.Signature`.

        Parameters
        ----------
        sig : `~inspect.Signature`
            Function signature from which to create `~astropy.units.UnitSpec`
            (or subclass of kind ``uspec_cls``)

        Returns
        -------
        uspecs : dict[str, `~astropy.units.UnitSpec`]
            Key is parameter name. If parameter is the variable positional
            argument the name is prepended with a '*'. Similarly, the variable
            keyword parameter's name is prepended with a '**'.
            The value is the corresponding `~astropy.units.UnitSpec`.
        """
        uspecs = {}

        # 1) parse (kw)args
        for param in sig.parameters.values():
            name, spec = self._make_uspec_from_parameter(param, action=self.action)
            if spec is not None:  # spec is None if no target
                uspecs[name] = spec

        # 2) parse return annotation
        param = inspect.Parameter("return", 0, default=inspect.Parameter.empty,
                                  annotation=sig.return_annotation)
        _, spec = self._make_uspec_from_parameter(param, action=self.return_action)
        uspecs["return"] = spec if spec is not None else NullUnitSpec()
        # always wrap in UnitSpec for code simplicity in wrapper

        return uspecs

    def _make_uspec_from_parameter(self, param, action):
        """Make `~astropy.units.UnitSpec` from `~inspect.Parameter`

        Parameters
        ----------
        param : `~inspect.Parameter`
        action : str
            The `~astropy.units.UnitSpec` 'action'

        Returns
        -------
        name : str
        uspec : `~astropy.units.UnitSpec`
        """
        # name is parameter name, unless variable positional or keyword,
        # in which case it is prepended by '*' or '**', respectively
        name = (param.name if param.kind not in (2, 4)
                else ("*", "**")[param.kind // 2 - 1] + param.name)
        dec_key = param.name if name != "return" else "return_"

        # annotation and decorator arg
        dnote = self.decorator_kwargs.get(dec_key, ...)  # ... = pass-thru
        anote = param.annotation

        # parse target:
        # 1) skip if set target=False in decorator
        # 2) decorator overrides annotations, unless ... (Ellipsis)
        # 3) parse the annotation, skipping if empty
        if dnote is False:  # skip if set parameter to False in decorator
            name = target = None
        elif dnote is not ...:  # unless pass-thru, prioritize decorator arguments
            # decorator arguments must be correct, so error if not.
            target = _parse_target(dnote, error=True)
        elif anote is inspect.Parameter.empty:  # no dec arg nor annotation
            name = target = None
        else:  # annotation (becomes `None` if not recognized)
            target = _parse_annotation(anote)

        # target -> Unitspec, composite thereof, or None (no target)
        if isinstance(target, UnitSpecBase) or target is None:
            # don't check the uspec if already a UnitSpec. This allows for
            # the 'action' to be ignored.  # TODO? option to enforce 'action'?
            uspec = target
        elif isinstance(target, list):  # multiple targets -> CompoundUnitSpec
            uspec = np.sum([UnitSpec(t, action=action, bound=self.bound,
                                     strict_dimensionless=self.strict_dimensionless)
                            for t in target])
            # note that UnitSpec does not convert, so if any 't' in target
            # was a UnitSpec, 'action' did not apply
        else:
            uspec = UnitSpec(target, action=action, bound=self.bound,
                             strict_dimensionless=self.strict_dimensionless)

        # If the default value is None, allow the argument value to be None
        # even if there is a target unit.
        if param.default is None and target is not None:
            uspec = uspec + NullUnitSpec()

        return name, uspec

    def _update_wrapper(self, wrapper, equivalencies=None, bound=None, strict_dimensionless=None):
        """Update (in-place) a QuantityInput wrapped functions' wrapper.

        Parameters
        ----------
        wrapper : callable

        equivalencies : list or None, optional
            New equivalencies for wrapper. Only set if not None.
        bound : bool or None, optional
            New flag for wrapper. Only set if not None.
        strict_dimensionless : bool or None, optional
            New flag for wrapper. Only set if not None.
        """
        # get and store updated signature
        wrapped_sig = inspect.signature(wrapper.__wrapped__)
        wrapper.__wrapped__.__signature__ = wrapped_sig
        # process signature for uspecs
        wrapper.__units_specs__["uspecs"] = self._make_uspecs_from_signature(wrapped_sig)
        # optionally update equivalencies and strict_dimensionless
        if equivalencies is not None:
            wrapper.__units_specs__["equivalencies"] = equivalencies

        if bound is not None:
            for uspec in wrapper.__units_specs__["uspecs"].values():
                uspec.bound = bound

        if strict_dimensionless is not None:
            # set on all relevant UnitSpecs
            for uspec in wrapper.__units_specs__["uspecs"].values():
                if hasattr(uspec, "strict_dimensionless"):
                    uspec.strict_dimensionless = strict_dimensionless


quantity_input = QuantityInput.as_decorator


# -------------------------------------------------------------------

class DeQuantityInput(QuantityInput):
    """QuantityInput but strips units beforehand and reinstates after.

    Parameters
    ----------
    func
    strict_dimensionless
    equivalencies
    parameters
    **kwargs
    """

    def __init__(self, func=None, strict_dimensionless=False,
                 equivalencies=(), parameters={}, **kwargs):
        super().__init__(func=func, strict_dimensionless=strict_dimensionless,
                         equivalencies=equivalencies, parameters=parameters,
                         action="to value", return_action="from value",
                         **kwargs)


dequantity_input = DeQuantityInput.as_decorator

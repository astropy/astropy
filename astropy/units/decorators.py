# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['quantity_input']

from numbers import Number
from collections.abc import Sequence
import inspect

import numpy as np

from astropy.utils.decorators import wraps
from .core import (Unit, UnitBase, UnitsError, add_enabled_equivalencies,
                   dimensionless_unscaled)
from .physical import _unit_physical_mapping


def _get_allowed_units(targets):
    """
    From a list of target units (either as strings or unit objects) and physical
    types, return a list of Unit objects.
    """

    allowed_units = []
    for target in targets:

        try:  # unit passed in as a string
            target_unit = Unit(target)

        except ValueError:

            try:  # See if the function writer specified a physical type
                physical_type_id = _unit_physical_mapping[target]

            except KeyError:  # Function argument target is invalid
                raise ValueError(f"Invalid unit or physical type '{target}'.")

            # get unit directly from physical type id
            target_unit = Unit._from_physical_type_id(physical_type_id)

        allowed_units.append(target_unit)

    return allowed_units


def _validate_arg_value(param_name, func_name, arg, targets, equivalencies,
                        strict_dimensionless=False):
    """
    Validates the object passed in to the wrapped function, ``arg``, with target
    unit or physical type, ``target``.
    """

    if len(targets) == 0:
        return

    allowed_units = _get_allowed_units(targets)

    # If dimensionless is an allowed unit and the argument is unit-less,
    #   allow numbers or numpy arrays with numeric dtypes
    if (dimensionless_unscaled in allowed_units and not strict_dimensionless
            and not hasattr(arg, "unit")):
        if isinstance(arg, Number):
            return

        elif (isinstance(arg, np.ndarray)
              and np.issubdtype(arg.dtype, np.number)):
            return

    for allowed_unit in allowed_units:
        try:
            is_equivalent = arg.unit.is_equivalent(allowed_unit,
                                                   equivalencies=equivalencies)

            if is_equivalent:
                break

        except AttributeError:  # Either there is no .unit or no .is_equivalent
            if hasattr(arg, "unit"):
                error_msg = ("a 'unit' attribute without an 'is_equivalent' "
                             "method")
            else:
                error_msg = "no 'unit' attribute"

            raise TypeError(f"Argument '{param_name}' to function '{func_name}'"
                            f" has {error_msg}. You should pass in an astropy "
                            "Quantity instead.")

    else:
        error_msg = (f"Argument '{param_name}' to function '{func_name}' must "
                     "be in units convertible to")
        if len(targets) > 1:
            targ_names = ", ".join([f"'{str(targ)}'" for targ in targets])
            raise UnitsError(f"{error_msg} one of: {targ_names}.")
        else:
            raise UnitsError(f"{error_msg} '{str(targets[0])}'.")


class QuantityInput:

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

        The checking of arguments inside variable arguments to a function is not
        supported (i.e. \*arg or \**kwargs).

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
        if func is not None and not kwargs:
            return self(func)
        else:
            return self

    def __init__(self, func=None, strict_dimensionless=False, **kwargs):
        self.equivalencies = kwargs.pop('equivalencies', [])
        self.decorator_kwargs = kwargs
        self.strict_dimensionless = strict_dimensionless

    def __call__(self, wrapped_function):

        # Extract the function signature for the function we are wrapping.
        wrapped_signature = inspect.signature(wrapped_function)

        # Define a new function to return in place of the wrapped one
        @wraps(wrapped_function)
        def wrapper(*func_args, **func_kwargs):
            # Bind the arguments to our new function to the signature of the original.
            bound_args = wrapped_signature.bind(*func_args, **func_kwargs)

            # Iterate through the parameters of the original signature
            for param in wrapped_signature.parameters.values():
                # We do not support variable arguments (*args, **kwargs)
                if param.kind in (inspect.Parameter.VAR_KEYWORD,
                                  inspect.Parameter.VAR_POSITIONAL):
                    continue

                # Catch the (never triggered) case where bind relied on a default value.
                if (param.name not in bound_args.arguments
                        and param.default is not param.empty):
                    bound_args.arguments[param.name] = param.default

                # Get the value of this parameter (argument to new function)
                arg = bound_args.arguments[param.name]

                # Get target unit or physical type, either from decorator kwargs
                #   or annotations
                if param.name in self.decorator_kwargs:
                    targets = self.decorator_kwargs[param.name]
                    is_annotation = False
                else:
                    targets = param.annotation
                    is_annotation = True

                # If the targets is empty, then no target units or physical
                #   types were specified so we can continue to the next arg
                if targets is inspect.Parameter.empty:
                    continue

                # If the argument value is None, and the default value is None,
                #   pass through the None even if there is a target unit
                if arg is None and param.default is None:
                    continue

                # Here, we check whether multiple target unit/physical type's
                #   were specified in the decorator/annotation, or whether a
                #   single string (unit or physical type) or a Unit object was
                #   specified
                if (isinstance(targets, str)
                        or not isinstance(targets, Sequence)):
                    valid_targets = [targets]

                # Check for None in the supplied list of allowed units and, if
                #   present and the passed value is also None, ignore.
                elif None in targets:
                    if arg is None:
                        continue
                    else:
                        valid_targets = [t for t in targets if t is not None]

                else:
                    valid_targets = targets

                # If we're dealing with an annotation, skip all the targets that
                #    are not strings or subclasses of Unit. This is to allow
                #    non unit related annotations to pass through
                if is_annotation:
                    valid_targets = [t for t in valid_targets
                                     if isinstance(t, (str, UnitBase))]

                # Now we loop over the allowed units/physical types and validate
                #   the value of the argument:
                _validate_arg_value(param.name, wrapped_function.__name__,
                                    arg, valid_targets, self.equivalencies,
                                    self.strict_dimensionless)

            # Call the original function with any equivalencies in force.
            with add_enabled_equivalencies(self.equivalencies):
                return_ = wrapped_function(*func_args, **func_kwargs)

            valid_empty = (inspect.Signature.empty, None)
            if wrapped_signature.return_annotation not in valid_empty:
                return return_.to(wrapped_signature.return_annotation)
            else:
                return return_

        return wrapper


quantity_input = QuantityInput.as_decorator

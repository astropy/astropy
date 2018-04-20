# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['quantity_input']

from ..utils.decorators import wraps
from ..utils.compat import funcsigs
from ..utils.misc import isiterable

from .core import Unit, UnitsError, add_enabled_equivalencies
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
                raise ValueError("Invalid unit or physical type '{0}'."
                                 .format(target))

            # get unit directly from physical type id
            target_unit = Unit._from_physical_type_id(physical_type_id)

        allowed_units.append(target_unit)

    return allowed_units


def _validate_arg_value(param_name, func_name, arg, targets, equivalencies):
    """
    Validates the object passed in to the wrapped function, ``arg``, with target
    unit or physical type, ``target``.
    """

    allowed_units = _get_allowed_units(targets)

    for allowed_unit in allowed_units:
        try:
            is_equivalent = arg.unit.is_equivalent(allowed_unit,
                                                   equivalencies=equivalencies)

            if is_equivalent:
                break

        except AttributeError:  # Either there is no .unit or no .is_equivalent
            if hasattr(arg, "unit"):
                error_msg = "a 'unit' attribute without an 'is_equivalent' method"
            else:
                error_msg = "no 'unit' attribute"

            raise TypeError("Argument '{0}' to function '{1}' has {2}. "
                  "You may want to pass in an astropy Quantity instead."
                     .format(param_name, func_name, error_msg))

    else:
        if len(targets) > 1:
            raise UnitsError("Argument '{0}' to function '{1}' must be in units"
                             " convertible to one of: {2}."
                             .format(param_name, func_name,
                                     [str(targ) for targ in targets]))
        else:
            raise UnitsError("Argument '{0}' to function '{1}' must be in units"
                             " convertible to '{2}'."
                             .format(param_name, func_name,
                                     str(targets[0])))


class QuantityInput(object):

    @classmethod
    def as_decorator(cls, func=None, **kwargs):
        r"""
        A decorator for validating the units of arguments to functions.

        Unit specifications can be provided as keyword arguments to the decorator,
        or by using Python 3's function annotation syntax. Arguments to the decorator
        take precedence over any function annotations present.

        A `~astropy.units.UnitsError` will be raised if the unit attribute of
        the argument is not equivalent to the unit specified to the decorator
        or in the annotation.
        If the argument has no unit attribute, i.e. it is not a Quantity object, a
        `ValueError` will be raised.

        Where an equivalency is specified in the decorator, the function will be
        executed with that equivalency in force.

        Notes
        -----

        The checking of arguments inside variable arguments to a function is not
        supported (i.e. \*arg or \**kwargs).

        Examples
        --------

        Python 2 and 3::

            import astropy.units as u
            @u.quantity_input(myangle=u.arcsec)
            def myfunction(myangle):
                return myangle**2

        Python 3 only:

        .. code-block:: python3

            import astropy.units as u
            @u.quantity_input
            def myfunction(myangle: u.arcsec):
                return myangle**2

        Also in Python 3 you can specify a return value annotation, which will
        cause the function to always return a `~astropy.units.Quantity` in that
        unit.

        .. code-block:: python3

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

    def __init__(self, func=None, **kwargs):
        self.equivalencies = kwargs.pop('equivalencies', [])
        self.decorator_kwargs = kwargs

    def __call__(self, wrapped_function):

        # Extract the function signature for the function we are wrapping.
        wrapped_signature = funcsigs.signature(wrapped_function)

        # Define a new function to return in place of the wrapped one
        @wraps(wrapped_function)
        def wrapper(*func_args, **func_kwargs):
            # Bind the arguments to our new function to the signature of the original.
            bound_args = wrapped_signature.bind(*func_args, **func_kwargs)

            # Iterate through the parameters of the original signature
            for param in wrapped_signature.parameters.values():
                # We do not support variable arguments (*args, **kwargs)
                if param.kind in (funcsigs.Parameter.VAR_KEYWORD,
                                  funcsigs.Parameter.VAR_POSITIONAL):
                    continue

                # Catch the (never triggered) case where bind relied on a default value.
                if param.name not in bound_args.arguments and param.default is not param.empty:
                    bound_args.arguments[param.name] = param.default

                # Get the value of this parameter (argument to new function)
                arg = bound_args.arguments[param.name]

                # Get target unit or physical type, either from decorator kwargs
                #   or annotations
                if param.name in self.decorator_kwargs:
                    targets = self.decorator_kwargs[param.name]
                else:
                    targets = param.annotation

                # If the targets is empty, then no target units or physical
                #   types were specified so we can continue to the next arg
                if targets is funcsigs.Parameter.empty:
                    continue

                # If the argument value is None, and the default value is None,
                #   pass through the None even if there is a target unit
                if arg is None and param.default is None:
                    continue

                # Here, we check whether multiple target unit/physical type's
                #   were specified in the decorator/annotation, or whether a
                #   single string (unit or physical type) or a Unit object was
                #   specified
                if isinstance(targets, str) or not isiterable(targets):
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

                # Now we loop over the allowed units/physical types and validate
                #   the value of the argument:
                _validate_arg_value(param.name, wrapped_function.__name__,
                                    arg, valid_targets, self.equivalencies)

            # Call the original function with any equivalencies in force.
            with add_enabled_equivalencies(self.equivalencies):
                return_ = wrapped_function(*func_args, **func_kwargs)
            if wrapped_signature.return_annotation not in (funcsigs.Signature.empty, None):
                return return_.to(wrapped_signature.return_annotation)
            else:
                return return_

        return wrapper


quantity_input = QuantityInput.as_decorator

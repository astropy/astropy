# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['quantity_input']

from ..utils.decorators import wraps
from ..utils.compat import funcsigs

from .core import (Unit, UnitsError, add_enabled_equivalencies,
                   get_current_unit_registry)
from .physical import _unit_physical_mapping

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

                # Get target unit, either from decorator kwargs or annotations
                if param.name in self.decorator_kwargs:
                    target_unit = self.decorator_kwargs[param.name]
                else:
                    target_unit = param.annotation

                # If the target unit is empty, then no unit was specified so we
                # move past it
                if target_unit is not funcsigs.Parameter.empty:

                    if isinstance(target_unit, str):

                        try: # unit passed in as a string
                            target_unit = Unit(target_unit)
                            str_target_unit = target_unit.to_string()
                        except ValueError:
                            # user specified a physical type instead of a unit
                            try:
                                physical_type_id = _unit_physical_mapping[target_unit]
                                str_target_unit = target_unit
                            except KeyError:
                                raise ValueError("Invalid unit of physical type '{0}'."
                                                 .format(target_unit))

                            ureg = get_current_unit_registry()
                            target_units = ureg._by_physical_type[physical_type_id]
                            target_unit = next(iter(target_units)) # get first valid unit from set

                    else:
                        str_target_unit = target_unit.to_string()

                    try:
                        equivalent = arg.unit.is_equivalent(target_unit,
                                                  equivalencies=self.equivalencies)

                        if not equivalent:
                            raise UnitsError("Argument '{0}' to function '{1}'"
                                             " must be in units convertible to"
                                             " '{2}'.".format(param.name,
                                                              wrapped_function.__name__,
                                                              str_target_unit))

                    # Either there is no .unit or no .is_equivalent
                    except AttributeError:
                        if hasattr(arg, "unit"):
                            error_msg = "a 'unit' attribute without an 'is_equivalent' method"
                        else:
                            error_msg = "no 'unit' attribute"
                        raise TypeError("Argument '{0}' to function '{1}' has {2}. "
                              "You may want to pass in an astropy Quantity instead."
                                 .format(param.name, wrapped_function.__name__, error_msg))

            # Call the original function with any equivalencies in force.
            with add_enabled_equivalencies(self.equivalencies):
                return_ = wrapped_function(*func_args, **func_kwargs)
            if wrapped_signature.return_annotation is not funcsigs.Signature.empty:
                return return_.to(wrapped_signature.return_annotation)
            else:
                return return_

        return wrapper

quantity_input = QuantityInput.as_decorator

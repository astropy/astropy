# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['quantity_input']

from ..utils.decorators import wraps
from ..utils.compat import funcsigs

from .core import UnitsError, add_enabled_equivalencies

class QuantityInput(object):

    @classmethod
    def as_decorator(cls, func=None, **kwargs):
        """
        A decorator for validating the units of arguments to functions.

        Unit specifications can be provided as keyword arguments to the decorator,
        or by using Python 3's function annotation syntax. Arguments to the decorator
        take precidence over any function annotations present.

        A `~astropy.units.UnitsError` will be raised if the unit attribute of
        the argument is not equivalent to the unit specified to the decorator
        or in the annotation.
        If the argument has no unit attribute, i.e. it is not a Quantity object, a
        `~exceptions.ValueError` will be raised.

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

        Python 3 only::

            import astropy.units as u
            @u.quantity_input
            def myfunction(myangle: u.arcsec):
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
                    try:
                        equivalent = arg.unit.is_equivalent(target_unit,
                                                  equivalencies=self.equivalencies)

                        if not equivalent:
                            raise UnitsError("Argument '{0}' to function '{1}'"
                                             " must be in units convertable to"
                                             " '{2}'.".format(param.name,
                                                     wrapped_function.__name__,
                                                     target_unit.to_string()))

                    # Either there is no .unit or no .is_equivalent
                    except AttributeError:
                        if hasattr(arg, "unit"):
                            error_msg = "a 'unit' attribute without an 'is_equivalent' method"
                        else:
                            error_msg = "no 'unit' attribute"
                        raise TypeError("Argument '{0}' to function has '{1}' {2}. "
                              "You may want to pass in an astropy Quantity instead."
                                 .format(param.name, wrapped_function.__name__, error_msg))

            # Call the original function with any equivalencies in force.
            with add_enabled_equivalencies(self.equivalencies):
                return wrapped_function(*func_args, **func_kwargs)

        return wrapper

quantity_input = QuantityInput.as_decorator

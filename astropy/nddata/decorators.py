# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import warnings

from ..utils import wraps
from ..utils.exceptions import AstropyUserWarning
from ..utils.compat.funcsigs import signature

from .nddata import NDData

__all__ = ['support_nddata']


def support_nddata(_func=None, accepts=NDData, repack=False, returns=None):
    """
    Decorator to split NDData properties into function arguments.

    This is a decorator to allow functions to take NDData objects as their
    first arguments and split up the properties into kwargs as required by the
    function. For example, if you consider the following function::

        def downsample(data, wcs=None):
            # downsample data and optionally WCS here
            pass

    This function takes a Numpy array for the data, and some WCS information
    with the ``data`` keyword argument. However, you might have an NDData
    instance that has the ``wcs`` property set and you would like to be able to
    call the function with ``downsample(my_nddata)`` and have the WCS
    information, if present, automatically be passed to the ``wcs`` keyword
    argument.

    This decorator can be used to make this possible::

        @support_nddata
        def downsample(data, wcs=None):
            # downsample data and optionally WCS here
            pass

    This function can now either be called as before, specifying the data and
    WCS separately, or an NDData instance can be passed to the ``data``
    argument.

    The restrictions on functions to use this function are:

    * The first positional argument should be ``data`` and take a Numpy array.

    * The following arguments can optionally be specified in the function
      signature, but if they are specified they should be keyword arguments:
      ``uncertainty``, ``mask``, ``meta``, ``unit``, and ``wcs``. If
      you are making use of this decorator, you should be prepared for these
      keyword arguments to be set to the properties of the NDData object (if
      present).

    The behavior of the decorator is to check through the NDData properties and
    if they are set, it checks if the function accepts them as keyword
    arguments. If an NDData property is set but cannot be passed to a keyword
    argument, a warning is emitted to tell the user that the NDData property in
    question will not be used by the function (to ensure that they know when
    e.g. uncertainties cannot be used).

    If the user passes an NDData object *and* explicitly sets a keyword
    argument that is one of the valid NDData properties, a warning is emitted
    to inform the user that the explicitly specified value will take priority.
    """

    if returns is not None and not repack:
        raise ValueError('returns should only be set if repack=True')

    if returns is None and repack:
        raise ValueError('returns should be set if repack=True')

    def support_nddata_decorator(func):

        # Find out args and kwargs
        sig = signature(func)
        func_args = []
        func_kwargs = []
        for param in sig.parameters.values():
            if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
                raise ValueError("func may not have *args or **kwargs")
            elif param.default == param.empty:
                func_args.append(param.name)
            else:
                func_kwargs.append(param.name)

        # First argument should be data
        if len(func_args) == 0 or func_args[0] != 'data':
            raise ValueError("Can only wrap functions whose first positional argument is `data`")

        supported_properties = ['uncertainty', 'mask', 'meta', 'unit', 'wcs']

        @wraps(func)
        def wrapper(data, *args, **kwargs):

            unpack = isinstance(data, accepts)
            input_data = data

            if not unpack and isinstance(data, NDData):
                raise TypeError("Only NDData sub-classes that inherit from {0}"
                                " can be used by this function".format(accepts.__name__))

            # If data is an NDData instance, we can try and find properties that
            # can be passed as kwargs.
            if unpack:

                ignored = []

                # We loop over a list of pre-defined properties
                for prop in supported_properties:

                    # We only need to do something if the property exists on the
                    # NDData object
                    if hasattr(data, prop):
                        value = getattr(data, prop)
                        if (prop == 'meta' and len(value) > 0) or (prop != 'meta' and value is not None):
                            if prop in func_kwargs:
                                if prop in kwargs and kwargs[prop] is not None:
                                    warnings.warn("Property {0} has been passed explicitly and as an "
                                                  "NDData property, using explicitly specified value".format(prop),
                                                  AstropyUserWarning)
                                else:
                                    kwargs[prop] = value
                            else:
                                ignored.append(prop)

                if ignored:
                    warnings.warn("The following attributes were set on the data object, "
                                  "but will be ignored by the function: " + ", ".join(ignored),
                                  AstropyUserWarning)

                # Finally, replace data by the data itself
                data = data.data

            result = func(data, *args, **kwargs)

            if unpack:

                if repack:
                    if len(returns) > 1 and len(returns) != len(result):
                        raise ValueError("Function did not return the expected number of arguments")
                    elif len(returns) == 1:
                        result = [result]
                    return input_data.__class__(**dict(zip(returns, result)))
                else:
                    return result

            else:

                return result

        return wrapper

    # If _func is set, this means that the decorator was used without
    # parameters so we have to return the result of the support_nddata_decorator
    # decorator rather than the decorator itself
    if _func is not None:
        return support_nddata_decorator(_func)
    else:
        return support_nddata_decorator

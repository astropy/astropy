import inspect
import warnings

from ..utils.exceptions import AstropyUserWarning

from .nddata import NDData

__all__ = ['expand_nddata_args']

# TODO: preserve signature and docstring
# TODO: write tests
# TODO: write docstring and documentation


def expand_nddata_args(func):
    """
    This is a decorator to allow functions to take NDData objects as their
    first arguments and split up the properties into kwargs as required by the
    function.
    """

    # Find out args and kwargs
    wrapped_argspec = inspect.getargspec(func)

    # Find out the args and kwargs
    if wrapped_argspec.defaults:
        func_args = wrapped_argspec.args[:-len(wrapped_argspec.defaults)]
        func_kwargs = wrapped_argspec.args[len(func_args):]
    else:
        func_args = wrapped_argspec.args
        func_kwargs = []

    # For now, let's be strict to simplify things
    if func_args != ['data']:
        raise ValueError("Can only wrap functions that have a single positional data argument")

    supported_properties = ['wcs', 'unit', 'uncertainty', 'mask', 'meta']

    def wrapper(data, **kwargs):

        # If data is an NDData instance, we can try and find properties that
        # can be passed as kwargs.
        if isinstance(data, NDData):

            ignored = []

            # We loop over a list of pre-defined properties
            for prop in supported_properties:

                # We only need to do something if the property exists on the
                # NDData object
                if hasattr(data, prop):
                    value = getattr(data, prop)
                    if value is not None:
                        if prop in func_kwargs:
                            if prop in kwargs and kwargs[prop] is not None:
                                warnings.warn("Property {0} has been passed explicitly and as an "
                                              "NDData property, using explicitly specified value".format(prop),
                                              AstropyUserWarning)
                            kwargs[prop] = value
                        else:
                            ignored.append(prop)

            if ignored:
                warnings.warn("The following attributes were set on the data object, "
                              "but will be ignored by the function: " + ", ".join(ignored),
                              AstropyUserWarning)

            # Finally, replace data by the data itself
            data = data.data

        return func(data, **kwargs)

    return wrapper

# Licensed under a 3-clause BSD style license - see LICENSE.rst


import warnings
from copy import deepcopy
from functools import wraps
from inspect import signature
from itertools import islice

from astropy.utils.exceptions import AstropyUserWarning

from .nddata import NDData

__all__ = ["support_nddata"]


# All supported properties are optional except "data" which is mandatory!
SUPPORTED_PROPERTIES = ["data", "uncertainty", "mask", "meta", "unit", "wcs", "flags"]


def support_nddata(
    _func=None,
    accepts=NDData,
    repack=False,
    returns=None,
    keeps=None,
    **attribute_argument_mapping,
):
    """Decorator to wrap functions that could accept an NDData instance with
    its properties passed as function arguments.

    Parameters
    ----------
    _func : callable, None, optional
        The function to decorate or ``None`` if used as factory. The first
        positional argument should be ``data`` and take a numpy array. It is
        possible to overwrite the name, see ``attribute_argument_mapping``
        argument.
        Default is ``None``.

    accepts : class, optional
        The class or subclass of ``NDData`` that should be unpacked before
        calling the function.
        Default is ``NDData``

    repack : bool, optional
        Should be ``True`` if the return should be converted to the input
        class again after the wrapped function call.
        Default is ``False``.

        .. note::
           Must be ``True`` if either one of ``returns`` or ``keeps``
           is specified.

    returns : iterable, None, optional
        An iterable containing strings which returned value should be set
        on the class. For example if a function returns data and mask, this
        should be ``['data', 'mask']``. If ``None`` assume the function only
        returns one argument: ``'data'``.
        Default is ``None``.

        .. note::
           Must be ``None`` if ``repack=False``.

    keeps : iterable. None, optional
        An iterable containing strings that indicate which values should be
        copied from the original input to the returned class. If ``None``
        assume that no attributes are copied.
        Default is ``None``.

        .. note::
           Must be ``None`` if ``repack=False``.

    attribute_argument_mapping :
        Keyword parameters that optionally indicate which function argument
        should be interpreted as which attribute on the input. By default
        it assumes the function takes a ``data`` argument as first argument,
        but if the first argument is called ``input`` one should pass
        ``support_nddata(..., data='input')`` to the function.

    Returns
    -------
    decorator_factory or decorated_function : callable
        If ``_func=None`` this returns a decorator, otherwise it returns the
        decorated ``_func``.

    Notes
    -----
    If properties of ``NDData`` are set but have no corresponding function
    argument a Warning is shown.

    If a property is set of the ``NDData`` are set and an explicit argument is
    given, the explicitly given argument is used and a Warning is shown.

    The supported properties are:

    - ``mask``
    - ``unit``
    - ``wcs``
    - ``meta``
    - ``uncertainty``
    - ``flags``

    Examples
    --------
    This function takes a Numpy array for the data, and some WCS information
    with the ``wcs`` keyword argument::

        def downsample(data, wcs=None):
            # downsample data and optionally WCS here
            pass

    However, you might have an NDData instance that has the ``wcs`` property
    set and you would like to be able to call the function with
    ``downsample(my_nddata)`` and have the WCS information, if present,
    automatically be passed to the ``wcs`` keyword argument.

    This decorator can be used to make this possible::

        @support_nddata
        def downsample(data, wcs=None):
            # downsample data and optionally WCS here
            pass

    This function can now either be called as before, specifying the data and
    WCS separately, or an NDData instance can be passed to the ``data``
    argument.
    """
    if (returns is not None or keeps is not None) and not repack:
        raise ValueError("returns or keeps should only be set if repack=True.")
    elif returns is None and repack:
        raise ValueError("returns should be set if repack=True.")
    else:
        # Use empty lists for returns and keeps so we don't need to check
        # if any of those is None later on.
        if returns is None:
            returns = []
        if keeps is None:
            keeps = []

    # Short version to avoid the long variable name later.
    attr_arg_map = attribute_argument_mapping
    if any(keep in returns for keep in keeps):
        raise ValueError("cannot specify the same attribute in `returns` and `keeps`.")
    all_returns = returns + keeps

    def support_nddata_decorator(func):
        func_sig = signature(func)
        func_parameters = func_sig.parameters
        for param in func_parameters.values():
            if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
                raise ValueError("func may not have *args or **kwargs.")

        data_arg_name = attr_arg_map.get("data", "data")
        if data_arg_name not in func_parameters.keys():
            raise ValueError(
                f"Can only wrap a function with a {data_arg_name} argument."
            )

        @wraps(func)
        def wrapper(*args, **kwargs):
            bound_args = func_sig.bind_partial(*args, **kwargs)
            data = bound_args.arguments[data_arg_name]
            unpack = isinstance(data, accepts)
            input_data = data
            ignored = []
            if not unpack and isinstance(data, NDData):
                raise TypeError(
                    f"Only NDData sub-classes that inherit from {accepts.__name__}"
                    " can be used by this function"
                )

            # If data is an NDData instance, we can try and find properties
            # that can be passed as kwargs.
            if unpack:
                # We loop over a list of pre-defined properties
                for prop in islice(SUPPORTED_PROPERTIES, 1, None):
                    # We only need to do something if the property exists on
                    # the NDData object
                    try:
                        value = getattr(data, prop)
                    except AttributeError:
                        continue
                    # Skip if the property exists but is None or empty.
                    if prop == "meta" and not value:
                        continue
                    elif value is None:
                        continue
                    # Warn if the property is set but not used by the function.
                    propmatch = attr_arg_map.get(prop, prop)
                    if propmatch not in func_parameters:
                        ignored.append(prop)
                        continue

                    # Check if the property was explicitly given and issue a
                    # Warning if it is.
                    if propmatch in bound_args.arguments:
                        # If it's in the func_args it's trivial but if it was
                        # in the func_kwargs we need to compare it to the
                        # default.
                        # Comparison to the default is done by comparing their
                        # identity, this works because defaults in function
                        # signatures are only created once and always reference
                        # the same item.
                        # FIXME: Python interns some values, for example the
                        # integers from -5 to 255 (any maybe some other types
                        # as well). In that case the default is
                        # indistinguishable from an explicitly passed kwarg
                        # and it won't notice that and use the attribute of the
                        # NDData.
                        warn_conflict = False
                        try:
                            if (
                                bound_args.arguments[propmatch]
                                is not func_parameters[propmatch].default
                            ):
                                warn_conflict = True
                        except ValueError as exc:
                            if (
                                "The truth value of an array with more than one element "
                                "is ambiguous." in str(exc)
                            ):
                                warn_conflict = False
                            else:
                                raise

                        if warn_conflict:
                            warnings.warn(
                                "Property {} has been passed explicitly and "
                                "as an NDData property{}, using explicitly "
                                "specified value"
                                "".format(
                                    propmatch, "" if prop == propmatch else " " + prop
                                ),
                                AstropyUserWarning,
                            )
                            continue
                    # Otherwise use the property as input for the function.
                    bound_args.arguments[propmatch] = value
                # Finally, replace data by the data attribute
                bound_args.arguments[data_arg_name] = data.data

                if ignored:
                    warnings.warn(
                        "The following attributes were set on the "
                        "data object, but will be ignored by the "
                        "function: " + ", ".join(ignored),
                        AstropyUserWarning,
                    )

            # Need to apply defaults so that any positional only arguments are
            # passed as positional properly
            bound_args.apply_defaults()
            result = func(*bound_args.args, **bound_args.kwargs)

            if unpack and repack:
                # If there are multiple required returned arguments make sure
                # the result is a tuple (because we don't want to unpack
                # numpy arrays or compare their length, never!) and has the
                # same length.
                if len(returns) > 1:
                    if not isinstance(result, tuple) or len(returns) != len(result):
                        raise ValueError(
                            "Function did not return the expected number of arguments."
                        )
                elif len(returns) == 1:
                    result = [result]
                if keeps is not None:
                    for keep in keeps:
                        result.append(deepcopy(getattr(input_data, keep)))
                resultdata = result[all_returns.index("data")]
                resultkwargs = {
                    ret: res for ret, res in zip(all_returns, result) if ret != "data"
                }
                return input_data.__class__(resultdata, **resultkwargs)
            else:
                return result

        return wrapper

    # If _func is set, this means that the decorator was used without
    # parameters so we have to return the result of the
    # support_nddata_decorator decorator rather than the decorator itself
    if _func is not None:
        return support_nddata_decorator(_func)
    else:
        return support_nddata_decorator

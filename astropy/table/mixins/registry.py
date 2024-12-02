# This module handles the definition of mixin 'handlers' which are functions
# that given an arbitrary object (e.g. a dask array) will return an object that
# can be used as a mixin column. This is useful because it means that users can
# then add objects to tables that are not formally mixin columns and where
# adding an info attribute is beyond our control.

__all__ = ["MixinRegistryError", "get_mixin_handler", "register_mixin_handler"]

# The internal dictionary of handlers maps fully qualified names of classes
# to a function that can take an object and return a mixin-compatible object.
_handlers = {}


class MixinRegistryError(Exception):
    pass


def register_mixin_handler(fully_qualified_name, handler, force=False):
    """
    Register a mixin column 'handler'.

    A mixin column handler is a function that given an arbitrary Python object,
    will return an object with the .info attribute that can then be used as a
    mixin column (this can be e.g. a copy of the object with a new attribute,
    a subclass instance, or a wrapper class - this is left up to the handler).

    The handler will be used on classes that have an exactly matching fully
    qualified name.

    Parameters
    ----------
    fully_qualified_name : str
        The fully qualified name of the class that the handler can operate on,
        such as e.g. ``dask.array.core.Array``.
    handler : func
        The handler function.
    force : bool, optional
        Whether to overwrite any previous handler if there is already one for
        the same fully qualified name.
    """
    if fully_qualified_name not in _handlers or force:
        _handlers[fully_qualified_name] = handler
    else:
        raise MixinRegistryError(
            f"Handler for class {fully_qualified_name} is already defined"
        )


def get_mixin_handler(obj):
    """
    Given an arbitrary object, return the matching mixin handler (if any).

    Parameters
    ----------
    obj : object or str
        The object to find a mixin handler for, or a fully qualified name.

    Returns
    -------
    handler : None or func
        Then matching handler, if found, or `None`
    """
    if isinstance(obj, str):
        return _handlers.get(obj)
    else:
        return _handlers.get(
            obj.__class__.__module__ + "." + obj.__class__.__name__, None
        )


# Add built-in handlers to registry. Note that any third-party package imports
# required by the handlers should go inside the handler function to delay
# the imports until they are actually needed.


def dask_handler(arr):
    from astropy.table.mixins.dask import as_dask_column

    return as_dask_column(arr)


register_mixin_handler("dask.array.core.Array", dask_handler)

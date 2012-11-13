# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains astropy's tools to manage optional dependencies.
"""

def requires_optional_dependencies(*dependencymodnames):
    """
    A decorator to signify that a function

    Parameters
    ----------
    dependencymodnames : str
        A comma-seperated list of

    Returns
    -------
    deco
        A function to be called on another function signifying it

    Examples
    --------

    The primary use case is as a decorator for a function::

        @requires_optional_dependencies('scipy, matplotlib')
        def foo():
            import scipy
            from matplotlib import pyplot as plt

            ...

    It can also be used on methods::

        class AstroFoo(object):

            @requires_optional_dependencies('scipy', 'sympy')
            def __init__(self):

                ...

    """
    from functools import wraps

    if len(dependencymodnames) == 1:
        if callable(dependencymodnames[0]):
            raise TypeError('dependcy list must be provided to requires_optional_dependencies')
        dependencymodnames = [d.strip() for d in dependencymodnames[0].split(',')]

    # This is the decorator that will be applied to the function.
    def real_decorator(fcn):
        # This is the hook that you pass through on the way into
        # the decorated function.
        @wraps(fcn)
        def fn_hook(*args, **kwargs):

            # If we have not successfully imported the optional
            # packages yet, import them.
            if fn_hook.optional_imports_needed:
                missing_packages = []
                glo = globals()
                for x in fn_hook.optional_imports_needed:
                    # If we don't have a global with that name,
                    # try to import it.
                    if not x in glo:
                        try:
                            glo[x] = __import__(x)
                        except ImportError:
                            missing_packages.append(x)

                    # If there is an error with the import, raise
                    # the exception here -- we never get to the
                    # decorated function, AND we do not remember
                    # importing anything.  That means we will try
                    # again if we get called again.
                    if len(missing_packages) > 0:
                        msg = 'Optional packages missing: '
                        raise ImportError(msg + (' '.join(missing_packages)))

                    # Remember that we succeeeded in the imports.
                    fn_hook.optional_imports_unchecked = False

                # finally call through to the real function
                return fcn(*args, **kwargs)

        fn_hook.optional_imports_needed = tuple(dependencymodnames)
        fn_hook.optional_imports_unchecked = True

        return fn_hook

    return real_decorator


def find_all_optional_dependencies(pkgornm=None):
    """ Given a root package name or package, this function walks
    through all the subpackages and modules and searches for optional
    dependencies specified by way of the `requires_optional_dependencies`
    decorator.

    .. note::
        This will import all of the package and subpackage, but
        will *not* fail on ImportErrors by design - instead it will
        silently skip those packages.


    Parameters
    ----------
    pkgornm : module, str, or None
        The package (as a module object) or name of the package to
        search. If None, it determines the package to search as the root
        package of the function where this function is called.
    """
    from pkgutil import get_loader, walk_packages
    from types import ModuleType
    from inspect import isclass

    from .misc import find_current_module

    if pkgornm is None:
        pkgornm = find_current_module(1).__name__.split('.')[0]

    if isinstance(pkgornm, basestring):
        package = get_loader(pkgornm).load_module(pkgornm)
    elif isinstance(pkgornm, ModuleType) and '__init__' in pkgornm.__file__:
        package = pkgornm
    else:
        msg = 'find_all_optional_dependencies was not given a package/package name'
        raise TypeError(msg)

    def do_check(modorcls):
        """
        Actually check for the opdeps - this is a function to allow recursion
        """
        s = set()
        for k, v in modorcls.__dict__.iteritems():
            if isclass(v):  # need to check methods of classes
                s = s.union(do_check(v))
            elif hasattr(v, 'optional_imports_needed'):
                s = s.union(v.optional_imports_needed)
        return s

    opdeps = set()
    for imper, nm, ispkg in walk_packages(package.__path__, package.__name__ + '.'):
        imper.find_module(nm)
        try:
            mod = __import__(nm)
            opdeps = opdeps.union(do_check(mod))
        except ImportError:
            pass

    return opdeps

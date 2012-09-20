# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains astropy's tools to manage optional dependencies.
"""


def requires_optional_dependencies(*module_name_list):
    """ A docstring
    """

    if len(module_name_list) == 1:
        module_name_list.split(',')

    # This is the decorator that will be applied to the function.
    def real_decorator(fcn):
        # This is the hook that you pass through on the way into
        # the decorated function.
        def fn_hook(*args, **kwargs):

            # If we have not successfully imported the optional
            # packages yet, import them.
            if fn_hook.option_imports_needed:
                missing_packages = []
                glo = globals()
                for x in module_name_list:
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
                    fn_hook.option_imports_needed = False

                # finally call through to the real function
                return fcn(*args, **kwargs)

        fn_hook.option_imports_needed = True
        return fn_hook

    return real_decorator



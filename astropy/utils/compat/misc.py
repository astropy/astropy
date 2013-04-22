# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Simple utility functions and bug fixes for compatibility with all supported
versions of Python.  This module should generally not be used directly, as
everything in `__all__` will be imported into `astropy.utils.compat` and can
be accessed from there.

Includes the following fixes:

* The `inspect.getmodule` function does not always work in Python 3.1 and 3.2.
  This package includes a function `inspect_getmodule` that will simply be an
  alias to `inspect.getmodule` if the stdlib version is correct, but for
  versions of python with the bug, it uses an internal patched version.

"""

from __future__ import absolute_import

from functools import wraps
from sys import version_info

__all__ = ['inspect_getmodule', 'invalidate_caches', 'override__dir__']


def _patched_getmodule(object, _filename=None):
    """Return the module an object was defined in, or None if not found.

    This replicates the functionality of the stdlib `inspect.getmodule`
    function but includes a fix for a bug present in Python 3.1 and 3.2.
    """
    #these imports mock up what would otherwise have been in inspect
    import sys
    import os
    from inspect import modulesbyfile, _filesbymodname, getabsfile, ismodule

    if ismodule(object):
        return object
    if hasattr(object, '__module__'):
        return sys.modules.get(object.__module__)
    # Try the filename to modulename cache
    if _filename is not None and _filename in modulesbyfile:
        return sys.modules.get(modulesbyfile[_filename])
    # Try the cache again with the absolute file name
    try:
        file = getabsfile(object, _filename)
    except TypeError:
        return None
    if file in modulesbyfile:
        return sys.modules.get(modulesbyfile[file])
    # Update the filename to module name cache and check yet again
    # Copy sys.modules in order to cope with changes while iterating
    # This is where the fix is made - the adding of the "list" call:
    for modname, module in list(sys.modules.items()):
        if ismodule(module) and hasattr(module, '__file__'):
            f = module.__file__
            if f == _filesbymodname.get(modname, None):
                # Have already mapped this module, so skip it
                continue
            _filesbymodname[modname] = f
            f = getabsfile(module)
            # Always map to the name the module knows itself by
            modulesbyfile[f] = modulesbyfile[
                os.path.realpath(f)] = module.__name__
    if file in modulesbyfile:
        return sys.modules.get(modulesbyfile[file])
    # Check the main module
    main = sys.modules['__main__']
    if not hasattr(object, '__name__'):
        return None
    if hasattr(main, object.__name__):
        mainobject = getattr(main, object.__name__)
        if mainobject is object:
            return main
    # Check builtins
    builtin = sys.modules['builtins']
    if hasattr(builtin, object.__name__):
        builtinobject = getattr(builtin, object.__name__)
        if builtinobject is object:
            return builtin

inspect_getmodule = None
"""
An alias to `inspect.getmodule`, or a patched version that replicates the
functionality with a bugfix for Python 3.1 and 3.2.
"""

#This assigns the stdlib inspect.getmodule to the variable name
#`inspect_getmodule` if it's not buggy, and uses the matched version if it is.
if version_info[0] < 3 or version_info[1] > 2:
    #in 2.x everythig is fine, as well as >=3.3
    from inspect import getmodule as inspect_getmodule
else:
    inspect_getmodule = _patched_getmodule


# Python 3.3's importlib caches filesystem reads for faster imports in the
# general case. But sometimes it's necessary to manually invalidate those
# caches so that the import system can pick up new generated files.  See
# https://github.com/astropy/astropy/issues/820
if version_info[:2] >= (3, 3):
    from importlib import invalidate_caches
else:
    invalidate_caches = lambda: None


def override__dir__(f):
    """
    When overriding a __dir__ method on an object, you often want to
    include the "standard" members on the object as well.  This
    decorator takes care of that automatically, and all the wrapped
    function needs to do is return a list of the "special" members
    that wouldn't be found by the normal Python means.

    Example
    -------

    @override__dir__
    def __dir__(self):
        return ['special_method1', 'special_method2']
    """
    if version_info[:2] < (3, 3):
        # There was no straightforward way to do this until Python 3.3, so
        # we have this complex monstrosity
        @wraps(f)
        def override__dir__wrapper(self):
            members = set()
            for cls in self.__class__.mro():
                members.update(dir(cls))
            members.update(self.__dict__.iterkeys())
            members.update(f(self))
            return sorted(members)
    else:
        @wraps(f)
        def override__dir__wrapper(self):
            members = set(super(self.__class__, self).__dir__())
            members.update(f(self))
            return sorted(members)

    return override__dir__wrapper

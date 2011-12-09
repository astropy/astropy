# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains smallish general-purpose utilities that don't have a
clear other module to live in.

This module should not be used directly, as everything in `__all__`` is
imported into `astropy.utils`
"""

__all__ = ['find_current_module']


def find_current_module(depth=1, finddiff=False):
    """ Determines the module/package this function is called from.

    Parameters
    ----------
    depth : int
        Specifies how far back to go in the call stack (0-indexed, so that
        passing in 0 gives back `astropy.utils.misc`).
    finddiff : bool
        If True, once the module at `depth` is determined, a search will be
        performed up the call stack until a *different* module is found
        from the one at `depth`.

    Returns
    -------
    mod : module or None
        The module object or None if the package cannot be found. The name of
        the module is available as the ``__name__`` attribute of the returned
        object (if it isn't None).

    Examples
    --------
    The examples below assume that there are two modules in a package named
    `pkg`. ``mod1.py``:

        def find1():
            from astropy.utils import find_current_module
            print find_current_module(1).__name__
        def find2():
            from astropy.utils import find_current_module
            cmod = find_current_module(2)
            if cmod is None:
                print 'None'
            else:
                print cmod.__name__
        def find_diff():
            from astropy.utils import find_current_module
            print find_current_module(0,True).__name__

    ``mod2.py``:

        def find():
            from .mod1 import find2
            find2()

    With these modules in place, the following occurs:

        >>> from pkg import mod1, mod2
        >>> from astropy.utils import find_current_module
        >>> mod1.find1()
        pkg.mod1
        >>> mod1.find2()
        None
        >>> mod2.find()
        pkg.mod2
        >>> find_current_module(0)
        <module 'astropy.utils.misc' from 'astropy/utils/misc.py'>
        >>> mod1.find_diff()
        pkg.mod1

    """
    from inspect import currentframe

    # using a patched version of getmodule because the py 3.1 and 3.2 stdlib
    # is broken if the list of modules changes during import
    from .compat import inspect_getmodule

    frm = currentframe()
    for i in range(depth):
        frm = frm.f_back
        if frm is None:
            return None

    if finddiff:
        currmod = inspect_getmodule(frm)
        while frm:
            frmb = frm.f_back
            modb = inspect_getmodule(frmb)
            if modb is not currmod:
                return modb
            frm = frmb
    else:
        return inspect_getmodule(frm)

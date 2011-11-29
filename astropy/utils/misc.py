# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains smallish general-purpose utilities that don't have a 
clear other module to live in.

This module should not be used directly, as everything in `__all__`` is 
imported into `astropy.utils`
"""

__all__ = ['find_current_module']

def find_current_module(depth=1):
    """ Determines the module/package this function is called from.  
    
    Parameters
    ----------
    depth : int
        Specifies how far back to go in the call stack. e.g. 0 returns the
        `astropy.utils.general` module, 1 returns the module of this function's
        caller, 2 retreives the caller's caller, etc.
        
    Returns
    -------
    mod : module or None
        The module object or None if the package cannot be found. The name of the
        module is available as the ``__name__`` attribute of the returned object
        (if it isn't None).
    """
    from inspect import currentframe,getmodule
    from sys import version_info
    
    # using a patched version of getmodule because the py 3.1 and 3.2 stdlib
    # is broken if the list of modules changes during import
    from .compat import inspect_getmodule
    
    frm = currentframe()
    for i in range(depth):
        frm = frm.f_back
        if frm is None:
            return None
    return inspect_getmodule(frm)
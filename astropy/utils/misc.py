# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains smallish general-purpose utilities that don't have a
clear other module to live in.

This module should not be used directly, as everything in `__all__`` is
imported into `astropy.utils`
"""
from __future__ import absolute_import

__all__ = ['find_current_module', 'fnpickle', 'fnunpickle']


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
    `pkg`. ``mod1.py``::

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

    ``mod2.py``::

        def find():
            from .mod1 import find2
            find2()

    With these modules in place, the following occurs::

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


def fnunpickle(fileorname, number=0, usecPickle=True):
    """ Unpickle pickled objects from a specified file and return the contents.

    Parameters
    ----------
    fileorname : str or `file`-like
        The file from which to unpickle objects.
    number : int
        If 0, a single object will be returned (the first in the file). If >0,
        this specifies the number of objects to be unpickled, and a list will be
        returned with exactly that many objects. If <0, all objects in the file
        will be unpickled and returned as a list.
    usecPickle : bool
        If True, the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster). This only applies for python 2.x.

    Raises
    ------
    EOFError
        If `number` is >0 and there are fewer than `number` objects in the
        pickled file.

    Returns
    -------
    contents : obj or list
        If `number` is 0, this is a individual object - the first one unpickled
        from the file. Otherwise, it is a list of objects unpickled from the
        file.

    """
    import sys

    if usecPickle and sys.version_info[0] < 3:
        import cPickle as pickle
    else:
        import pickle

    if isinstance(fileorname, basestring):
        f = open(fileorname, 'r')
        close = True
    else:
        f = fileorname
        close = False

    try:
        if number > 0: #get that number
            res = []
            for i in range(number):
                res.append(pickle.load(f))
        elif number < 0: #get all objects
            res = []
            eof = False
            while not eof:
                try:
                    res.append(pickle.load(f))
                except EOFError:
                    eof = True
        else:  # number==0
            res = pickle.load(f)
    finally:
        if close:
            f.close()

    return res


def fnpickle(object, fileorname, usecPickle=True, protocol=None, append=False):
    """Pickle an object to a specified file.

    Parameters
    ----------
    object
        The python object to pickle.
    fileorname : str or `file`-like
        The file into which the `object` should be pickled.
    usecPickle : bool
        If True, the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster). This only applies for python 2.x.
    protocol : int or None
        Pickle protocol to use - see the :mod:`pickle` module for details on
        these options. If None, the most recent protocol will be used.
    append : bool
        If True, the object is appended to the end of the file, otherwise the
        file will be overwritten (if a file object is given instead of a
        file name, this has no effect).

    """
    import sys

    if usecPickle and sys.version_info.major < 3:
        import cPickle as pickle
    else:
        import pickle

    if protocol is None:
        protocol = pickle.HIGHEST_PROTOCOL

    if isinstance(fileorname, basestring):
        f = open(fileorname, 'a' if append else 'w')
        close = True
    else:
        f = fileorname
        close = False

    try:
        pickle.dump(object, f, protocol=protocol)
    finally:
        if close:
            f.close()

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A "grab bag" of smallish general-purpose utilities that don't have a
clear other module/pakcage to live in.
"""

from __future__ import absolute_import

import collections
import functools
import sys
import textwrap
import warnings

__all__ = ['find_current_module', 'fnpickle', 'fnunpickle', 'isiterable',
           'deprecated', 'lazyproperty']


def find_current_module(depth=1, finddiff=False):
    """ Determines the module/package from which this function is called.

    This function has two modes, determined by the `finddiff` option. it
    will either simply go the requested number of frames up the call
    stack (if `finddiff` is False), or it will go up the call stack until
    it reaches a module that is *not* in a specified set.

    Parameters
    ----------
    depth : int
        Specifies how far back to go in the call stack (0-indexed, so that
        passing in 0 gives back `astropy.utils.misc`).
    finddiff : bool or list
        If False, the returned `mod` will just be `depth` frames up from
        the current frame. Otherwise, the function will start at a frame
        `depth` up from current, and continue up the call stack to the
        first module that is *different* from those in the provided list.
        In this case, `finddiff` can be a list of modules or modules
        names. Alternatively, it can be True, which will use the module
        `depth` call stack frames up as the module the returned module
        most be different from.

    Returns
    -------
    mod : module or None
        The module object or None if the package cannot be found. The name of
        the module is available as the ``__name__`` attribute of the returned
        object (if it isn't None).

    Raises
    ------
    ValueError
        If `finddiff` is a list with an invalid entry.

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
    from inspect import currentframe, ismodule

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
        if finddiff is True:
            diffmods = [currmod]
        else:
            diffmods = []
            for fd in finddiff:
                if ismodule(fd):
                    diffmods.append(fd)
                elif isinstance(fd, basestring):
                    diffmods.append(__import__(fd))
                elif fd is True:
                    diffmods.append(currmod)
                else:
                    raise ValueError('invalid entry in finddiff')

        while frm:
            frmb = frm.f_back
            modb = inspect_getmodule(frmb)
            if modb not in diffmods:
                return modb
            frm = frmb
    else:
        return inspect_getmodule(frm)


def find_mod_objs(modname, onlylocals=False):
    """ Returns all the public attributes of a module referenced by name.

    .. note::
        The returned list *not* include subpackages or modules of
        `modname`,nor does it include private attributes (those that
        beginwith '_' or are not in `__all__`).

    Parameters
    ----------
    modname : str
        The name of the module to search.
    onlylocals : bool
        If True, only attributes that are either members of `modname` OR one of
        its modules or subpackages will be included.

    Returns
    -------
    localnames : list of str
        A list of the names of the attributes as they are named in the
        module `modname` .
    fqnames : list of str
        A list of the full qualified names of the attributes (e.g.,
        ``astropy.utils.misc.find_mod_objs``). For attributes that are
        simple variables, this is based on the local name, but for
        functions or classes it can be different if they are actually
        defined elsewhere and just referenced in `modname`.
    objs : list of objects
        A list of the actual attributes themselves (in the same order as
        the other arguments)

    """
    from inspect import ismodule

    __import__(modname)
    mod = sys.modules[modname]

    if hasattr(mod, '__all__'):
        pkgitems = [(k, getattr(mod, k)) for k in mod.__all__]
    else:
        pkgitems = [(k, getattr(mod, k)) for k in dir(mod) if k[0] != '_']

    #filter out modules and pull the names and objs out
    localnames = [k for k, v in pkgitems if not ismodule(v)]
    objs = [v for k, v in pkgitems if not ismodule(v)]

    #fully qualified names can be determined from the object's module
    fqnames = []
    for obj, lnm in zip(objs, localnames):
        if hasattr(obj, '__module__') and hasattr(obj, '__name__'):
            fqnames.append(obj.__module__ + '.' + obj.__name__)
        else:
            fqnames.append(modname + '.' + lnm)

    if onlylocals:
        valids = [fqn.startswith(modname) for fqn in fqnames]
        localnames = [e for i, e in enumerate(localnames) if valids[i]]
        fqnames = [e for i, e in enumerate(fqnames) if valids[i]]
        objs = [e for i, e in enumerate(objs) if valids[i]]

    return localnames, fqnames, objs


def fnunpickle(fileorname, number=0, usecPickle=True):
    """ Unpickle pickled objects from a specified file and return the contents.

    Parameters
    ----------
    fileorname : str or `file`-like
        The file name or file from which to unpickle objects. If a file object,
        it should have been opened in binary mode.
    number : int
        If 0, a single object will be returned (the first in the file). If >0,
        this specifies the number of objects to be unpickled, and a list will
        be returned with exactly that many objects. If <0, all objects in the
        file will be unpickled and returned as a list.
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

    if usecPickle and sys.version_info[0] < 3:  # pragma: py2
        import cPickle as pickle
    else:
        import pickle

    if isinstance(fileorname, basestring):
        f = open(fileorname, 'rb')
        close = True
    else:
        f = fileorname
        close = False

    try:
        if number > 0:  # get that number
            res = []
            for i in range(number):
                res.append(pickle.load(f))
        elif number < 0:  # get all objects
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
        The filename or file into which the `object` should be pickled. If a
        file object, it should have been opened in binary mode.
    usecPickle : bool
        If True (default), the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster). This only applies for python 2.x.
    protocol : int or None
        Pickle protocol to use - see the :mod:`pickle` module for details on
        these options. If None, the most recent protocol will be used.
    append : bool
        If True, the object is appended to the end of the file, otherwise the
        file will be overwritten (if a file object is given instead of a
        file name, this has no effect).

    """

    if usecPickle and sys.version_info[0] < 3:  # pragma: py2
        import cPickle as pickle
    else:
        import pickle

    if protocol is None:
        protocol = pickle.HIGHEST_PROTOCOL

    if isinstance(fileorname, basestring):
        f = open(fileorname, 'ab' if append else 'wb')
        close = True
    else:
        f = fileorname
        close = False

    try:
        pickle.dump(object, f, protocol=protocol)
    finally:
        if close:
            f.close()


def isiterable(obj):
    """Returns `True` if the given object is iterable."""

    if isinstance(obj, collections.Iterable):
        return True

    try:
        iter(obj)
        return True
    except TypeError:
        return False


class lazyproperty(object):
    """
    Works similarly to property(), but computes the value only once.

    This essentially memoizes the value of the property by storing the result
    of its computation in the ``__dict__`` of the object instance.  This is
    useful for computing the value of some property that should otherwise be
    invariant.  For example::

        >>> class LazyTest(object):
        ...     @lazyproperty
        ...     def complicated_property(self):
        ...         print 'Computing the value for complicated_property..."
        ...         return 42
        ...
        >>> lt = LazyTest()
        >>> lt.complicated_property
        Computing the value for complicated_property...
        42
        >>> lt.complicated_property
        42

    If a setter for this property is defined, it will still be possible to
    manually update the value of the property, if that capability is desired.

    Adapted from the recipe at
    http://code.activestate.com/recipes/363602-lazy-property-evaluation
    """

    def __init__(self, fget, fset=None, fdel=None, doc=None):
        self._fget = fget
        self._fset = fset
        self._fdel = fdel
        if doc is None:
            self.__doc__ = fget.__doc__
        else:
            self.__doc__ = doc

    def __get__(self, obj, owner=None):
        if obj is None:
            return self
        key = self._fget.func_name
        if key not in obj.__dict__:
            val = self._fget(obj)
            obj.__dict__[key] = val
            return val
        else:
            return obj.__dict__[key]

    def __set__(self, obj, val):
        if self._fset:
            self._fset(obj, val)
        obj.__dict__[self._fget.func_name] = val

    def __delete__(self, obj):
        if self._fdel:
            self._fdel(obj)
        key = self._fget.func_name
        if key in obj.__dict__:
            del obj.__dict__[key]

    def getter(self, fget):
        return self.__ter(fget, 0)

    def setter(self, fset):
        return self.__ter(fset, 1)

    def deleter(self, fdel):
        return self.__ter(fdel, 2)

    def __ter(self, f, arg):
        args = [self._fget, self._fset, self._fdel, self.__doc__]
        args[arg] = f
        cls_ns = sys._getframe(1).f_locals
        for k, v in cls_ns.iteritems():
            if v is self:
                property_name = k
                break

        cls_ns[property_name] = lazyproperty(*args)

        return cls_ns[property_name]


# TODO: Provide a class deprecation marker as well.
def deprecated(since, message='', name='', alternative='', pending=False):
    """
    Used to mark a function as deprecated.

    To mark an attribute as deprecated, replace that attribute with a
    depcrecated property.

    Parameters
    ------------
    since : str
        The release at which this API became deprecated.  This is required.

    message : str, optional
        Override the default deprecation message.  The format specifier
        %(func)s may be used for the name of the function, and %(alternative)s
        may be used in the deprecation message to insert the name of an
        alternative to the deprecated function.

    name : str, optional
        The name of the deprecated function; if not provided the name is
        automatically determined from the passed in function, though this is
        useful in the case of renamed functions, where the new function is just
        assigned to the name of the deprecated function.  For example::

            def new_function():
                ...
            oldFunction = new_function

    alternative : str, optional
        An alternative function that the user may use in place of the
        deprecated function.  The deprecation warning will tell the user about
        this alternative if provided.

    pending : bool, optional
        If True, uses a PendingDeprecationWarning instead of a
        DeprecationWarning.
    """

    def deprecate(func, message=message, name=name, alternative=alternative,
                  pending=pending):
        if isinstance(func, classmethod):
            try:
                func = func.__func__
            except AttributeError:
                # classmethods in Python2.6 and below lack the __func__
                # attribute so we need to hack around to get it
                method = func.__get__(None, object)
                if hasattr(method, '__func__'):
                    func = method.__func__
                elif hasattr(method, 'im_func'):
                    func = method.im_func
                else:
                    # Nothing we can do really...  just return the original
                    # classmethod
                    return func
            is_classmethod = True
        else:
            is_classmethod = False

        if not name:
            name = func.__name__

        altmessage = ''
        if not message or type(message) == type(deprecate):
            if pending:
                message = ('The %(func)s function will be deprecated in a '
                           'future version.')
            else:
                message = ('The %(func)s function is deprecated and may '
                           'be removed in a future version.')
            if alternative:
                altmessage = '\n        Use %s instead.' % alternative

        message = ((message % {'func': name, 'alternative': alternative}) +
                   altmessage)

        @functools.wraps(func)
        def deprecated_func(*args, **kwargs):
            if pending:
                category = PendingDeprecationWarning
            else:
                category = DeprecationWarning

            warnings.warn(message, category, stacklevel=2)

            return func(*args, **kwargs)

        old_doc = deprecated_func.__doc__
        if not old_doc:
            old_doc = ''
        old_doc = textwrap.dedent(old_doc).strip('\n')
        altmessage = altmessage.strip()
        if not altmessage:
            altmessage = message.strip()
        new_doc = (('\n.. deprecated:: %(since)s'
                    '\n    %(message)s\n\n' %
                    {'since': since, 'message': altmessage.strip()}) + old_doc)
        if not old_doc:
            # This is to prevent a spurious 'unexected unindent' warning from
            # docutils when the original docstring was blank.
            new_doc += r'\ '

        deprecated_func.__doc__ = new_doc

        if is_classmethod:
            deprecated_func = classmethod(deprecated_func)
        return deprecated_func

    if type(message) == type(deprecate):
        return deprecate(message)

    return deprecate

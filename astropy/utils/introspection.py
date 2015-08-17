# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions related to Python runtime introspection."""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import inspect
import sys
import types

from ..extern import six


__all__ = ['resolve_name', 'minversion', 'find_current_module',
           'isinstancemethod']


__doctest_skip__ = ['find_current_module']


def resolve_name(name, *additional_parts):
    """Resolve a name like ``module.object`` to an object and return it.

    This ends up working like ``from module import object`` but is easier
    to deal with than the `__import__` builtin and supports digging into
    submodules.

    Parameters
    ----------

    name : `str`
        A dotted path to a Python object--that is, the name of a function,
        class, or other object in a module with the full path to that module,
        including parent modules, separated by dots.  Also known as the fully
        qualified name of the object.

    additional_parts : iterable, optional
        If more than one positional arguments are given, those arguments are
        automatically dotted together with ``name``.

    Examples
    --------

    >>> resolve_name('astropy.utils.introspection.resolve_name')
    <function resolve_name at 0x...>
    >>> resolve_name('astropy', 'utils', 'introspection', 'resolve_name')
    <function resolve_name at 0x...>

    Raises
    ------
    `ImportError`
        If the module or named object is not found.
    """

    additional_parts = '.'.join(additional_parts)

    if additional_parts:
        name = name + '.' + additional_parts

    # Note: On python 2 these must be str objects and not unicode
    parts = [str(part) for part in name.split('.')]

    if len(parts) == 1:
        # No dots in the name--just a straight up module import
        cursor = 1
        fromlist=[]
    else:
        cursor = len(parts) - 1
        fromlist = [parts[-1]]

    module_name = parts[:cursor]

    while cursor > 0:
        try:
            ret = __import__(str('.'.join(module_name)), fromlist=fromlist)
            break
        except ImportError:
            if cursor == 0:
                raise
            cursor -= 1
            module_name = parts[:cursor]
            fromlist = [parts[cursor]]
            ret = ''

    for part in parts[cursor:]:
        try:
            ret = getattr(ret, part)
        except AttributeError:
            raise ImportError(name)

    return ret


def minversion(module, version, inclusive=True, version_path='__version__'):
    """
    Returns `True` if the specified Python module satisfies a minimum version
    requirement, and `False` if not.

    By default this uses `pkg_resources.parse_version` to do the version
    comparison if available.  Otherwise it falls back on
    `distutils.version.LooseVersion`.

    Parameters
    ----------

    module : module or `str`
        An imported module of which to check the version, or the name of
        that module (in which case an import of that module is attempted--
        if this fails `False` is returned).

    version : `str`
        The version as a string that this module must have at a minimum (e.g.
        ``'0.12'``).

    inclusive : `bool`
        The specified version meets the requirement inclusively (i.e. ``>=``)
        as opposed to strictly greater than (default: `True`).

    version_path : `str`
        A dotted attribute path to follow in the module for the version.
        Defaults to just ``'__version__'``, which should work for most Python
        modules.

    Examples
    --------

    >>> import astropy
    >>> minversion(astropy, '0.4.4')
    True
    """

    if isinstance(module, types.ModuleType):
        module_name = module.__name__
    elif isinstance(module, six.string_types):
        module_name = module
        try:
            module = resolve_name(module_name)
        except ImportError:
            return False
    else:
        raise ValueError('module argument must be an actual imported '
                         'module, or the import name of the module; '
                         'got {0!r}'.format(module))

    if '.' not in version_path:
        have_version = getattr(module, version_path)
    else:
        have_version = resolve_name('.'.join([module.__name__, version_path]))

    try:
        from pkg_resources import parse_version
    except ImportError:
        from distutils.version import LooseVersion as parse_version

    if inclusive:
        return parse_version(have_version) >= parse_version(version)
    else:
        return parse_version(have_version) > parse_version(version)


def find_current_module(depth=1, finddiff=False):
    """
    Determines the module/package from which this function is called.

    This function has two modes, determined by the ``finddiff`` option. it
    will either simply go the requested number of frames up the call
    stack (if ``finddiff`` is False), or it will go up the call stack until
    it reaches a module that is *not* in a specified set.

    Parameters
    ----------
    depth : int
        Specifies how far back to go in the call stack (0-indexed, so that
        passing in 0 gives back `astropy.utils.misc`).
    finddiff : bool or list
        If False, the returned ``mod`` will just be ``depth`` frames up from
        the current frame. Otherwise, the function will start at a frame
        ``depth`` up from current, and continue up the call stack to the
        first module that is *different* from those in the provided list.
        In this case, ``finddiff`` can be a list of modules or modules
        names. Alternatively, it can be True, which will use the module
        ``depth`` call stack frames up as the module the returned module
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
        If ``finddiff`` is a list with an invalid entry.

    Examples
    --------
    The examples below assume that there are two modules in a package named
    ``pkg``. ``mod1.py``::

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

    frm = inspect.currentframe()
    for i in range(depth):
        frm = frm.f_back
        if frm is None:
            return None

    if finddiff:
        currmod = inspect.getmodule(frm)
        if finddiff is True:
            diffmods = [currmod]
        else:
            diffmods = []
            for fd in finddiff:
                if inspect.ismodule(fd):
                    diffmods.append(fd)
                elif isinstance(fd, six.string_types):
                    diffmods.append(__import__(fd))
                elif fd is True:
                    diffmods.append(currmod)
                else:
                    raise ValueError('invalid entry in finddiff')

        while frm:
            frmb = frm.f_back
            modb = inspect.getmodule(frmb)
            if modb not in diffmods:
                return modb
            frm = frmb
    else:
        return inspect.getmodule(frm)


def find_mod_objs(modname, onlylocals=False):
    """ Returns all the public attributes of a module referenced by name.

    .. note::
        The returned list *not* include subpackages or modules of
        ``modname``, nor does it include private attributes (those that
        begin with '_' or are not in `__all__`).

    Parameters
    ----------
    modname : str
        The name of the module to search.
    onlylocals : bool or list of str
        If `True`, only attributes that are either members of ``modname`` OR
        one of its modules or subpackages will be included. If it is a list
        of strings, those specify the possible packages that will be
        considered "local".

    Returns
    -------
    localnames : list of str
        A list of the names of the attributes as they are named in the
        module ``modname`` .
    fqnames : list of str
        A list of the full qualified names of the attributes (e.g.,
        ``astropy.utils.introspection.find_mod_objs``). For attributes that are
        simple variables, this is based on the local name, but for functions or
        classes it can be different if they are actually defined elsewhere and
        just referenced in ``modname``.
    objs : list of objects
        A list of the actual attributes themselves (in the same order as
        the other arguments)

    """

    __import__(modname)
    mod = sys.modules[modname]

    if hasattr(mod, '__all__'):
        pkgitems = [(k, mod.__dict__[k]) for k in mod.__all__]
    else:
        pkgitems = [(k, mod.__dict__[k]) for k in dir(mod) if k[0] != '_']

    # filter out modules and pull the names and objs out
    ismodule = inspect.ismodule
    localnames = [k for k, v in pkgitems if not ismodule(v)]
    objs = [v for k, v in pkgitems if not ismodule(v)]

    # fully qualified names can be determined from the object's module
    fqnames = []
    for obj, lnm in zip(objs, localnames):
        if hasattr(obj, '__module__') and hasattr(obj, '__name__'):
            fqnames.append(obj.__module__ + '.' + obj.__name__)
        else:
            fqnames.append(modname + '.' + lnm)

    if onlylocals:
        if onlylocals is True:
            onlylocals = [modname]
        valids = [any([fqn.startswith(nm) for nm in onlylocals]) for fqn in fqnames]
        localnames = [e for i, e in enumerate(localnames) if valids[i]]
        fqnames = [e for i, e in enumerate(fqnames) if valids[i]]
        objs = [e for i, e in enumerate(objs) if valids[i]]

    return localnames, fqnames, objs


# Note: I would have preferred call this is_instancemethod, but this naming is
# for consistency with other functions in the `inspect` module
def isinstancemethod(cls, obj):
    """
    Returns `True` if the given object is an instance method of the class
    it is defined on (as opposed to a `staticmethod` or a `classmethod`).

    This requires both the class the object is a member of as well as the
    object itself in order to make this determination.

    Parameters
    ----------
    cls : `type`
        The class on which this method was defined.
    obj : `object`
        A member of the provided class (the membership is not checked directly,
        but this function will always return `False` if the given object is not
        a member of the given class).

    Examples
    --------
    >>> from astropy.extern import six
    >>> class MetaClass(type):
    ...     def a_classmethod(cls): pass
    ...
    >>> @six.add_metaclass(MetaClass)
    ... class MyClass(object):
    ...     __metaclass__ = MetaClass
    ...     def an_instancemethod(self): pass
    ...     @classmethod
    ...     def another_classmethod(cls): pass
    ...     @staticmethod
    ...     def a_staticmethod(): pass
    ...
    >>> isinstancemethod(MyClass, MyClass.a_classmethod)
    False
    >>> isinstancemethod(MyClass, MyClass.another_classmethod)
    False
    >>> isinstancemethod(MyClass, MyClass.a_staticmethod)
    False
    >>> isinstancemethod(MyClass, MyClass.an_instancemethod)
    True
    """

    return _isinstancemethod(cls, obj)


if six.PY3:
    def _isinstancemethod(cls, obj):
        if not isinstance(obj, types.FunctionType):
            return False

        # Unfortunately it seems the easiest way to get to the original
        # staticmethod object is to look in the class's __dict__, though we
        # also need to look up the MRO in case the method is not in the given
        # class's dict
        name = obj.__name__
        for basecls in cls.mro():  # This includes cls
            if name in basecls.__dict__:
                return not isinstance(basecls.__dict__[name], staticmethod)

        # This shouldn't happen, though this is the most sensible response if
        # it does.
        raise AttributeError(name)
else:
    def _isinstancemethod(cls, obj):
        return isinstance(obj, types.MethodType) and obj.im_class is cls


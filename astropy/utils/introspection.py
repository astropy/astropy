# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions related to Python runtime introspection."""

from __future__ import annotations

import importlib
import inspect
import os
import sys
from importlib import metadata
from importlib.metadata import packages_distributions
from typing import TYPE_CHECKING

from packaging.version import Version

from .decorators import deprecated

if TYPE_CHECKING:
    from types import FrameType, ModuleType
    from typing import Literal

__all__ = ["resolve_name", "minversion", "find_current_module", "isinstancemethod"]

__doctest_skip__ = ["find_current_module"]


def resolve_name(name: str, *additional_parts: str) -> object:
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
    additional_parts = ".".join(additional_parts)

    if additional_parts:
        name = name + "." + additional_parts

    parts = name.split(".")

    if len(parts) == 1:
        # No dots in the name--just a straight up module import
        cursor = 1
        fromlist = []
    else:
        cursor = len(parts) - 1
        fromlist = [parts[-1]]

    module_name = parts[:cursor]

    while cursor > 0:
        try:
            ret = __import__(".".join(module_name), fromlist=fromlist)
            break
        except ImportError:
            if cursor == 0:
                raise
            cursor -= 1
            module_name = parts[:cursor]
            fromlist = [parts[cursor]]
            ret = ""

    for part in parts[cursor:]:
        try:
            ret = getattr(ret, part)
        except AttributeError:
            raise ImportError(name)

    return ret


def minversion(module: ModuleType | str, version: str, inclusive: bool = True) -> bool:
    """
    Returns `True` if the specified Python module satisfies a minimum version
    requirement, and `False` if not.

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

    Examples
    --------
    >>> import astropy
    >>> minversion(astropy, '0.4.4')
    True
    """
    if inspect.ismodule(module):
        module_name = module.__name__
        module_version = getattr(module, "__version__", None)
    elif isinstance(module, str):
        module_name = module
        module_version = None
        try:
            module = resolve_name(module_name)
        except ImportError:
            return False
    else:
        raise ValueError(
            "module argument must be an actual imported "
            "module, or the import name of the module; "
            f"got {repr(module)}"
        )

    if module_version is None:
        try:
            module_version = metadata.version(module_name)
        except metadata.PackageNotFoundError:
            # Maybe the distribution name is different from package name.
            # Calling packages_distributions is costly so we do it only
            # if necessary, as only a few packages don't have the same
            # distribution name.
            dist_names = packages_distributions()
            module_version = metadata.version(dist_names[module_name][0])

    if inclusive:
        return Version(module_version) >= Version(version)
    else:
        return Version(module_version) > Version(version)


def find_current_module(
    depth: int = 1, finddiff: bool | list[Literal[True] | str | ModuleType] = False
) -> ModuleType | None:
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
        currmod = _get_module_from_frame(frm)
        if finddiff is True:
            diffmods = [currmod]
        else:
            diffmods = []
            for fd in finddiff:
                if inspect.ismodule(fd):
                    diffmods.append(fd)
                elif isinstance(fd, str):
                    diffmods.append(importlib.import_module(fd))
                elif fd is True:
                    diffmods.append(currmod)
                else:
                    raise ValueError("invalid entry in finddiff")

        while frm:
            frmb = frm.f_back
            modb = _get_module_from_frame(frmb)
            if modb not in diffmods:
                return modb
            frm = frmb
    else:
        return _get_module_from_frame(frm)


def _get_module_from_frame(frm: FrameType) -> ModuleType | None:
    """Uses inspect.getmodule() to get the module that the current frame's
    code is running in.

    However, this does not work reliably for code imported from a zip file,
    so this provides a fallback mechanism for that case which is less
    reliable in general, but more reliable than inspect.getmodule() for this
    particular case.
    """
    mod = inspect.getmodule(frm)
    if mod is not None:
        return mod

    # Check to see if we're importing from a bundle file. First ensure that
    # __file__ is available in globals; this is cheap to check to bail out
    # immediately if this fails

    if "__file__" in frm.f_globals and "__name__" in frm.f_globals:
        filename = frm.f_globals["__file__"]

        # Using __file__ from the frame's globals and getting it into the form
        # of an absolute path name with .py at the end works pretty well for
        # looking up the module using the same means as inspect.getmodule

        if filename[-4:].lower() in (".pyc", ".pyo"):
            filename = filename[:-4] + ".py"
        filename = os.path.realpath(os.path.abspath(filename))
        if filename in inspect.modulesbyfile:
            return sys.modules.get(inspect.modulesbyfile[filename])

        # On Windows, inspect.modulesbyfile appears to have filenames stored
        # in lowercase, so we check for this case too.
        if filename.lower() in inspect.modulesbyfile:
            return sys.modules.get(inspect.modulesbyfile[filename.lower()])

    # Otherwise there are still some even trickier things that might be possible
    # to track down the module, but we'll leave those out unless we find a case
    # where it's really necessary.  So return None if the module is not found.
    return None


@deprecated(since="6.1")
def find_mod_objs(modname, onlylocals=False):
    """Returns all the public attributes of a module referenced by name.

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
    mod = resolve_name(modname)

    if hasattr(mod, "__all__"):
        pkgitems = [(k, mod.__dict__[k]) for k in mod.__all__]
    else:
        pkgitems = [(k, mod.__dict__[k]) for k in dir(mod) if k[0] != "_"]

    # filter out modules and pull the names and objs out
    ismodule = inspect.ismodule
    localnames = [k for k, v in pkgitems if not ismodule(v)]
    objs = [v for k, v in pkgitems if not ismodule(v)]

    # fully qualified names can be determined from the object's module
    fqnames = []
    for obj, lnm in zip(objs, localnames):
        if hasattr(obj, "__module__") and hasattr(obj, "__name__"):
            fqnames.append(obj.__module__ + "." + obj.__name__)
        else:
            fqnames.append(modname + "." + lnm)

    if onlylocals:
        if onlylocals is True:
            onlylocals = [modname]
        valids = [any(fqn.startswith(nm) for nm in onlylocals) for fqn in fqnames]
        localnames = [e for i, e in enumerate(localnames) if valids[i]]
        fqnames = [e for i, e in enumerate(fqnames) if valids[i]]
        objs = [e for i, e in enumerate(objs) if valids[i]]

    return localnames, fqnames, objs


# Note: I would have preferred call this is_instancemethod, but this naming is
# for consistency with other functions in the `inspect` module
@deprecated(since="6.1")
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
    """
    if not inspect.isfunction(obj):
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

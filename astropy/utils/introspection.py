# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions related to Python runtime introspection."""

import inspect
import os
import sys
from importlib import import_module, metadata
from importlib.metadata import packages_distributions
from types import FrameType, ModuleType
from typing import Literal

from packaging.version import Version

from .decorators import deprecated

__all__ = ["find_current_module", "minversion", "resolve_name"]

__doctest_skip__ = ["find_current_module"]


@deprecated(
    since="7.0", alternative="importlib (e.g. importlib.import_module for modules)"
)
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
    >>> import numpy
    >>> minversion(numpy, '1.21.0')
    True
    """
    if inspect.ismodule(module):
        module_name = module.__name__
        module_version = getattr(module, "__version__", None)
    elif isinstance(module, str):
        module_name = module
        module_version = None
        try:
            module = import_module(module_name)
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
    for _ in range(depth):
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
                    diffmods.append(import_module(fd))
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

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Simple utility functions and bug fixes for compatibility with all supported
versions of Python.  This module should generally not be used directly, as
everything in `__all__` will be imported into `astropy.utils.compat` and can
be accessed from there.
"""

import sys
import functools
from contextlib import suppress


__all__ = ['override__dir__', 'suppress',
           'possible_filename', 'namedtuple_asdict']


def possible_filename(filename):
    """
    Determine if the ``filename`` argument is an allowable type for a filename.

    In Python 3.3 use of non-unicode filenames on system calls such as
    `os.stat` and others that accept a filename argument was deprecated (and
    may be removed outright in the future).

    Therefore this returns `True` in all cases except for `bytes` strings in
    Windows.
    """

    if isinstance(filename, str):
        return True
    elif isinstance(filename, bytes):
        return not (sys.platform == 'win32')

    return False


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
    # http://bugs.python.org/issue12166

    @functools.wraps(f)
    def override__dir__wrapper(self):
        members = set(object.__dir__(self))
        members.update(f(self))
        return sorted(members)

    return override__dir__wrapper


def namedtuple_asdict(namedtuple):
    """
    The same as ``namedtuple._adict()``.

    Parameters
    ----------
    namedtuple : collections.namedtuple
    The named tuple to get the dict of
    """
    return namedtuple._asdict()

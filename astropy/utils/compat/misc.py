# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Simple utility functions and bug fixes for compatibility with all supported
versions of Python.  This module should generally not be used directly, as
everything in `__all__` will be imported into `astropy.utils.compat` and can
be accessed from there.

Includes the following fixes:

* The `contextlib.ignored` context manager, which is only available in Python
  3.4 or greater.

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six

import functools
import sys


__all__ = ['invalidate_caches', 'override__dir__', 'ignored',
           'possible_filename']


def possible_filename(filename):
    """
    Determine if the ``filename`` argument is an allowable type for a filename.

    In Python 3.3 use of non-unicode filenames on system calls such as
    `os.stat` and others that accept a filename argument was deprecated (and
    may be removed outright in the future).

    Therefore this returns `True` in all cases except for `bytes` strings in
    Windows on Python >= 3.3.
    """

    if isinstance(filename, six.text_type):
        return True
    elif isinstance(filename, six.binary_type):
        return not (sys.platform == 'win32' and
                    sys.version_info[:2] >= (3, 3))

    return False


# Python 3.3's importlib caches filesystem reads for faster imports in the
# general case. But sometimes it's necessary to manually invalidate those
# caches so that the import system can pick up new generated files.  See
# https://github.com/astropy/astropy/issues/820
if sys.version_info[:2] >= (3, 3):
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
    if sys.version_info[:2] < (3, 3):
        # There was no straightforward way to do this until Python 3.3, so
        # we have this complex monstrosity
        @functools.wraps(f)
        def override__dir__wrapper(self):
            members = set()
            for cls in self.__class__.mro():
                members.update(dir(cls))
            members.update(six.iterkeys(self.__dict__))
            members.update(f(self))
            return sorted(members)
    else:
        # http://bugs.python.org/issue12166

        @functools.wraps(f)
        def override__dir__wrapper(self):
            members = set(object.__dir__(self))
            members.update(f(self))
            return sorted(members)

    return override__dir__wrapper


try:
    from contextlib import ignored
except ImportError:
    from contextlib import contextmanager
    @contextmanager
    def ignored(*exceptions):
        """A context manager for ignoring exceptions.  Equivalent to::

            try:
                <body>
            except exceptions:
                pass

        Example::

            >>> import os
            >>> with ignored(OSError):
            ...     os.remove('file-that-does-not-exist')

        """

        try:
            yield
        except exceptions:
            pass

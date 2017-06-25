# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Simple utility functions and bug fixes for compatibility with all supported
versions of Python.  This module should generally not be used directly, as
everything in `__all__` will be imported into `astropy.utils.compat` and can
be accessed from there.

Includes the following fixes:

* The `contextlib.suppress` context manager, which is only available in Python
  3.4 or greater.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six

import functools
import sys

__all__ = ['invalidate_caches', 'override__dir__', 'suppress',
           'possible_filename', 'namedtuple_asdict']


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
    from contextlib import suppress
except ImportError:
    from contextlib import contextmanager

    @contextmanager
    def suppress(*exceptions):
        """A context manager for ignoring exceptions.  Equivalent to::

            try:
                <body>
            except exceptions:
                pass

        Example::

            >>> import os
            >>> with suppress(OSError):
            ...     os.remove('file-that-does-not-exist')

        """
        try:
            yield
        except exceptions:
            pass


# For unclear reasons, the `_asdict` method of namedtuple produces an empty
# dictionary if the namedtuple is a subclass of another namedtuple... But
# *only* in py 3.3.  >3.4 or 2.7 seem to work just fine.  So we provide this
# for compatibility as long as 3.3 is supported.
# Furthermore, in python 3.4.x except for 3.4.4, `_asdict` produces only a
# *partial* dictionary.  So we work around that case too.
if sys.version_info[0] == 3 and sys.version_info[:3] < (3, 4, 4):
    def namedtuple_asdict(namedtuple):
        """
        The same as ``namedtuple._adict()``, but fixed to work even when
        namedtuple is a subclass of another namedtuple

        Parameters
        ----------
        namedtuple : collections.namedtuple
            The named tuple to get the dict of
        """
        return {fi: getattr(namedtuple, fi) for fi in namedtuple._fields}
else:
    def namedtuple_asdict(namedtuple):
        """
        The same as ``namedtuple._adict()``.

        Parameters
        ----------
        namedtuple : collections.namedtuple
            The named tuple to get the dict of
        """
        return namedtuple._asdict()

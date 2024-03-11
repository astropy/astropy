# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Simple utility functions and bug fixes for compatibility with all supported
versions of Python.  This module should generally not be used directly, as
everything in `__all__` will be imported into `astropy.utils.compat` and can
be accessed from there.
"""

import functools

from astropy.utils.decorators import deprecated

__all__ = ["override__dir__"]


@deprecated(
    since="5.2",
    message="http://bugs.python.org/issue12166 is resolved, {alternative}.",
    alternative="use ``sorted(super().__dir__() + ...)`` instead",
)
def override__dir__(f):
    """
    When overriding a __dir__ method on an object, you often want to include the
    "standard" members on the object as well.  This decorator takes care of that
    automatically, and all the wrapped function needs to do is return a list of
    the "special" members that wouldn't be found by the normal Python means.

    Example
    -------

    Your class could define __dir__ as follows::

        @override__dir__
        def __dir__(self):
            return ['special_method1', 'special_method2']

    Notes
    -----
    This function was introduced because of http://bugs.python.org/issue12166,
    which has since been resolved by
    http://hg.python.org/cpython/rev/8f403199f999. Now, the best way to
    customize ``__dir__`` is to use ``super``.
    ::

        def __dir__(self):
            added = {'special_method1', 'special_method2'}
            return sorted(set(super().__dir__()) | added)
    """
    # http://bugs.python.org/issue12166

    @functools.wraps(f)
    def override__dir__wrapper(self):
        members = set(object.__dir__(self))
        members.update(f(self))
        return sorted(members)

    return override__dir__wrapper

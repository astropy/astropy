# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Supporting definitions for the Python regression tests."""

import contextlib
import sys
import warnings
import re


__all__ = ["check_py3k_warnings"]


class WarningsRecorder(object):
    """Convenience wrapper for the warnings list returned on
       entry to the warnings.catch_warnings() context manager.
    """
    def __init__(self, warnings_list):
        self._warnings = warnings_list
        self._last = 0

    def __getattr__(self, attr):
        if len(self._warnings) > self._last:
            return getattr(self._warnings[-1], attr)
        elif attr in warnings.WarningMessage._WARNING_DETAILS:
            return None
        raise AttributeError("%r has no attribute %r" % (self, attr))

    @property
    def warnings(self):
        return self._warnings[self._last:]

    def reset(self):
        self._last = len(self._warnings)


def _filterwarnings(filters, quiet=False):
    """Catch the warnings, then check if all the expected
    warnings have been raised and re-raise unexpected warnings.
    If 'quiet' is True, only re-raise the unexpected warnings.
    """
    # Clear the warning registry of the calling module
    # in order to re-raise the warnings.
    frame = sys._getframe(2)
    registry = frame.f_globals.get('__warningregistry__')
    if registry:
        registry.clear()
    with warnings.catch_warnings(record=True) as w:
        # Set filter "always" to record all warnings.  Because
        # test_warnings swap the module, we need to look up in
        # the sys.modules dictionary.
        sys.modules['warnings'].simplefilter("always")
        yield WarningsRecorder(w)
    # Filter the recorded warnings
    reraise = [warning.message for warning in w]
    missing = []
    for msg, cat in filters:
        seen = False
        for exc in reraise[:]:
            message = str(exc)
            # Filter out the matching messages
            if (re.match(msg, message, re.I) and
                issubclass(exc.__class__, cat)):
                seen = True
                reraise.remove(exc)
        if not seen and not quiet:
            # This filter caught nothing
            missing.append((msg, cat.__name__))
    if reraise:
        raise AssertionError("unhandled warning %r" % reraise[0])
    if missing:
        raise AssertionError("filter (%r, %s) did not catch any warning" %
                             missing[0])


@contextlib.contextmanager
def check_py3k_warnings(*filters, **kwargs):
    """Context manager to silence py3k warnings.

    Accept 2-tuples as positional arguments:
        ("message regexp", WarningCategory)

    Optional argument:
     - if 'quiet' is True, it does not fail if a filter catches nothing
        (default False)

    Without argument, it defaults to:
        check_py3k_warnings(("", DeprecationWarning), quiet=False)
    """
    if sys.py3kwarning:
        if not filters:
            filters = (("", DeprecationWarning),)
    else:
        # It should not raise any py3k warning
        filters = ()
    return _filterwarnings(filters, kwargs.get('quiet'))

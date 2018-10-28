# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


import contextlib
import warnings
from ...tests.helper import catch_warnings


@contextlib.contextmanager
def ignore_non_integer_warning():
    # We need to ignore this warning on Scipy < 0.14.
    # When our minimum version of Scipy is bumped up, this can be
    # removed.
    with catch_warnings():
        warnings.filterwarnings(
            "always", "using a non-integer number instead of an integer "
            "will result in an error in the future", DeprecationWarning)
        yield

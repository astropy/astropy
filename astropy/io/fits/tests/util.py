# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test utility functions."""

import warnings

from ....tests.helper import catch_warnings


class ignore_warnings(catch_warnings):
    def __enter__(self):
        retval = super(ignore_warnings, self).__enter__()
        warnings.simplefilter('ignore')
        return retval

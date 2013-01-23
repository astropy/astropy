# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test utility functions."""

import warnings


class ignore_warnings(warnings.catch_warnings):
    def __enter__(self):
        retval = super(ignore_warnings, self).__enter__()
        warnings.simplefilter('ignore')
        return retval

"""Test utility functions."""

import warnings


class ignore_warnings(catch_warnings):
    def __enter__(self):
        retval = super(ignore_warnings, self).__enter__()
        warnings.simplefilter('ignore')
        return retval

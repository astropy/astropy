# Licensed under a 3-clause BSD style license - see LICENSE.rst


# This module is not a test module, but is used as part of the
# test_run_after_2to3 test in test_run_tests.py
def test_run_after_2to3():
    try:
        1 / 0
    except ZeroDivisionError, e:
        pass

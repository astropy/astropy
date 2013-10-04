# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# test helper.run_tests function

from .. import helper
from ... import _get_test_runner


# run_tests should raise ValueError when asked to run on a module it can't find
def test_module_not_found():
    with helper.pytest.raises(ValueError):
        _get_test_runner().run_tests('fake.module')


# run_tests should raise ValueError when passed an invalid pastebin= option
def test_pastebin_keyword():
    with helper.pytest.raises(ValueError):
        _get_test_runner().run_tests(pastebin='not_an_option')


# tests that tests are only run in Python 3 out of the 2to3'd build (otherwise
# a syntax error would occur)
try:
    from .run_after_2to3 import test_run_after_2to3
except SyntaxError:
    def test_run_after_2to3():
        helper.pytest.fail("Not running the 2to3'd tests!")

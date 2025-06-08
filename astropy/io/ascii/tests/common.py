# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from pathlib import Path

__all__ = [
    "setup_function",
    "teardown_function",
]

CWD = Path.cwd()
TEST_DIR = Path(__file__).parent


def setup_function(function):
    os.chdir(TEST_DIR)


def teardown_function(function):
    os.chdir(CWD)


def assert_equal_splitlines(arg1, arg2):
    __tracebackhide__ = True
    assert arg1.splitlines() == arg2.splitlines()

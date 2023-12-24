# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

__all__ = [
    "setup_function",
    "teardown_function",
]

CWD = os.getcwd()
TEST_DIR = os.path.dirname(__file__)


def setup_function(function):
    os.chdir(TEST_DIR)


def teardown_function(function):
    os.chdir(CWD)

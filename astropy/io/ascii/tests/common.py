# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import numpy as np

__all__ = [
    "assert_equal",
    "assert_almost_equal",
    "assert_true",
    "setup_function",
    "teardown_function",
]

CWD = os.getcwd()
TEST_DIR = os.path.dirname(__file__)


def setup_function(function):
    os.chdir(TEST_DIR)


def teardown_function(function):
    os.chdir(CWD)


# Compatibility functions to convert from nose to pytest
def assert_equal(a, b):
    assert a == b


def assert_almost_equal(a, b, **kwargs):
    assert np.allclose(a, b, **kwargs)


def assert_true(a):
    assert a

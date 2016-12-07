# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import

import os


import numpy as np


__all__ = ['raises', 'assert_equal', 'assert_almost_equal',
           'assert_true', 'setup_function', 'teardown_function',
           'has_isnan']

CWD = os.getcwd()
TEST_DIR = os.path.dirname(__file__)

has_isnan = True
try:
    from math import isnan  # pylint: disable=W0611
except ImportError:
    try:
        from numpy import isnan  # pylint: disable=W0611
    except ImportError:
        has_isnan = False
        print('Tests requiring isnan will fail')


def setup_function(function):
    os.chdir(TEST_DIR)


def teardown_function(function):
    os.chdir(CWD)


# Compatibility functions to convert from nose to py.test
def assert_equal(a, b):
    assert a == b


def assert_almost_equal(a, b, **kwargs):
    assert np.allclose(a, b, **kwargs)


def assert_true(a):
    assert a


def make_decorator(func):
    """
    Wraps a test decorator so as to properly replicate metadata
    of the decorated function, including nose's additional stuff
    (namely, setup and teardown).
    """
    def decorate(newfunc):
        if hasattr(func, 'compat_func_name'):
            name = func.compat_func_name
        else:
            name = func.__name__
        newfunc.__dict__ = func.__dict__
        newfunc.__doc__ = func.__doc__
        newfunc.__module__ = func.__module__
        if not hasattr(newfunc, 'compat_co_firstlineno'):
            try:
                newfunc.compat_co_firstlineno = func.func_code.co_firstlineno
            except AttributeError:
                newfunc.compat_co_firstlineno = func.__code__.co_firstlineno
        try:
            newfunc.__name__ = name
        except TypeError:
            # can't set func name in 2.3
            newfunc.compat_func_name = name
        return newfunc
    return decorate


def raises(*exceptions):
    """Test must raise one of expected exceptions to pass.

    Example use::

      @raises(TypeError, ValueError)
      def test_raises_type_error():
          raise TypeError("This test passes")

      @raises(Exception)
      def test_that_fails_by_passing():
          pass

    If you want to test many assertions about exceptions in a single test,
    you may want to use `assert_raises` instead.
    """
    valid = ' or '.join([e.__name__ for e in exceptions])

    def decorate(func):
        name = func.__name__

        def newfunc(*arg, **kw):
            try:
                func(*arg, **kw)
            except exceptions:
                pass
            else:
                message = "{}() did not raise {}".format(name, valid)
                raise AssertionError(message)
        newfunc = make_decorator(func)(newfunc)
        return newfunc
    return decorate

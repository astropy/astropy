import inspect

from ...extern import six
from ...tests.helper import pytest


def test_argparse():
    from ..compat import argparse


@pytest.mark.skipif("sys.version_info >= (3, 4)")
def test_wraps():
    """
    Tests the compatibility replacement for functools.wraps which supports
    argument preservation across all supported Python versions.
    """

    from ..compat import wraps

    def foo(a, b, c=1, d=2, e=3, **kwargs):
        """A test function."""

        return a, b, c, d, e, kwargs

    @wraps(foo)
    def bar(*args, **kwargs):
        return ('test',) + foo(*args, **kwargs)

    expected = ('test', 1, 2, 3, 4, 5, {'f': 6, 'g': 7})
    assert bar(1, 2, 3, 4, 5, f=6, g=7) == expected
    assert bar.__name__ == 'foo'
    assert bar.__doc__ == "A test function."

    if hasattr(foo, '__qualname__'):
        assert bar.__qualname__ == foo.__qualname__

    if six.PY2:
        argspec = inspect.getargspec(bar)
        assert argspec.keywords == 'kwargs'
    else:
        argspec = inspect.getfullargspec(bar)
        assert argspec.varkw == 'kwargs'

    assert argspec.args == ['a', 'b', 'c', 'd', 'e']
    assert argspec.defaults == (1, 2, 3)

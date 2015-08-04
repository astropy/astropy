# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import functools
import inspect
import pickle

from ..decorators import (deprecated_attribute, deprecated, wraps,
                          sharedmethod, classproperty)
from ..exceptions import AstropyDeprecationWarning
from ...extern import six
from ...tests.helper import pytest, catch_warnings


def test_wraps():
    """
    Tests the compatibility replacement for functools.wraps which supports
    argument preservation across all supported Python versions.
    """

    def foo(a, b, c=1, d=2, e=3, **kwargs):
        """A test function."""

        return a, b, c, d, e, kwargs

    @wraps(foo)
    def bar(*args, **kwargs):
        return ('test',) + foo(*args, **kwargs)

    expected = ('test', 1, 2, 3, 4, 5, {'f': 6, 'g': 7})
    assert bar(1, 2, 3, 4, 5, f=6, g=7) == expected
    assert bar.__name__ == 'foo'

    if foo.__doc__ is not None:
        # May happen if using optimized opcode
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


def test_wraps_exclude_names():
    """
    Test the optional ``exclude_names`` argument to the wraps decorator.
    """

    # This particular test demonstrates wrapping an instance method
    # as a function and excluding the "self" argument:

    class TestClass(object):
        def method(self, a, b, c=1, d=2, **kwargs):
            return (a, b, c, d, kwargs)

    test = TestClass()

    @wraps(test.method, exclude_args=('self',))
    def func(*args, **kwargs):
        return test.method(*args, **kwargs)

    argspec = inspect.getargspec(func)
    assert argspec.args == ['a', 'b', 'c', 'd']

    assert func('a', 'b', e=3) == ('a', 'b', 1, 2, {'e': 3})


def test_wraps_keep_orig_name():
    """
    Test that when __name__ is excluded from the ``assigned`` argument
    to ``wrap`` that the function being wrapped keeps its original name.

    Regression test for https://github.com/astropy/astropy/pull/4016
    """

    def foo():
        pass

    assigned = list(functools.WRAPPER_ASSIGNMENTS)
    assigned.remove('__name__')

    def bar():
        pass

    orig_bar = bar

    bar = wraps(foo, assigned=assigned)(bar)

    assert bar is not orig_bar
    assert bar.__name__ == 'bar'


def test_deprecated_attribute():
    class DummyClass:
        def __init__(self):
            self._foo = 42

        def set_private(self):
            self._foo = 100

        foo = deprecated_attribute('foo', '0.2')

    dummy = DummyClass()

    with catch_warnings(AstropyDeprecationWarning) as w:
        x = dummy.foo

    assert len(w) == 1
    assert str(w[0].message) == ("The foo attribute is deprecated and may be "
                                 "removed in a future version.")

    with catch_warnings() as w:
        dummy.set_private()

    assert len(w) == 0


# This needs to be defined outside of the test function, because we
# want to try to pickle it.
@deprecated('100.0')
class TestA(object):
    """
    This is the class docstring.
    """
    def __init__(self):
        """
        This is the __init__ docstring
        """
        pass


class TestMeta(type):
    metaclass_attr = 1


@deprecated('100.0')
@six.add_metaclass(TestMeta)
class TestB(object):
    pass


def test_deprecated_class():
    orig_A = TestA.__bases__[0]

    # The only thing that should be different about the new class
    # is __doc__, __init__, __bases__ and __subclasshook__.
    for x in dir(orig_A):
        if x not in ('__doc__', '__init__', '__bases__', '__dict__',
                     '__subclasshook__'):
            assert getattr(TestA, x) == getattr(orig_A, x)

    with catch_warnings(AstropyDeprecationWarning) as w:
        TestA()

    assert len(w) == 1
    if TestA.__doc__ is not None:
        assert 'function' not in TestA.__doc__
        assert 'deprecated' in TestA.__doc__
        assert 'function' not in TestA.__init__.__doc__
        assert 'deprecated' in TestA.__init__.__doc__

    # Make sure the object is picklable
    pickle.dumps(TestA)


def test_deprecated_class_with_super():
    """
    Regression test for an issue where classes that used `super()` in their
    ``__init__`` did not actually call the correct class's ``__init__`` in the
    MRO.
    """

    @deprecated('100.0')
    class TestB(object):
        def __init__(self, a, b):
            super(TestB, self).__init__()

    with catch_warnings(AstropyDeprecationWarning) as w:
        TestB(1, 2)

    assert len(w) == 1
    if TestB.__doc__ is not None:
        assert 'function' not in TestB.__doc__
        assert 'deprecated' in TestB.__doc__
        assert 'function' not in TestB.__init__.__doc__
        assert 'deprecated' in TestB.__init__.__doc__


def test_deprecated_class_with_custom_metaclass():
    """
    Regression test for an issue where deprecating a class with a metaclass
    other than type did not restore the metaclass properly (at least not
    on Python 3).
    """

    with catch_warnings(AstropyDeprecationWarning) as w:
        TestB()

    assert len(w) == 1
    assert type(TestB) is TestMeta
    assert TestB.metaclass_attr == 1


def test_deprecated_static_and_classmethod():
    """
    Regression test for issue introduced by
    https://github.com/astropy/astropy/pull/2811 and mentioned also here:
    https://github.com/astropy/astropy/pull/2580#issuecomment-51049969
    where it appears that deprecated staticmethods didn't work on Python 2.6.
    """

    class A(object):
        """Docstring"""

        @deprecated('1.0')
        @staticmethod
        def B():
            pass

        @deprecated('1.0')
        @classmethod
        def C(cls):
            pass

    with catch_warnings(AstropyDeprecationWarning) as w:
        A.B()

    assert len(w) == 1
    if A.__doc__ is not None:
        assert 'deprecated' in A.B.__doc__

    with catch_warnings(AstropyDeprecationWarning) as w:
        A.C()

    assert len(w) == 1
    if A.__doc__ is not None:
        assert 'deprecated' in A.C.__doc__


@pytest.mark.skipif('six.PY3')
def test_sharedmethod_imfunc():
    """
    Test that the im_func of a sharedmethod always points to the correct
    underlying function.

    This only applies to Python 2 as Python 3 does not have an im_func
    attribute on methods.
    """

    # The original function
    def foo(self): pass
    actual_foo = foo

    class Bar(object):
        foo = sharedmethod(actual_foo)

    assert Bar.foo.im_func is actual_foo
    assert Bar().foo.im_func is actual_foo

    # Now test the case where there the metaclass has a separate
    # implementation
    def foo(cls): pass
    actual_foo_2 = foo

    class MetaBar(type):
        foo = actual_foo_2

    class Bar(object):
        __metaclass__ = MetaBar

        foo = sharedmethod(actual_foo)

    assert Bar.foo.im_func is actual_foo_2
    assert Bar().foo.im_func is actual_foo

    # Finally, test case where the metaclass also has an attribute called
    # 'foo', but it is not a method (hence sharedmethod should ignore it)
    class MetaBar(type):
        foo = None

    class Bar(object):
        __metaclass__ = MetaBar

        foo = sharedmethod(actual_foo)

    assert Bar.foo.im_func is actual_foo
    assert Bar().foo.im_func is actual_foo


def test_sharedmethod_reuse_on_subclasses():
    """
    Regression test for an issue where sharedmethod would bind to one class
    for all time, causing the same method not to work properly on other
    subclasses of that class.

    It has the same problem when the same sharedmethod is called on different
    instances of some class as well.
    """

    class AMeta(type):
        def foo(cls):
            return cls.x

    six.add_metaclass(AMeta)
    class A(object):
        x = 3

        def __init__(self, x):
            self.x = x

        @sharedmethod
        def foo(self):
            return self.x

    a1 = A(1)
    a2 = A(2)

    assert a1.foo() == 1
    assert a2.foo() == 2

    # Similar test now, but for multiple subclasses using the same sharedmethod
    # as a classmethod
    assert A.foo() == 3

    class B(A):
        x = 5

    assert B.foo() == 5


def test_classproperty_docstring():
    """
    Tests that the docstring is set correctly on classproperties.

    This failed previously due to a bug in Python that didn't always
    set __doc__ properly on instances of property subclasses.
    """

    class A(object):
        # Inherits docstring from getter
        @classproperty
        def foo(cls):
            """The foo."""

            return 1

    assert A.__dict__['foo'].__doc__ == "The foo."

    class B(object):
        # Use doc passed to classproperty constructor
        def _get_foo(cls): return 1

        foo = classproperty(_get_foo, doc="The foo.")

    assert B.__dict__['foo'].__doc__ == "The foo."

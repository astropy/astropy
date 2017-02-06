# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import functools
import inspect
import pickle

from ..decorators import (deprecated_attribute, deprecated, wraps,
                          sharedmethod, classproperty,
                          format_doc, deprecated_renamed_argument)
from ..exceptions import AstropyDeprecationWarning, AstropyUserWarning
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

    if six.PY2:
        argspec = inspect.getargspec(func)
    else:
        argspec = inspect.getfullargspec(func)
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
class TA(object):
    """
    This is the class docstring.
    """
    def __init__(self):
        """
        This is the __init__ docstring
        """
        pass


class TMeta(type):
    metaclass_attr = 1


@deprecated('100.0')
@six.add_metaclass(TMeta)
class TB(object):
    pass


def test_deprecated_class():
    orig_A = TA.__bases__[0]

    # The only thing that should be different about the new class
    # is __doc__, __init__, __bases__ and __subclasshook__.
    # and __init_subclass__ for Python 3.6+.
    for x in dir(orig_A):
        if x not in ('__doc__', '__init__', '__bases__', '__dict__',
                     '__subclasshook__', '__init_subclass__'):
            assert getattr(TA, x) == getattr(orig_A, x)

    with catch_warnings(AstropyDeprecationWarning) as w:
        TA()

    assert len(w) == 1
    if TA.__doc__ is not None:
        assert 'function' not in TA.__doc__
        assert 'deprecated' in TA.__doc__
        assert 'function' not in TA.__init__.__doc__
        assert 'deprecated' in TA.__init__.__doc__

    # Make sure the object is picklable
    pickle.dumps(TA)


def test_deprecated_class_with_super():
    """
    Regression test for an issue where classes that used `super()` in their
    ``__init__`` did not actually call the correct class's ``__init__`` in the
    MRO.
    """

    @deprecated('100.0')
    class TB(object):
        def __init__(self, a, b):
            super(TB, self).__init__()

    with catch_warnings(AstropyDeprecationWarning) as w:
        TB(1, 2)

    assert len(w) == 1
    if TB.__doc__ is not None:
        assert 'function' not in TB.__doc__
        assert 'deprecated' in TB.__doc__
        assert 'function' not in TB.__init__.__doc__
        assert 'deprecated' in TB.__init__.__doc__


def test_deprecated_class_with_custom_metaclass():
    """
    Regression test for an issue where deprecating a class with a metaclass
    other than type did not restore the metaclass properly (at least not
    on Python 3).
    """

    with catch_warnings(AstropyDeprecationWarning) as w:
        TB()

    assert len(w) == 1
    assert type(TB) is TMeta
    assert TB.metaclass_attr == 1


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


def test_deprecated_argument():
    # Tests the decorator with function, method, staticmethod and classmethod.

    class Test(object):

        @classmethod
        @deprecated_renamed_argument('clobber', 'overwrite', '1.3')
        def test1(cls, overwrite):
            return overwrite

        @staticmethod
        @deprecated_renamed_argument('clobber', 'overwrite', '1.3')
        def test2(overwrite):
            return overwrite

        @deprecated_renamed_argument('clobber', 'overwrite', '1.3')
        def test3(self, overwrite):
            return overwrite

    @deprecated_renamed_argument('clobber', 'overwrite', '1.3', relax=False)
    def test1(overwrite):
        return overwrite

    for method in [Test().test1, Test().test2, Test().test3, test1]:
        # As positional argument only
        assert method(1) == 1

        # As new keyword argument
        assert method(overwrite=1) == 1

        # Using the deprecated name
        with catch_warnings(AstropyDeprecationWarning) as w:
            assert method(clobber=1) == 1
            assert len(w) == 1
            assert '1.3' in str(w[0].message)

        # Using both. Both keyword
        with pytest.raises(TypeError):
            method(clobber=2, overwrite=1)
        # One positional, one keyword
        with pytest.raises(TypeError):
            method(1, clobber=2)


def test_deprecated_argument_in_kwargs():
    # To rename an argument that is consumed by "kwargs" the "arg_in_kwargs"
    # parameter is used.
    @deprecated_renamed_argument('clobber', 'overwrite', '1.3',
                                 arg_in_kwargs=True)
    def test(**kwargs):
        return kwargs['overwrite']

    # As positional argument only
    with pytest.raises(TypeError):
        test(1)

    # As new keyword argument
    assert test(overwrite=1) == 1

    # Using the deprecated name
    with catch_warnings(AstropyDeprecationWarning) as w:
        assert test(clobber=1) == 1
        assert len(w) == 1
        assert '1.3' in str(w[0].message)

    # Using both. Both keyword
    with pytest.raises(TypeError):
        test(clobber=2, overwrite=1)
    # One positional, one keyword
    with pytest.raises(TypeError):
        test(1, clobber=2)


def test_deprecated_argument_relaxed():
    # Relax turns the TypeError if both old and new keyword are used into
    # a warning.
    @deprecated_renamed_argument('clobber', 'overwrite', '1.3', relax=True)
    def test(overwrite):
        return overwrite

    # As positional argument only
    assert test(1) == 1

    # As new keyword argument
    assert test(overwrite=1) == 1

    # Using the deprecated name
    with catch_warnings(AstropyDeprecationWarning) as w:
        assert test(clobber=1) == 1
        assert len(w) == 1
        assert '1.3' in str(w[0].message)

    # Using both. Both keyword
    with catch_warnings(AstropyUserWarning) as w:
        assert test(clobber=2, overwrite=1) == 1
        assert len(w) == 1

    # One positional, one keyword
    with catch_warnings(AstropyUserWarning) as w:
        assert test(1, clobber=2) == 1
        assert len(w) == 1


def test_deprecated_argument_pending():
    # Relax turns the TypeError if both old and new keyword are used into
    # a warning.
    @deprecated_renamed_argument('clobber', 'overwrite', '1.3', pending=True)
    def test(overwrite):
        return overwrite

    # As positional argument only
    assert test(1) == 1

    # As new keyword argument
    assert test(overwrite=1) == 1

    # Using the deprecated name
    with catch_warnings(AstropyUserWarning, AstropyDeprecationWarning) as w:
        assert test(clobber=1) == 1
        assert len(w) == 0

    # Using both. Both keyword
    with catch_warnings(AstropyUserWarning, AstropyDeprecationWarning) as w:
        assert test(clobber=2, overwrite=1) == 1
        assert len(w) == 0

    # One positional, one keyword
    with catch_warnings(AstropyUserWarning, AstropyDeprecationWarning) as w:
        assert test(1, clobber=2) == 1
        assert len(w) == 0


def test_deprecated_argument_multi_deprecation():
    @deprecated_renamed_argument(['x', 'y', 'z'], ['a', 'b', 'c'],
                                 [1.3, 1.2, 1.3], relax=True)
    def test(a, b, c):
        return a, b, c

    with catch_warnings(AstropyDeprecationWarning) as w:
        assert test(x=1, y=2, z=3) == (1, 2, 3)
        assert len(w) == 3

    # Make sure relax is valid for all arguments
    with catch_warnings(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, b=3) == (1, 3, 3)
        assert len(w) == 1

    with catch_warnings(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, a=3) == (3, 2, 3)
        assert len(w) == 1

    with catch_warnings(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, c=5) == (1, 2, 5)
        assert len(w) == 1


def test_deprecated_argument_multi_deprecation_2():
    @deprecated_renamed_argument(['x', 'y', 'z'], ['a', 'b', 'c'],
                                 [1.3, 1.2, 1.3], relax=[True, True, False])
    def test(a, b, c):
        return a, b, c

    with catch_warnings(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, b=3) == (1, 3, 3)
        assert len(w) == 1

    with catch_warnings(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, a=3) == (3, 2, 3)
        assert len(w) == 1

    with pytest.raises(TypeError):
        assert test(x=1, y=2, z=3, c=5) == (1, 2, 5)


def test_deprecated_argument_not_allowed_use():
    # If the argument is supposed to be inside the kwargs one needs to set the
    # arg_in_kwargs parameter. Without it it raises a TypeError.
    with pytest.raises(TypeError):
        @deprecated_renamed_argument('clobber', 'overwrite', '1.3')
        def test1(**kwargs):
            return kwargs['overwrite']

    # Cannot replace "*args".
    with pytest.raises(TypeError):
        @deprecated_renamed_argument('overwrite', 'args', '1.3')
        def test2(*args):
            return args

    # Cannot replace "**kwargs".
    with pytest.raises(TypeError):
        @deprecated_renamed_argument('overwrite', 'kwargs', '1.3')
        def test3(**kwargs):
            return kwargs


@pytest.mark.skipif('not six.PY2')
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


def test_format_doc_stringInput_simple():
    # Simple tests with string input

    docstring_fail = ''

    # Raises an valueerror if input is empty
    with pytest.raises(ValueError):
        @format_doc(docstring_fail)
        def testfunc_fail():
            pass

    docstring = 'test'

    # A first test that replaces an empty docstring
    @format_doc(docstring)
    def testfunc_1():
        pass
    assert inspect.getdoc(testfunc_1) == docstring

    # Test that it replaces an existing docstring
    @format_doc(docstring)
    def testfunc_2():
        '''not test'''
        pass
    assert inspect.getdoc(testfunc_2) == docstring


def test_format_doc_stringInput_format():
    # Tests with string input and formatting

    docstring = 'yes {0} no {opt}'

    # Raises an indexerror if not given the formatted args and kwargs
    with pytest.raises(IndexError):
        @format_doc(docstring)
        def testfunc1():
            pass

    # Test that the formatting is done right
    @format_doc(docstring, '/', opt='= life')
    def testfunc2():
        pass
    assert inspect.getdoc(testfunc2) == 'yes / no = life'

    # Test that we can include the original docstring

    docstring2 = 'yes {0} no {__doc__}'

    @format_doc(docstring2, '/')
    def testfunc3():
        '''= 2 / 2 * life'''
        pass
    assert inspect.getdoc(testfunc3) == 'yes / no = 2 / 2 * life'


def test_format_doc_objectInput_simple():
    # Simple tests with object input

    def docstring_fail():
        pass

    # Self input while the function has no docstring raises an error
    with pytest.raises(ValueError):
        @format_doc(docstring_fail)
        def testfunc_fail():
            pass

    def docstring0():
        '''test'''
        pass

    # A first test that replaces an empty docstring
    @format_doc(docstring0)
    def testfunc_1():
        pass
    assert inspect.getdoc(testfunc_1) == inspect.getdoc(docstring0)

    # Test that it replaces an existing docstring
    @format_doc(docstring0)
    def testfunc_2():
        '''not test'''
        pass
    assert inspect.getdoc(testfunc_2) == inspect.getdoc(docstring0)


def test_format_doc_objectInput_format():
    # Tests with object input and formatting

    def docstring():
        '''test {0} test {opt}'''
        pass

    # Raises an indexerror if not given the formatted args and kwargs
    with pytest.raises(IndexError):
        @format_doc(docstring)
        def testfunc_fail():
            pass

    # Test that the formatting is done right
    @format_doc(docstring, '+', opt='= 2 * test')
    def testfunc2():
        pass
    assert inspect.getdoc(testfunc2) == 'test + test = 2 * test'

    # Test that we can include the original docstring

    def docstring2():
        '''test {0} test {__doc__}'''
        pass

    @format_doc(docstring2, '+')
    def testfunc3():
        '''= 4 / 2 * test'''
        pass
    assert inspect.getdoc(testfunc3) == 'test + test = 4 / 2 * test'


def test_format_doc_selfInput_simple():
    # Simple tests with self input

    # Self input while the function has no docstring raises an error
    with pytest.raises(ValueError):
        @format_doc(None)
        def testfunc_fail():
            pass

    # Test that it keeps an existing docstring
    @format_doc(None)
    def testfunc_1():
        '''not test'''
        pass
    assert inspect.getdoc(testfunc_1) == 'not test'


def test_format_doc_selfInput_format():
    # Tests with string input which is '__doc__' (special case) and formatting

    # Raises an indexerror if not given the formatted args and kwargs
    with pytest.raises(IndexError):
        @format_doc(None)
        def testfunc_fail():
            '''dum {0} dum {opt}'''
            pass

    # Test that the formatting is done right
    @format_doc(None, 'di', opt='da dum')
    def testfunc1():
        '''dum {0} dum {opt}'''
        pass
    assert inspect.getdoc(testfunc1) == 'dum di dum da dum'

    # Test that we cannot recursively insert the original documentation

    @format_doc(None, 'di')
    def testfunc2():
        '''dum {0} dum {__doc__}'''
        pass
    assert inspect.getdoc(testfunc2) == 'dum di dum '


def test_format_doc_onMethod():
    # Check if the decorator works on methods too, to spice it up we try double
    # decorator
    docstring = 'what we do {__doc__}'

    class TestClass(object):
        @format_doc(docstring)
        @format_doc(None, 'strange.')
        def test_method(self):
            '''is {0}'''
            pass

    assert inspect.getdoc(TestClass.test_method) == 'what we do is strange.'


#@pytest.mark.skipif('six.PY2')
def test_format_doc_onClass():
    # Check if the decorator works on classes too
    docstring = 'what we do {__doc__} {0}{opt}'

    @format_doc(docstring, 'strange', opt='.')
    class TestClass(object):
        '''is'''
        pass

    assert inspect.getdoc(TestClass) == 'what we do is strange.'

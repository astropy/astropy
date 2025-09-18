# Licensed under a 3-clause BSD style license - see LICENSE.rst

import concurrent.futures
import inspect
import pickle
import sys

import pytest

from astropy.tests.helper import _skip_docstring_tests_with_optimized_python
from astropy.utils.decorators import (
    classproperty,
    deprecated,
    deprecated_attribute,
    deprecated_renamed_argument,
    format_doc,
    future_keyword_only,
    lazyproperty,
    sharedmethod,
)
from astropy.utils.exceptions import (
    AstropyDeprecationWarning,
    AstropyPendingDeprecationWarning,
    AstropyUserWarning,
)


class NewDeprecationWarning(AstropyDeprecationWarning):
    """
    New Warning subclass to be used to test the deprecated decorator's
    ``warning_type`` parameter.
    """


class NewPendingDeprecationWarning(AstropyPendingDeprecationWarning):
    """
    New Warning subclass to be used to test the deprecated decorator's
    ``pending_warning_type`` parameter.
    """


def test_deprecated_attribute():
    class DummyClass:
        def __init__(self):
            self.other = [42]
            self._foo = 42
            self._bar = 4242
            self._message = "42"
            self._pending = {42}
            self._pending_custom = {42}

        foo = deprecated_attribute("foo", "0.2")

        bar = deprecated_attribute("bar", "0.2", warning_type=NewDeprecationWarning)

        alternative = deprecated_attribute("alternative", "0.2", alternative="other")

        message = deprecated_attribute("message", "0.2", message="MSG")

        pending = deprecated_attribute("pending", "0.2", pending=True)

        pending_custom = deprecated_attribute(
            "pending_custom",
            "0.2",
            pending=True,
            pending_warning_type=NewPendingDeprecationWarning,
        )

    dummy = DummyClass()

    default_msg = (
        r"^The {} attribute is deprecated and may be removed in a future version\.$"
    )

    # Test getters and setters.
    msg = default_msg.format("foo")
    with pytest.warns(AstropyDeprecationWarning, match=msg) as w:
        assert dummy.foo == 42
    assert len(w) == 1
    with pytest.warns(AstropyDeprecationWarning, match=msg):
        dummy.foo = 24
    # Handling ``_foo`` should not cause deprecation warnings.
    assert dummy._foo == 24
    dummy._foo = 13
    assert dummy._foo == 13

    msg = default_msg.format("bar")
    with pytest.warns(NewDeprecationWarning, match=msg) as w:
        assert dummy.bar == 4242
    assert len(w) == 1
    with pytest.warns(NewDeprecationWarning, match=msg):
        dummy.bar = 2424

    with pytest.warns(AstropyDeprecationWarning, match="^MSG$"):
        assert dummy.message == "42"
    with pytest.warns(AstropyDeprecationWarning, match="^MSG$"):
        dummy.message = "24"

    msg = default_msg.format("alternative")[:-1] + r"\n        Use other instead\.$"
    with pytest.warns(AstropyDeprecationWarning, match=msg):
        assert dummy.alternative == [42]
    with pytest.warns(AstropyDeprecationWarning, match=msg):
        dummy.alternative = [24]
    # ``other`` is not deprecated.
    assert dummy.other == [24]
    dummy.other = [31]

    msg = r"^The pending attribute will be deprecated in a future version\.$"
    with pytest.warns(AstropyPendingDeprecationWarning, match=msg):
        assert dummy.pending == {42}
    with pytest.warns(AstropyPendingDeprecationWarning, match=msg):
        dummy.pending = {24}

    msg = r"^The pending_custom attribute will be deprecated in a future version\.$"
    with pytest.warns(NewPendingDeprecationWarning, match=msg):
        assert dummy.pending_custom == {42}
    with pytest.warns(NewPendingDeprecationWarning, match=msg):
        dummy.pending_custom = {24}


# This needs to be defined outside of the test function, because we
# want to try to pickle it.
@deprecated("100.0")
class TA:
    """
    This is the class docstring.
    """

    def __init__(self):
        """
        This is the __init__ docstring
        """


class TMeta(type):
    metaclass_attr = 1


@deprecated("100.0")
class TB(metaclass=TMeta):
    pass


@deprecated("100.0", warning_type=NewDeprecationWarning)
class TC:
    """
    This class has the custom warning.
    """


def test_deprecated_class():
    orig_A = TA.__bases__[0]

    # The only thing that should be different about the new class
    # is __doc__, __init__, __bases__, __subclasshook__,
    # and __init_subclass__.
    for x in dir(orig_A):
        if x not in (
            "__doc__",
            "__init__",
            "__bases__",
            "__dict__",
            "__subclasshook__",
            "__init_subclass__",
        ):
            assert getattr(TA, x) == getattr(orig_A, x)

    with pytest.warns(AstropyDeprecationWarning) as w:
        TA()
    assert len(w) == 1
    if TA.__doc__ is not None:
        assert "function" not in TA.__doc__
        assert "deprecated" in TA.__doc__
        assert "function" not in TA.__init__.__doc__
        assert "deprecated" in TA.__init__.__doc__

    # Test that the ``__deprecated__`` attribute is set
    # See https://peps.python.org/pep-0702/ for more information
    assert (
        TA.__deprecated__
        == "The TA class is deprecated and may be removed in a future version."
    )

    # Make sure the object is picklable
    pickle.dumps(TA)

    with pytest.warns(NewDeprecationWarning) as w:
        TC()
    assert len(w) == 1


def test_deprecated_class_with_new_method():
    """
    Test that a class with __new__ method still works even if it accepts
    additional arguments.
    This previously failed because the deprecated decorator would wrap objects
    __init__ which takes no arguments.
    """

    @deprecated("1.0")
    class A:
        def __new__(cls, a):
            return super().__new__(cls)

    # Creating an instance should work but raise a DeprecationWarning
    with pytest.warns(AstropyDeprecationWarning) as w:
        A(1)
    assert len(w) == 1

    @deprecated("1.0")
    class B:
        def __new__(cls, a):
            return super().__new__(cls)

        def __init__(self, a):
            pass

    # Creating an instance should work but raise a DeprecationWarning
    with pytest.warns(AstropyDeprecationWarning) as w:
        B(1)
    assert len(w) == 1


def test_deprecated_class_with_super():
    """
    Regression test for an issue where classes that used ``super()`` in their
    ``__init__`` did not actually call the correct class's ``__init__`` in the
    MRO.
    """

    @deprecated("100.0")
    class TB:
        def __init__(self, a, b):
            super().__init__()

    with pytest.warns(AstropyDeprecationWarning) as w:
        TB(1, 2)
    assert len(w) == 1
    if TB.__doc__ is not None:
        assert "function" not in TB.__doc__
        assert "deprecated" in TB.__doc__
        assert "function" not in TB.__init__.__doc__
        assert "deprecated" in TB.__init__.__doc__


def test_deprecated_class_with_custom_metaclass():
    """
    Regression test for an issue where deprecating a class with a metaclass
    other than type did not restore the metaclass properly.
    """

    with pytest.warns(AstropyDeprecationWarning) as w:
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

    class A:
        """Docstring"""

        @deprecated("1.0")
        @staticmethod
        def B():
            pass

        @deprecated("1.0")
        @classmethod
        def C(cls):
            pass

    with pytest.warns(AstropyDeprecationWarning) as w:
        A.B()
    assert len(w) == 1
    if A.__doc__ is not None:
        assert "deprecated" in A.B.__doc__

    # Test that the ``__deprecated__`` attribute is set
    # See https://peps.python.org/pep-0702/ for more information
    assert (
        A.B.__deprecated__
        == "The B method is deprecated and may be removed in a future version."
    )
    # And that it is not set on the original function which doesn't have
    # the deprecation warning.
    assert not hasattr(A.B.__wrapped__, "__deprecated__")

    with pytest.warns(AstropyDeprecationWarning) as w:
        A.C()
    assert len(w) == 1
    if A.__doc__ is not None:
        assert "deprecated" in A.C.__doc__

    # Test that the ``__deprecated__`` attribute is set
    # See https://peps.python.org/pep-0702/ for more information
    assert (
        A.C.__deprecated__
        == "The C method is deprecated and may be removed in a future version."
    )
    # And that it is not set on the original function which doesn't have
    # the deprecation warning.
    assert not hasattr(A.C.__wrapped__, "__deprecated__")


def test_deprecated_argument():
    # Tests the decorator with function, method, staticmethod and classmethod.

    class Test:
        @classmethod
        @deprecated_renamed_argument("clobber", "overwrite", "1.3")
        def test1(cls, overwrite):
            return overwrite

        @staticmethod
        @deprecated_renamed_argument("clobber", "overwrite", "1.3")
        def test2(overwrite):
            return overwrite

        @deprecated_renamed_argument("clobber", "overwrite", "1.3")
        def test3(self, overwrite):
            return overwrite

        @deprecated_renamed_argument(
            "clobber", "overwrite", "1.3", warning_type=NewDeprecationWarning
        )
        def test4(self, overwrite):
            return overwrite

    @deprecated_renamed_argument("clobber", "overwrite", "1.3", relax=False)
    def test1(overwrite):
        return overwrite

    for method in [Test().test1, Test().test2, Test().test3, Test().test4, test1]:
        # As positional argument only
        assert method(1) == 1

        # As new keyword argument
        assert method(overwrite=1) == 1

        # Using the deprecated name
        with pytest.warns(AstropyDeprecationWarning, match=r"1\.3") as w:
            assert method(clobber=1) == 1
        assert len(w) == 1
        assert "test_decorators.py" in str(w[0].filename)
        if method.__name__ == "test4":
            assert issubclass(w[0].category, NewDeprecationWarning)

        # Using both. Both keyword
        with pytest.raises(TypeError), pytest.warns(AstropyDeprecationWarning):
            method(clobber=2, overwrite=1)
        # One positional, one keyword
        with pytest.raises(TypeError), pytest.warns(AstropyDeprecationWarning):
            method(1, clobber=2)


def test_deprecated_argument_custom_message():
    @deprecated_renamed_argument("foo", "bar", "4.0", message="Custom msg")
    def test(bar=0):
        pass

    with pytest.warns(AstropyDeprecationWarning, match="Custom msg"):
        test(foo=0)


def test_deprecated_argument_in_kwargs():
    # To rename an argument that is consumed by "kwargs" the "arg_in_kwargs"
    # parameter is used.
    @deprecated_renamed_argument("clobber", "overwrite", "1.3", arg_in_kwargs=True)
    def test(**kwargs):
        return kwargs["overwrite"]

    # As positional argument only
    with pytest.raises(TypeError):
        test(1)

    # As new keyword argument
    assert test(overwrite=1) == 1

    # Using the deprecated name
    with pytest.warns(AstropyDeprecationWarning, match=r"1\.3") as w:
        assert test(clobber=1) == 1
    assert len(w) == 1
    assert "test_decorators.py" in str(w[0].filename)

    # Using both. Both keyword
    with pytest.raises(TypeError), pytest.warns(AstropyDeprecationWarning):
        test(clobber=2, overwrite=1)
    # One positional, one keyword
    with pytest.raises(TypeError), pytest.warns(AstropyDeprecationWarning):
        test(1, clobber=2)


def test_deprecated_argument_relaxed():
    # Relax turns the TypeError if both old and new keyword are used into
    # a warning.
    @deprecated_renamed_argument("clobber", "overwrite", "1.3", relax=True)
    def test(overwrite):
        return overwrite

    # As positional argument only
    assert test(1) == 1

    # As new keyword argument
    assert test(overwrite=1) == 1

    # Using the deprecated name
    with pytest.warns(AstropyDeprecationWarning, match=r"1\.3") as w:
        assert test(clobber=1) == 1
    assert len(w) == 1

    # Using both. Both keyword
    with pytest.warns(AstropyDeprecationWarning), pytest.warns(AstropyUserWarning) as w:
        assert test(clobber=2, overwrite=1) == 1
    assert len(w) == 2
    assert '"clobber" was deprecated' in str(w[0].message)
    assert '"clobber" and "overwrite" keywords were set' in str(w[1].message)

    # One positional, one keyword
    with pytest.warns(AstropyDeprecationWarning), pytest.warns(AstropyUserWarning) as w:
        assert test(1, clobber=2) == 1
    assert len(w) == 2
    assert '"clobber" was deprecated' in str(w[0].message)
    assert '"clobber" and "overwrite" keywords were set' in str(w[1].message)


def test_deprecated_argument_pending():
    # Relax turns the TypeError if both old and new keyword are used into
    # a warning.
    @deprecated_renamed_argument("clobber", "overwrite", "1.3", pending=True)
    def test(overwrite):
        return overwrite

    # As positional argument only
    assert test(1) == 1

    # As new keyword argument
    assert test(overwrite=1) == 1

    # Using the deprecated name
    assert test(clobber=1) == 1

    # Using both. Both keyword
    assert test(clobber=2, overwrite=1) == 1

    # One positional, one keyword
    assert test(1, clobber=2) == 1


def test_deprecated_argument_multi_deprecation():
    @deprecated_renamed_argument(
        ["x", "y", "z"], ["a", "b", "c"], [1.3, 1.2, 1.3], relax=True
    )
    def test(a, b, c):
        return a, b, c

    with pytest.warns(AstropyDeprecationWarning) as w:
        assert test(x=1, y=2, z=3) == (1, 2, 3)
    assert len(w) == 3

    # Make sure relax is valid for all arguments
    with pytest.warns(AstropyDeprecationWarning), pytest.warns(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, b=3) == (1, 3, 3)
    assert len(w) == 4

    with pytest.warns(AstropyDeprecationWarning), pytest.warns(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, a=3) == (3, 2, 3)
    assert len(w) == 4

    with pytest.warns(AstropyDeprecationWarning), pytest.warns(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, c=5) == (1, 2, 5)
    assert len(w) == 4


def test_deprecated_argument_multi_deprecation_2():
    @deprecated_renamed_argument(
        ["x", "y", "z"], ["a", "b", "c"], [1.3, 1.2, 1.3], relax=[True, True, False]
    )
    def test(a, b, c):
        return a, b, c

    with pytest.warns(AstropyDeprecationWarning), pytest.warns(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, b=3) == (1, 3, 3)
    assert len(w) == 4

    with pytest.warns(AstropyDeprecationWarning), pytest.warns(AstropyUserWarning) as w:
        assert test(x=1, y=2, z=3, a=3) == (3, 2, 3)
    assert len(w) == 4

    with pytest.raises(TypeError), pytest.warns(AstropyDeprecationWarning):
        assert test(x=1, y=2, z=3, c=5) == (1, 2, 5)


def test_deprecated_argument_not_allowed_use():
    # If the argument is supposed to be inside the kwargs one needs to set the
    # arg_in_kwargs parameter. Without it it raises a TypeError.
    with pytest.raises(TypeError):

        @deprecated_renamed_argument("clobber", "overwrite", "1.3")
        def test1(**kwargs):
            return kwargs["overwrite"]

    # Cannot replace "*args".
    with pytest.raises(TypeError):

        @deprecated_renamed_argument("overwrite", "args", "1.3")
        def test2(*args):
            return args

    # Cannot replace "**kwargs".
    with pytest.raises(TypeError):

        @deprecated_renamed_argument("overwrite", "kwargs", "1.3")
        def test3(**kwargs):
            return kwargs


def test_deprecated_argument_remove():
    @deprecated_renamed_argument("x", None, "2.0", alternative="astropy.y")
    def test(dummy=11, x=3):
        return dummy, x

    with pytest.warns(AstropyDeprecationWarning, match=r"Use astropy\.y instead") as w:
        assert test(x=1) == (11, 1)
    assert len(w) == 1

    with pytest.warns(AstropyDeprecationWarning) as w:
        assert test(x=1, dummy=10) == (10, 1)
    assert len(w) == 1

    with pytest.warns(AstropyDeprecationWarning, match=r"Use astropy\.y instead"):
        test(121, 1) == (121, 1)  # noqa: B015

    assert test() == (11, 3)
    assert test(121) == (121, 3)
    assert test(dummy=121) == (121, 3)


def test_future_kwo_invalid_since():
    with pytest.raises(
        ValueError,
        match=(
            "Expected names and since values with "
            r"identical length. Got len\(names\)=3 "
            r"and len\(since\)=2"
        ),
    ):
        future_keyword_only(["a", "b", "c"], since=["7.1", "7.1"])(lambda a, b, c: ...)


def test_future_kwo_invalid_names():
    with pytest.raises(
        ValueError,
        match=(
            "The following arguments cannot be marked as future keyword-only "
            "because they were not found in the decorated function's signature: a"
        ),
    ):
        future_keyword_only(["a"], since=["7.1"])(lambda b: ...)


@pytest.mark.parametrize(
    "input_func",
    [
        pytest.param(lambda a, /: ..., id="arg-is-positional-only"),
        pytest.param(lambda *, a: ..., id="arg-is-keyword-only"),
        pytest.param(lambda *a: ..., id="arg-is-var-positional"),
        pytest.param(lambda **a: ..., id="arg-is-var-keyword"),
    ],
)
def test_future_kwo_not_pos_or_kw(input_func):
    with pytest.raises(
        ValueError,
        match=(
            "The following arguments cannot be marked as future keyword-only "
            "because they are not currently positional-or-keyword: a"
        ),
    ):
        future_keyword_only(["a"], since=["7.1"])(input_func)


def test_future_kwo_broken_pos_allowed():
    with pytest.raises(
        ValueError,
        match=(
            "The following positionally-allowed arguments were not marked as "
            "future keyword-only and would be broken under future keyword-only "
            r"requirements: c, d"
        ),
    ):
        future_keyword_only(["b"], since=["7.1"])(lambda a, b, c, d=None: ...)


def test_future_kwo_warn_if_and_only_if_needed():
    func = future_keyword_only(["b", "c", "d"], since=["7.1", "7.1", "7.2"])(
        lambda a, b, c=None, d=None: ...
    )
    with pytest.warns(
        AstropyDeprecationWarning,
        match=(
            "The following arguments were received positionally, which "
            "will be disallowed in a future release: b, c, d\n"
            "Pass them as keywords to suppress this warning. "
            r"\(deprecated since 7\.1, 7\.1, 7\.2, respectively\)"
        ),
    ):
        func(1, 2, 3, 4)

    # also check the warning if all arguments were marked at the same time
    with pytest.warns(
        AstropyDeprecationWarning,
        match=(
            "The following arguments were received positionally, which "
            "will be disallowed in a future release: b, c\n"
            "Pass them as keywords to suppress this warning. "
            r"\(deprecated since 7\.1\)"
        ),
    ):
        func(1, 2, 3, d=4)

    # check no warning is emitted if future constraints are already honored
    func(1, b=2, c=3, d=4)
    func(1, b=2, c=3)
    func(1, b=2)


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

    class A:
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

    class A:
        # Inherits docstring from getter
        @classproperty
        def foo(cls):
            """The foo."""

            return 1

    expected_doc = "The foo." if sys.flags.optimize < 2 else None
    assert inspect.getdoc(A.__dict__["foo"]) == expected_doc

    class B:
        # Use doc passed to classproperty constructor
        def _get_foo(cls):
            return 1

        foo = classproperty(_get_foo, doc="The foo.")

    # we should *always* get a string back by setting the doc argument.
    # As of Python 3.13, this is in line with how the builtin @property decorator
    # interacts with PYTHONOPTIMIZE=2
    assert inspect.getdoc(B.__dict__["foo"]) == "The foo."


@pytest.mark.slow
def test_classproperty_lazy_threadsafe(fast_thread_switching):
    """
    Test that a class property with lazy=True is thread-safe.
    """
    workers = 8
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        # This is testing for race conditions, so try many times in the
        # hope that we'll get the timing right.
        for _ in range(10000):

            class A:
                @classproperty(lazy=True)
                def foo(cls):
                    nonlocal calls
                    calls += 1
                    return object()

            # Have all worker threads query in parallel
            calls = 0
            futures = [executor.submit(lambda: A.foo) for i in range(workers)]
            # Check that only one call happened and they all received it
            values = [future.result() for future in futures]
            assert calls == 1
            assert values[0] is not None
            assert values == [values[0]] * workers


@pytest.mark.slow
def test_lazyproperty_threadsafe(fast_thread_switching):
    """
    Test thread safety of lazyproperty.
    """
    # This test is generally similar to test_classproperty_lazy_threadsafe
    # above. See there for comments.

    class A:
        def __init__(self):
            self.calls = 0

        @lazyproperty
        def foo(self):
            self.calls += 1
            return object()

    workers = 8
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        for _ in range(10000):
            a = A()
            futures = [executor.submit(lambda: a.foo) for i in range(workers)]
            values = [future.result() for future in futures]
            assert a.calls == 1
            assert a.foo is not None
            assert values == [a.foo] * workers


@_skip_docstring_tests_with_optimized_python
def test_format_doc_stringInput_simple():
    # Simple tests with string input
    docstring = "test"

    # A first test that replaces an empty docstring
    @format_doc(docstring)
    def testfunc_1():
        pass

    assert inspect.getdoc(testfunc_1) == "test"

    # Test that it replaces an existing docstring
    @format_doc(docstring)
    def testfunc_2():
        """not test"""

    assert inspect.getdoc(testfunc_2) == "test"


@_skip_docstring_tests_with_optimized_python
def test_format_doc_stringInput_format():
    # Tests with string input and formatting

    docstring = "yes {0} no {opt}"

    # Test that the formatting is done right
    @format_doc(docstring, "/", opt="= life")
    def testfunc2():
        pass

    assert inspect.getdoc(testfunc2) == "yes / no = life"

    # Test that we can include the original docstring

    docstring2 = "yes {0} no {__doc__}"

    @format_doc(docstring2, "/")
    def testfunc3():
        """= 2 / 2 * life"""

    assert inspect.getdoc(testfunc3) == "yes / no = 2 / 2 * life"


def test_format_doc_objectInput_simple():
    # Simple tests with object input

    def docstring0():
        """test"""

    # A first test that replaces an empty docstring
    @format_doc(docstring0)
    def testfunc_1():
        pass

    assert inspect.getdoc(testfunc_1) == inspect.getdoc(docstring0)

    # Test that it replaces an existing docstring
    @format_doc(docstring0)
    def testfunc_2():
        """not test"""

    assert inspect.getdoc(testfunc_2) == inspect.getdoc(docstring0)


@_skip_docstring_tests_with_optimized_python
def test_format_doc_objectInput_format():
    # Tests with object input and formatting

    def docstring():
        """test {0} test {opt}"""

    # Test that the formatting is done right
    @format_doc(docstring, "+", opt="= 2 * test")
    def testfunc2():
        pass

    assert inspect.getdoc(testfunc2) == "test + test = 2 * test"

    # Test that we can include the original docstring

    def docstring2():
        """test {0} test {__doc__}"""

    @format_doc(docstring2, "+")
    def testfunc3():
        """= 4 / 2 * test"""

    assert inspect.getdoc(testfunc3) == "test + test = 4 / 2 * test"


@_skip_docstring_tests_with_optimized_python
def test_format_doc_selfInput_simple():
    # Simple tests with self input

    # Test that it keeps an existing docstring
    @format_doc(None)
    def testfunc_1():
        """not test"""

    assert inspect.getdoc(testfunc_1) == "not test"


@_skip_docstring_tests_with_optimized_python
def test_format_doc_selfInput_format():
    # Tests with string input which is '__doc__' (special case) and formatting

    # Test that the formatting is done right
    @format_doc(None, "di", opt="da dum")
    def testfunc1():
        """dum {0} dum {opt}"""

    assert inspect.getdoc(testfunc1) == "dum di dum da dum"

    # Test that we cannot recursively insert the original documentation

    @format_doc(None, "di")
    def testfunc2():
        """dum {0} dum {__doc__}"""

    assert inspect.getdoc(testfunc2) == "dum di dum "


@_skip_docstring_tests_with_optimized_python
def test_format_doc_onMethod():
    # Check if the decorator works on methods too, to spice it up we try double
    # decorator
    docstring = "what we do {__doc__}"

    class TestClass:
        @format_doc(docstring)
        @format_doc(None, "strange.")
        def test_method(self):
            """is {0}"""

    assert inspect.getdoc(TestClass.test_method) == "what we do is strange."


@_skip_docstring_tests_with_optimized_python
def test_format_doc_onClass():
    # Check if the decorator works on classes too
    docstring = "what we do {__doc__} {0}{opt}"

    @format_doc(docstring, "strange", opt=".")
    class TestClass:
        """is"""

    assert inspect.getdoc(TestClass) == "what we do is strange."


@_skip_docstring_tests_with_optimized_python
@pytest.mark.parametrize(
    "docstring, expected_exception",
    [
        # Raises an valueerror if input is empty
        pytest.param("", ValueError, id="empty string"),
        # Raises an indexerror if not given the formatted args and kwargs
        pytest.param("yes {0} no {opt}", IndexError, id="missing args or kwargs (str)"),
        # Self input while the function has no docstring raises an error
        pytest.param(lambda: None, ValueError, id="function without a docstring"),
        # Self input while the function has no docstring raises an error
        pytest.param(None, ValueError, id="None"),
    ],
)
def test_format_doc_exceptions(docstring, expected_exception):
    with pytest.raises(expected_exception):

        @format_doc(docstring)
        def testfunc_fail():
            pass


@_skip_docstring_tests_with_optimized_python
def test_format_doc_indexerrors():
    def _FUNC_WITH_TEMPLATE_DOCSTRING():
        """test {0} test {opt}"""

    # Raises an indexerror if not given the formatted args and kwargs
    with pytest.raises(IndexError):

        @format_doc(_FUNC_WITH_TEMPLATE_DOCSTRING)
        def testfunc_fail():
            pass

    with pytest.raises(IndexError):
        format_doc(None)(_FUNC_WITH_TEMPLATE_DOCSTRING)

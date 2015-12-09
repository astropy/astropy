# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides the tools used to internally run the astropy test suite
from the installed astropy.  It makes use of the `pytest` testing framework.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import base64
import errno
import functools
import os
import sys
import types
import warnings
import zlib

from ..extern import six
from ..extern.six.moves import cPickle as pickle

try:
    # Import pkg_resources to prevent it from issuing warnings upon being
    # imported from within py.test.  See
    # https://github.com/astropy/astropy/pull/537 for a detailed explanation.
    import pkg_resources
except ImportError:
    pass

from .. import test
from ..utils.exceptions import (AstropyWarning,
                                AstropyDeprecationWarning,
                                AstropyPendingDeprecationWarning)


# For backward-compatibility with affiliated packages
from .runner import TestRunner

__all__ = ['raises', 'enable_deprecations_as_exceptions', 'remote_data',
           'treat_deprecations_as_exceptions', 'catch_warnings',
           'assert_follows_unicode_guidelines', 'quantity_allclose',
           'assert_quantity_allclose', 'check_pickling_recovery',
           'pickle_protocol', 'generic_recursive_equality_test']


if os.environ.get('ASTROPY_USE_SYSTEM_PYTEST') or '_pytest' in sys.modules:
    import pytest

else:
    from ..extern import pytest as extern_pytest

    if six.PY3:
        exec("def do_exec_def(co, loc): exec(co, loc)\n")
        extern_pytest.do_exec = do_exec_def

        unpacked_sources = extern_pytest.sources.encode("ascii")
        unpacked_sources = pickle.loads(
            zlib.decompress(base64.decodebytes(unpacked_sources)), encoding='utf-8')
    elif six.PY2:
        exec("def do_exec_def(co, loc): exec co in loc\n")
        extern_pytest.do_exec = do_exec_def

        unpacked_sources = pickle.loads(
            zlib.decompress(base64.decodestring(extern_pytest.sources)))

    importer = extern_pytest.DictImporter(unpacked_sources)
    sys.meta_path.insert(0, importer)

    pytest = importer.load_module(str('pytest'))


# Monkey-patch py.test to work around issue #811
# https://github.com/astropy/astropy/issues/811
from _pytest.assertion import rewrite as _rewrite
_orig_write_pyc = _rewrite._write_pyc


def _write_pyc_wrapper(*args):
    """Wraps the internal _write_pyc method in py.test to recognize
    PermissionErrors and just stop trying to cache its generated pyc files if
    it can't write them to the __pycache__ directory.

    When py.test scans for test modules, it actually rewrites the bytecode
    of each test module it discovers--this is how it manages to add extra
    instrumentation to the assert builtin.  Normally it caches these
    rewritten bytecode files--``_write_pyc()`` is just a function that handles
    writing the rewritten pyc file to the cache.  If it returns ``False`` for
    any reason py.test will stop trying to cache the files altogether.  The
    original function catches some cases, but it has a long-standing bug of
    not catching permission errors on the ``__pycache__`` directory in Python
    3.  Hence this patch.
    """

    try:
        return _orig_write_pyc(*args)
    except IOError as e:
        if e.errno == errno.EACCES:
            return False
_rewrite._write_pyc = _write_pyc_wrapper


# pytest marker to mark tests which get data from the web
remote_data = pytest.mark.remote_data


# This is for Python 2.x and 3.x compatibility.  distutils expects
# options to all be byte strings on Python 2 and Unicode strings on
# Python 3.
def _fix_user_options(options):
    def to_str_or_none(x):
        if x is None:
            return None
        return str(x)

    return [tuple(to_str_or_none(x) for x in y) for y in options]


def _save_coverage(cov, result, rootdir, testing_path):
    """
    This method is called after the tests have been run in coverage mode
    to cleanup and then save the coverage data and report.
    """
    from ..utils.console import color_print

    if result != 0:
        return

    # The coverage report includes the full path to the temporary
    # directory, so we replace all the paths with the true source
    # path. This means that the coverage line-by-line report will only
    # be correct for Python 2 code (since the Python 3 code will be
    # different in the build directory from the source directory as
    # long as 2to3 is needed). Therefore we only do this fix for
    # Python 2.x.
    if six.PY2:
        try:
            # Coverage 4.0: _harvest_data has been renamed to get_data, the
            # lines dict is private
            cov.get_data()
        except AttributeError:
            # Coverage < 4.0
            cov._harvest_data()
            lines = cov.data.lines
        else:
            lines = cov.data._lines

        for key in lines.keys():
            new_path = os.path.relpath(
                os.path.realpath(key),
                os.path.realpath(testing_path))
            new_path = os.path.abspath(
                os.path.join(rootdir, new_path))
            lines[new_path] = lines.pop(key)

    color_print('Saving coverage data in .coverage...', 'green')
    cov.save()

    color_print('Saving HTML coverage report in htmlcov...', 'green')
    cov.html_report(directory=os.path.join(rootdir, 'htmlcov'))


class raises(object):
    """
    A decorator to mark that a test should raise a given exception.
    Use as follows::

        @raises(ZeroDivisionError)
        def test_foo():
            x = 1/0

    This can also be used a context manager, in which case it is just
    an alias for the ``pytest.raises`` context manager (because the
    two have the same name this help avoid confusion by being
    flexible).
    """

    # pep-8 naming exception -- this is a decorator class
    def __init__(self, exc):
        self._exc = exc
        self._ctx = None

    def __call__(self, func):
        @functools.wraps(func)
        def run_raises_test(*args, **kwargs):
            pytest.raises(self._exc, func, *args, **kwargs)
        return run_raises_test

    def __enter__(self):
        self._ctx = pytest.raises(self._exc)
        return self._ctx.__enter__()

    def __exit__(self, *exc_info):
        return self._ctx.__exit__(*exc_info)


_deprecations_as_exceptions = False
_include_astropy_deprecations = True

def enable_deprecations_as_exceptions(include_astropy_deprecations=True):
    """
    Turn on the feature that turns deprecations into exceptions.
    """
    global _deprecations_as_exceptions
    _deprecations_as_exceptions = True

    global _include_astropy_deprecations
    _include_astropy_deprecations = include_astropy_deprecations


def treat_deprecations_as_exceptions():
    """
    Turn all DeprecationWarnings (which indicate deprecated uses of
    Python itself or Numpy, but not within Astropy, where we use our
    own deprecation warning class) into exceptions so that we find
    out about them early.

    This completely resets the warning filters and any "already seen"
    warning state.
    """
    # First, totally reset the warning state
    for module in list(six.itervalues(sys.modules)):
        # We don't want to deal with six.MovedModules, only "real"
        # modules.
        if (isinstance(module, types.ModuleType) and
            hasattr(module, '__warningregistry__')):
            del module.__warningregistry__

    if not _deprecations_as_exceptions:
        return

    warnings.resetwarnings()

    # Hide the next couple of DeprecationWarnings
    warnings.simplefilter('ignore', DeprecationWarning)
    # Here's the wrinkle: a couple of our third-party dependencies
    # (py.test and scipy) are still using deprecated features
    # themselves, and we'd like to ignore those.  Fortunately, those
    # show up only at import time, so if we import those things *now*,
    # before we turn the warnings into exceptions, we're golden.
    try:
        # A deprecated stdlib module used by py.test
        import compiler
    except ImportError:
        pass

    try:
        import scipy
    except ImportError:
        pass

    # Now, start over again with the warning filters
    warnings.resetwarnings()
    # Now, turn DeprecationWarnings into exceptions
    warnings.filterwarnings("error", ".*", DeprecationWarning)

    # Only turn astropy deprecation warnings into exceptions if requested
    if _include_astropy_deprecations:
        warnings.filterwarnings("error", ".*", AstropyDeprecationWarning)
        warnings.filterwarnings("error", ".*", AstropyPendingDeprecationWarning)

    if sys.version_info[:2] == (2, 6):
        # py.test's warning.showwarning does not include the line argument
        # on Python 2.6, so we need to explicitly ignore this warning.
        warnings.filterwarnings(
            "ignore",
            r"functions overriding warnings\.showwarning\(\) must support "
            r"the 'line' argument",
            DeprecationWarning)

    if sys.version_info[:2] >= (3, 4):
        # py.test reads files with the 'U' flag, which is now
        # deprecated in Python 3.4.
        warnings.filterwarnings(
            "ignore",
            r"'U' mode is deprecated",
            DeprecationWarning)

        # BeautifulSoup4 triggers a DeprecationWarning in stdlib's
        # html module.x
        warnings.filterwarnings(
            "ignore",
            r"The strict argument and mode are deprecated\.",
            DeprecationWarning)
        warnings.filterwarnings(
            "ignore",
            r"The value of convert_charrefs will become True in 3\.5\. "
            r"You are encouraged to set the value explicitly\.",
            DeprecationWarning)

    if sys.version_info[:2] >= (3, 5):
        # py.test raises this warning on Python 3.5.
        # This can be removed when fixed in py.test.
        # See https://github.com/pytest-dev/pytest/pull/1009
        warnings.filterwarnings(
            "ignore",
            r"inspect\.getargspec\(\) is deprecated, use "
            r"inspect\.signature\(\) instead",
            DeprecationWarning)


class catch_warnings(warnings.catch_warnings):
    """
    A high-powered version of warnings.catch_warnings to use for testing
    and to make sure that there is no dependence on the order in which
    the tests are run.

    This completely blitzes any memory of any warnings that have
    appeared before so that all warnings will be caught and displayed.

    ``*args`` is a set of warning classes to collect.  If no arguments are
    provided, all warnings are collected.

    Use as follows::

        with catch_warnings(MyCustomWarning) as w:
            do.something.bad()
        assert len(w) > 0
    """
    def __init__(self, *classes):
        super(catch_warnings, self).__init__(record=True)
        self.classes = classes

    def __enter__(self):
        warning_list = super(catch_warnings, self).__enter__()
        treat_deprecations_as_exceptions()
        if len(self.classes) == 0:
            warnings.simplefilter('always')
        else:
            warnings.simplefilter('ignore')
            for cls in self.classes:
                warnings.simplefilter('always', cls)
        return warning_list

    def __exit__(self, type, value, traceback):
        treat_deprecations_as_exceptions()


class ignore_warnings(catch_warnings):
    """
    This can be used either as a context manager or function decorator to
    ignore all warnings that occur within a function or block of code.

    An optional category option can be supplied to only ignore warnings of a
    certain category or categories (if a list is provided).
    """

    def __init__(self, category=None):
        super(ignore_warnings, self).__init__()

        if isinstance(category, type) and issubclass(category, Warning):
            self.category = [category]
        else:
            self.category = category

    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Originally this just reused self, but that doesn't work if the
            # function is called more than once so we need to make a new
            # context manager instance for each call
            with self.__class__(category=self.category):
                return func(*args, **kwargs)

        return wrapper

    def __enter__(self):
        retval = super(ignore_warnings, self).__enter__()
        if self.category is not None:
            for category in self.category:
                warnings.simplefilter('ignore', category)
        else:
            warnings.simplefilter('ignore')
        return retval


def assert_follows_unicode_guidelines(
        x, roundtrip=None):
    """
    Test that an object follows our Unicode policy.  See
    "Unicode guidelines" in the coding guidelines.

    Parameters
    ----------
    x : object
        The instance to test

    roundtrip : module, optional
        When provided, this namespace will be used to evaluate
        ``repr(x)`` and ensure that it roundtrips.  It will also
        ensure that ``__bytes__(x)`` and ``__unicode__(x)`` roundtrip.
        If not provided, no roundtrip testing will be performed.
    """
    from .. import conf
    from ..extern import six

    with conf.set_temp('unicode_output', False):
        bytes_x = bytes(x)
        unicode_x = six.text_type(x)
        repr_x = repr(x)

        assert isinstance(bytes_x, bytes)
        bytes_x.decode('ascii')
        assert isinstance(unicode_x, six.text_type)
        unicode_x.encode('ascii')
        assert isinstance(repr_x, six.string_types)
        if isinstance(repr_x, bytes):
            repr_x.decode('ascii')
        else:
            repr_x.encode('ascii')

        if roundtrip is not None:
            assert x.__class__(bytes_x) == x
            assert x.__class__(unicode_x) == x
            assert eval(repr_x, roundtrip) == x

    with conf.set_temp('unicode_output', True):
        bytes_x = bytes(x)
        unicode_x = six.text_type(x)
        repr_x = repr(x)

        assert isinstance(bytes_x, bytes)
        bytes_x.decode('ascii')
        assert isinstance(unicode_x, six.text_type)
        assert isinstance(repr_x, six.string_types)
        if isinstance(repr_x, bytes):
            repr_x.decode('ascii')
        else:
            repr_x.encode('ascii')

        if roundtrip is not None:
            assert x.__class__(bytes_x) == x
            assert x.__class__(unicode_x) == x
            assert eval(repr_x, roundtrip) == x


@pytest.fixture(params=[0, 1, -1])
def pickle_protocol(request):
    """
    Fixture to run all the tests for protocols 0 and 1, and -1 (most advanced).
    (Originally from astropy.table.tests.test_pickle)
    """
    return request.param


def generic_recursive_equality_test(a, b, class_history):
    """
    Check if the attributes of a and b are equal. Then,
    check if the attributes of the attributes are equal.
    """
    dict_a = a.__dict__
    dict_b = b.__dict__
    for key in dict_a:
        assert key in dict_b,\
          "Did not pickle {0}".format(key)
        if hasattr(dict_a[key], '__eq__'):
            eq = (dict_a[key] == dict_b[key])
            if '__iter__' in dir(eq):
                eq = (False not in eq)
            assert eq, "Value of {0} changed by pickling".format(key)

        if hasattr(dict_a[key], '__dict__'):
            if dict_a[key].__class__ in class_history:
                #attempt to prevent infinite recursion
                pass
            else:
                new_class_history = [dict_a[key].__class__]
                new_class_history.extend(class_history)
                generic_recursive_equality_test(dict_a[key],
                                                dict_b[key],
                                                new_class_history)


def check_pickling_recovery(original, protocol):
    """
    Try to pickle an object. If successful, make sure
    the object's attributes survived pickling and unpickling.
    """
    f = pickle.dumps(original, protocol=protocol)
    unpickled = pickle.loads(f)
    class_history = [original.__class__]
    generic_recursive_equality_test(original, unpickled,
                                    class_history)


def assert_quantity_allclose(actual, desired, rtol=1.e-7, atol=None,
                             **kwargs):
    """
    Raise an assertion if two objects are not equal up to desired tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.testing.assert_allclose`.
    """
    import numpy as np
    np.testing.assert_allclose(*_unquantify_allclose_arguments(actual, desired,
                                                               rtol, atol),
                               **kwargs)


def quantity_allclose(a, b, rtol=1.e-5, atol=None, **kwargs):
    """
    Returns True if two arrays are element-wise equal within a tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.allclose`.
    """
    import numpy as np
    return np.allclose(*_unquantify_allclose_arguments(a, b, rtol, atol),
                       **kwargs)


def _unquantify_allclose_arguments(actual, desired, rtol, atol):
    from .. import units as u

    actual = u.Quantity(actual, subok=True, copy=False)

    desired = u.Quantity(desired, subok=True, copy=False)
    try:
        desired = desired.to(actual.unit)
    except u.UnitsError:
        raise u.UnitsError("Units for 'desired' ({0}) and 'actual' ({1}) "
                           "are not convertible"
                           .format(desired.unit, actual.unit))

    if atol is None:
        # by default, we assume an absolute tolerance of 0
        atol = u.Quantity(0)
    else:
        atol = u.Quantity(atol, subok=True, copy=False)
        try:
            atol = atol.to(actual.unit)
        except u.UnitsError:
            raise u.UnitsError("Units for 'atol' ({0}) and 'actual' ({1}) "
                               "are not convertible"
                               .format(atol.unit, actual.unit))

    rtol =  u.Quantity(rtol, subok=True, copy=False)
    try:
        rtol = rtol.to(u.dimensionless_unscaled)
    except:
        raise u.UnitsError("`rtol` should be dimensionless")

    return actual.value, desired.value, rtol.value, atol.value

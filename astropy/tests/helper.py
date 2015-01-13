# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module prvoides the tools used to internally run the astropy test suite
from the installed astropy.  It makes use of the `pytest` testing framework.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..extern import six
from ..extern.six.moves import cPickle as pickle

import errno
import shlex
import sys
import base64
import zlib
import functools
import multiprocessing
import os
import subprocess
import shutil
import tempfile
import types
import warnings

try:
    # Import pkg_resources to prevent it from issuing warnings upon being
    # imported from within py.test.  See
    # https://github.com/astropy/astropy/pull/537 for a detailed explanation.
    import pkg_resources
except ImportError:
    pass

from distutils.core import Command

from .. import test
from ..utils.exceptions import (AstropyWarning,
                                AstropyDeprecationWarning,
                                AstropyPendingDeprecationWarning)
from ..config import configuration

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


class TestRunner(object):
    def __init__(self, base_path):
        self.base_path = base_path

    def run_tests(self, package=None, test_path=None, args=None, plugins=None,
                  verbose=False, pastebin=None, remote_data=False, pep8=False,
                  pdb=False, coverage=False, open_files=False, parallel=0,
                  docs_path=None, skip_docs=False, repeat=None):
        """
        The docstring for this method lives in astropy/__init__.py:test
        """
        try:
            get_ipython()
        except NameError:
            pass
        else:
            raise RuntimeError(
                "Running astropy tests inside of IPython is not supported.")

        if coverage:
            warnings.warn(
                "The coverage option is ignored on run_tests, since it "
                "can not be made to work in that context.  Use "
                "'python setup.py test --coverage' instead.",
                AstropyWarning)

        all_args = []

        if package is None:
            package_path = self.base_path
        else:
            package_path = os.path.join(self.base_path,
                                        package.replace('.', os.path.sep))

            if not os.path.isdir(package_path):
                raise ValueError('Package not found: {0}'.format(package))

        if docs_path is not None and not skip_docs:
            if package is not None:
                docs_path = os.path.join(
                    docs_path, package.replace('.', os.path.sep))
            if not os.path.exists(docs_path):
                warnings.warn(
                    "Can not test .rst docs, since docs path "
                    "({0}) does not exist.".format(docs_path))
                docs_path = None

        if test_path:
            base, ext = os.path.splitext(test_path)

            if ext in ('.rst', ''):
                if docs_path is None:
                    # This shouldn't happen from "python setup.py test"
                    raise ValueError(
                        "Can not test .rst files without a docs_path "
                        "specified.")

                abs_docs_path = os.path.abspath(docs_path)
                abs_test_path = os.path.abspath(
                    os.path.join(abs_docs_path, os.pardir, test_path))

                common = os.path.commonprefix((abs_docs_path, abs_test_path))

                if os.path.exists(abs_test_path) and common == abs_docs_path:
                    # Since we aren't testing any Python files within
                    # the astropy tree, we need to forcibly load the
                    # astropy py.test plugins, and then turn on the
                    # doctest_rst plugin.
                    all_args.extend(['-p', 'astropy.tests.pytest_plugins',
                                     '--doctest-rst'])
                    test_path = abs_test_path

            if not (os.path.isdir(test_path) or ext in ('.py', '.rst')):
                raise ValueError("Test path must be a directory or a path to "
                                 "a .py or .rst file")

            all_args.append(test_path)
        else:
            all_args.append(package_path)
            if docs_path is not None and not skip_docs:
                all_args.extend([docs_path, '--doctest-rst'])

        # add any additional args entered by the user
        if args is not None:
            all_args.extend(
                shlex.split(args, posix=not sys.platform.startswith('win')))

        # add verbosity flag
        if verbose:
            all_args.append('-v')

        # turn on pastebin output
        if pastebin is not None:
            if pastebin in ['failed', 'all']:
                all_args.append('--pastebin={0}'.format(pastebin))
            else:
                raise ValueError("pastebin should be 'failed' or 'all'")

        # run @remote_data tests
        if remote_data:
            all_args.append('--remote-data')

        if pep8:
            try:
                import pytest_pep8
            except ImportError:
                raise ImportError('PEP8 checking requires pytest-pep8 plugin: '
                                  'http://pypi.python.org/pypi/pytest-pep8')
            else:
                all_args.extend(['--pep8', '-k', 'pep8'])

        # activate post-mortem PDB for failing tests
        if pdb:
            all_args.append('--pdb')

        # check for opened files after each test
        if open_files:
            if parallel != 0:
                raise SystemError(
                    "open file detection may not be used in conjunction with "
                    "parallel testing.")

            try:
                import psutil
            except ImportError:
                raise SystemError(
                    "open file detection requested, but psutil package "
                    "is not installed.")

            all_args.append('--open-files')

            print("Checking for unclosed files")

        if parallel != 0:
            try:
                import xdist
            except ImportError:
                raise ImportError(
                    'Parallel testing requires the pytest-xdist plugin '
                    'https://pypi.python.org/pypi/pytest-xdist')

            try:
                parallel = int(parallel)
            except ValueError:
                raise ValueError(
                    "parallel must be an int, got {0}".format(parallel))

            if parallel < 0:
                parallel = multiprocessing.cpu_count()
            all_args.extend(['-n', six.text_type(parallel)])

        if repeat:
            all_args.append('--repeat={0}'.format(repeat))

        if six.PY2:
            all_args = [x.encode('utf-8') for x in all_args]

        # override the config locations to not make a new directory nor use
        # existing cache or config
        xdg_config_home = os.environ.get('XDG_CONFIG_HOME')
        xdg_cache_home = os.environ.get('XDG_CACHE_HOME')
        astropy_config = tempfile.mkdtemp('astropy_config')
        astropy_cache = tempfile.mkdtemp('astropy_cache')
        os.environ[str('XDG_CONFIG_HOME')] = str(astropy_config)
        os.environ[str('XDG_CACHE_HOME')] = str(astropy_cache)
        os.mkdir(os.path.join(os.environ['XDG_CONFIG_HOME'], 'astropy'))
        os.mkdir(os.path.join(os.environ['XDG_CACHE_HOME'], 'astropy'))
        # To fully force configuration reloading from a different file (in this
        # case our default one in a temp directory), clear the config object
        # cache.
        configuration._cfgobjs.clear()

        # This prevents cyclical import problems that make it
        # impossible to test packages that define Table types on their
        # own.
        from ..table import Table

        try:
            result = pytest.main(args=all_args, plugins=plugins)
        finally:
            shutil.rmtree(os.environ['XDG_CONFIG_HOME'])
            shutil.rmtree(os.environ['XDG_CACHE_HOME'])
            if xdg_config_home is not None:
                os.environ[str('XDG_CONFIG_HOME')] = xdg_config_home
            else:
                del os.environ['XDG_CONFIG_HOME']
            if xdg_cache_home is not None:
                os.environ[str('XDG_CACHE_HOME')] = xdg_cache_home
            else:
                del os.environ['XDG_CACHE_HOME']
            configuration._cfgobjs.clear()

        return result

    run_tests.__doc__ = test.__doc__


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
        d = cov.data
        cov._harvest_data()
        for key in d.lines.keys():
            new_path = os.path.relpath(
                os.path.realpath(key),
                os.path.realpath(testing_path))
            new_path = os.path.abspath(
                os.path.join(rootdir, new_path))
            d.lines[new_path] = d.lines.pop(key)

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

    This can also be used a context manager, in which case it is just an alias
    for the `pytest.raises` context manager (because the two have the same name
    this help avoid confusion by being flexible).
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
    if not _deprecations_as_exceptions:
        return

    # First, totally reset the warning state
    for module in list(six.itervalues(sys.modules)):
        # We don't want to deal with six.MovedModules, only "real"
        # modules.
        if (isinstance(module, types.ModuleType) and
            hasattr(module, '__warningregistry__')):
            del module.__warningregistry__

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
            "always",
            r"functions overriding warnings\.showwarning\(\) must support "
            r"the 'line' argument",
            DeprecationWarning)

    if sys.version_info[:2] >= (3, 4):
        # py.test reads files with the 'U' flag, which is now
        # deprecated in Python 3.4.
        warnings.filterwarnings(
            "always",
            r"'U' mode is deprecated",
            DeprecationWarning)


class catch_warnings(warnings.catch_warnings):
    """
    A high-powered version of warnings.catch_warnings to use for testing
    and to make sure that there is no dependence on the order in which
    the tests are run.

    This completely blitzes any memory of any warnings that have
    appeared before so that all warnings will be caught and displayed.

    *args is a set of warning classes to collect.  If no arguments are
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


def assert_follows_unicode_guidelines(
        x, roundtrip=None):
    """
    Test that an object follows our Unicode policy.  See
    "Unicode Policy" in the coding guidelines.

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


def assert_quantity_allclose(actual, desired, rtol=1.e-7, atol=0, err_msg='', verbose=True):
    """
    Raise an assertion if two objects are not equal up to desired tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.testing.assert_allclose`.
    """

    import numpy as np
    from .. import units as u

    actual = u.Quantity(actual, subok=True, copy=False)

    desired = u.Quantity(desired, subok=True, copy=False)
    try:
        desired = desired.to(actual.unit)
    except u.UnitsError:
        raise u.UnitsError("Units for 'desired' ({0}) and 'actual' ({1}) are not convertible".format(desired.unit, actual.unit))

    if atol == 0:
        atol = u.Quantity(0)
    else:
        atol = u.Quantity(atol, subok=True, copy=False)
        try:
            atol = atol.to(actual.unit)
        except u.UnitsError:
            raise u.UnitsError("Units for 'atol' ({0}) and 'actual' ({1}) are not convertible".format(atol.unit, actual.unit))

    rtol =  u.Quantity(rtol, subok=True, copy=False)
    try:
        rtol = rtol.to(u.dimensionless_unscaled)
    except:
        raise u.UnitsError("`rtol` should be dimensionless")

    np.testing.assert_allclose(actual.value, desired.value,
                               rtol=rtol.value, atol=atol.value,
                               err_msg=err_msg, verbose=verbose)

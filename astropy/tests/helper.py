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
from ..utils.exceptions import AstropyWarning
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
                  docs_path=None, skip_docs=False):
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
            if ext == '.py':
                test_path = os.path.abspath(test_path)
                all_args.append(test_path)
            elif ext == '.rst':
                if docs_path is None:
                    # This shouldn't happen from "python setup.py test"
                    raise ValueError(
                        "Can not test .rst files without a docs_path specified.")
                else:
                    # Since we aren't testing any Python files within
                    # the astropy tree, we need to forcibly load the
                    # astropy py.test plugins, and then turn on the
                    # doctest_rst plugin.
                    all_args.extend(['-p', 'astropy.tests.pytest_plugins', '--doctest-rst'])
                    test_path = os.path.join(docs_path, '..', test_path)
                    all_args.append(test_path)
            else:
                raise ValueError("Test file path must be to a .py or .rst file")
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
            try:
                subproc = subprocess.Popen(
                    ['lsof -F0 -n -p {0}'.format(os.getpid())],
                    shell=True, stdout=subprocess.PIPE)
                output = subproc.communicate()[0].strip()
            except subprocess.CalledProcessError:
                raise SystemError(
                    "open file detection requested, but could not "
                    "successfully run the 'lsof' command")

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


def enable_deprecations_as_exceptions():
    """
    Turn on the feature that turns deprecations into exceptions.
    """
    global _deprecations_as_exceptions
    _deprecations_as_exceptions = True


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


##############################################################################
# Note: the following class exists only for backward-compatibility purposes. #
#       It has been moved to the separate astropy-helpers package, located   #
#       at https://github.com/astropy/astropy-helpers. Any new development   #
#       or bug fixes should be done there                                    #
##############################################################################


class astropy_test(Command, object):
    user_options = [
        ('package=', 'P',
         "The name of a specific package to test, e.g. 'io.fits' or 'utils'.  "
         "If nothing is specified, all default tests are run."),
        ('test-path=', 't',
         'Specify a test location by path.  If a relative path to a '
         '.py file, it is relative to the built package.  If a relative '
         'path to a .rst file, it is relative to the docs directory '
         '(see --docs-path).  May also be an absolute path.'),
        ('verbose-results', 'V',
         'Turn on verbose output from pytest.'),
        ('plugins=', 'p',
         'Plugins to enable when running pytest.'),
        ('pastebin=', 'b',
         "Enable pytest pastebin output. Either 'all' or 'failed'."),
        ('args=', 'a',
         'Additional arguments to be passed to pytest.'),
        ('remote-data', 'R', 'Run tests that download remote data.'),
        ('pep8', '8',
         'Enable PEP8 checking and disable regular tests. '
         'Requires the pytest-pep8 plugin.'),
        ('pdb', 'd',
         'Start the interactive Python debugger on errors.'),
        ('coverage', 'c',
         'Create a coverage report. Requires the coverage package.'),
        ('open-files', 'o', 'Fail if any tests leave files open.'),
        ('parallel=', 'j',
         'Run the tests in parallel on the specified number of '
         'CPUs.  If negative, all the cores on the machine will be '
         'used.  Requires the pytest-xdist plugin.'),
        ('docs-path=', None,
         'The path to the documentation .rst files.  If not provided, and '
         'the current directory contains a directory called "docs", that '
         'will be used.'),
        ('skip-docs', None,
         "Don't test the documentation .rst files.")
    ]

    user_options = _fix_user_options(user_options)

    package_name = None

    def initialize_options(self):
        self.package = None
        self.test_path = None
        self.verbose_results = False
        self.plugins = None
        self.pastebin = None
        self.args = None
        self.remote_data = False
        self.pep8 = False
        self.pdb = False
        self.coverage = False
        self.open_files = False
        self.parallel = 0
        self.docs_path = None
        self.skip_docs = False

    def finalize_options(self):
        # Normally we would validate the options here, but that's handled in
        # run_tests
        pass

    def run(self):
        self.reinitialize_command('build', inplace=False)
        self.run_command('build')
        build_cmd = self.get_finalized_command('build')
        new_path = os.path.abspath(build_cmd.build_lib)

        if self.docs_path is None:
            if os.path.exists('docs'):
                self.docs_path = os.path.abspath('docs')

        # Copy the build to a temporary directory for the purposes of testing
        # - this avoids creating pyc and __pycache__ directories inside the
        # build directory
        tmp_dir = tempfile.mkdtemp(prefix='astropy-test-')
        testing_path = os.path.join(tmp_dir, os.path.basename(new_path))
        shutil.copytree(new_path, testing_path)
        shutil.copy('setup.cfg', testing_path)

        cmd_pre = ''
        cmd_post = ''

        try:
            if self.coverage:
                if self.parallel != 0:
                    raise ValueError(
                        "--coverage can not be used with --parallel")

                try:
                    import coverage
                except ImportError:
                    raise ImportError(
                        "--coverage requires that the coverage package is "
                        "installed.")

                # Don't use get_pkg_data_filename here, because it
                # requires importing astropy.config and thus screwing
                # up coverage results for those packages.
                coveragerc = os.path.join(
                    testing_path, self.package_name, 'tests', 'coveragerc')

                # We create a coveragerc that is specific to the version
                # of Python we're running, so that we can mark branches
                # as being specifically for Python 2 or Python 3
                with open(coveragerc, 'rb') as fd:
                    coveragerc_content = fd.read().decode('utf-8')
                if six.PY3:
                    ignore_python_version = '2'
                elif six.PY2:
                    ignore_python_version = '3'
                coveragerc_content = coveragerc_content.replace(
                    "{ignore_python_version}", ignore_python_version).replace(
                        "{packagename}", self.package_name)
                tmp_coveragerc = os.path.join(tmp_dir, 'coveragerc')
                with open(tmp_coveragerc, 'wb') as tmp:
                    tmp.write(coveragerc_content.encode('utf-8'))

                cmd_pre = (
                    'import coverage; '
                    'cov = coverage.coverage(data_file="{0}", config_file="{1}"); '
                    'cov.start();'.format(
                        os.path.abspath(".coverage"), tmp_coveragerc))
                cmd_post = (
                    'cov.stop(); '
                    'from astropy.tests.helper import _save_coverage; '
                    '_save_coverage(cov, result, "{0}", "{1}");'.format(
                        os.path.abspath('.'), testing_path))

            if six.PY3:
                set_flag = "import builtins; builtins._ASTROPY_TEST_ = True"
            elif six.PY2:
                set_flag = "import __builtin__; __builtin__._ASTROPY_TEST_ = True"

            cmd = ('{cmd_pre}{0}; import {1.package_name}, sys; result = ('
                   '{1.package_name}.test('
                   'package={1.package!r}, '
                   'test_path={1.test_path!r}, '
                   'args={1.args!r}, '
                   'plugins={1.plugins!r}, '
                   'verbose={1.verbose_results!r}, '
                   'pastebin={1.pastebin!r}, '
                   'remote_data={1.remote_data!r}, '
                   'pep8={1.pep8!r}, '
                   'pdb={1.pdb!r}, '
                   'open_files={1.open_files!r}, '
                   'parallel={1.parallel!r}, '
                   'docs_path={1.docs_path!r}, '
                   'skip_docs={1.skip_docs!r})); '
                   '{cmd_post}'
                   'sys.exit(result)')
            cmd = cmd.format(set_flag, self, cmd_pre=cmd_pre, cmd_post=cmd_post)

            # Run the tests in a subprocess--this is necessary since
            # new extension modules may have appeared, and this is the
            # easiest way to set up a new environment

            # On Python 3.x prior to 3.3, the creation of .pyc files
            # is not atomic.  py.test jumps through some hoops to make
            # this work by parsing import statements and carefully
            # importing files atomically.  However, it can't detect
            # when __import__ is used, so its carefulness still fails.
            # The solution here (admittedly a bit of a hack), is to
            # turn off the generation of .pyc files altogether by
            # passing the `-B` switch to `python`.  This does mean
            # that each core will have to compile .py file to bytecode
            # itself, rather than getting lucky and borrowing the work
            # already done by another core.  Compilation is an
            # insignificant fraction of total testing time, though, so
            # it's probably not worth worrying about.
            retcode = subprocess.call([sys.executable, '-B', '-c', cmd],
                                      cwd=testing_path, close_fds=False)
        finally:
            # Remove temporary directory
            shutil.rmtree(tmp_dir)

        raise SystemExit(retcode)

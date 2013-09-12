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

if os.environ.get('ASTROPY_USE_SYSTEM_PYTEST') or '_pytest' in sys.modules:
    import pytest

else:
    from ..extern import pytest as extern_pytest

    if sys.version_info >= (3, 0):
        exec("def do_exec_def(co, loc): exec(co, loc)\n")
        extern_pytest.do_exec = do_exec_def

        unpacked_sources = extern_pytest.sources.encode("ascii")
        unpacked_sources = pickle.loads(
            zlib.decompress(base64.decodebytes(unpacked_sources)))
    else:
        exec("def do_exec_def(co, loc): exec co in loc\n")
        extern_pytest.do_exec = do_exec_def

        unpacked_sources = pickle.loads(
            zlib.decompress(base64.decodestring(extern_pytest.sources)))

    importer = extern_pytest.DictImporter(unpacked_sources)
    sys.meta_path.insert(0, importer)

    # On Python 3.1, we need to forcibly import the py.test-"bundled"
    # argparse before importing py.test, since it isn't in the
    # standard library, and py.test's workaround doesn't appear to
    # work with the "absolute imports" of Python 3.x.
    if sys.version_info[0] == 3 and sys.version_info[1] <= 1:
        argparse = importer.load_module(str('argparse'))
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
                  pdb=False, coverage=False, open_files=False, parallel=0):
        """
        The docstring for this method lives in astropy/__init__.py:test
        """
        if package is None:
            package_path = self.base_path
        else:
            package_path = os.path.join(self.base_path,
                                        package.replace('.', os.path.sep))

            if not os.path.isdir(package_path):
                raise ValueError('Package not found: {0}'.format(package))

        if test_path:
            package_path = os.path.join(package_path,
                                        os.path.abspath(test_path))

        all_args = package_path

        # add any additional args entered by the user
        if args is not None:
            all_args += ' {0}'.format(args)

        # add verbosity flag
        if verbose:
            all_args += ' -v'

        # turn on pastebin output
        if pastebin is not None:
            if pastebin in ['failed', 'all']:
                all_args += ' --pastebin={0}'.format(pastebin)
            else:
                raise ValueError("pastebin should be 'failed' or 'all'")

        # run @remote_data tests
        if remote_data:
            all_args += ' --remote-data'

        if pep8:
            try:
                import pytest_pep8
            except ImportError:
                raise ImportError('PEP8 checking requires pytest-pep8 plugin: '
                                  'http://pypi.python.org/pypi/pytest-pep8')
            else:
                all_args += ' --pep8 -k pep8'

        # activate post-mortem PDB for failing tests
        if pdb:
            all_args += ' --pdb'

        if coverage:
            try:
                import pytest_cov
            except ImportError:
                raise ImportError(
                    'Coverage reporting requires pytest-cov plugin: '
                    'http://pypi.python.org/pypi/pytest-cov')
            else:
                # Don't use get_pkg_data_filename here, because it
                # requires importing astropy.config and thus screwing
                # up coverage results for those packages.
                coveragerc = os.path.join(
                    os.path.dirname(__file__), 'coveragerc')

                # We create a coveragerc that is specific to the version
                # of Python we're running, so that we can mark branches
                # as being specifically for Python 2 or Python 3
                with open(coveragerc, 'r') as fd:
                    coveragerc_content = fd.read()
                if sys.version_info[0] >= 3:
                    ignore_python_version = '2'
                else:
                    ignore_python_version = '3'
                coveragerc_content = coveragerc_content.replace(
                    "{ignore_python_version}", ignore_python_version)
                with tempfile.NamedTemporaryFile(delete=False) as tmp:
                    tmp.write(coveragerc_content.encode('utf-8'))

                all_args += (
                    ' --cov-report html --cov astropy'
                    ' --cov-config {0}'.format(tmp.name))

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

            all_args += ' --open-files'

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
            all_args += ' -n {0}'.format(parallel)

        if sys.version_info < (2, 7, 3):
            all_args = all_args.encode('utf-8')

        try:
            all_args = shlex.split(
                all_args, posix=not sys.platform.startswith('win'))

            result = pytest.main(args=all_args, plugins=plugins)
        finally:
            if coverage:
                if not tmp.closed:
                    tmp.close()
                os.remove(tmp.name)

        return result

    run_tests.__doc__ = test.__doc__


class astropy_test(Command, object):
    user_options = [
        (str('package='), str('P'),
         str("The name of a specific package to test, e.g. 'io.fits' or 'utils'.  "
             "If nothing is specified all default Astropy tests are run.")),
        (str('test-path='), str('t'),
         str('Specify a test location by path. Must be '
             'specified absolutely or relative to the current directory. '
             'May be a single file or directory.')),
        (str('verbose-results'), str('V'),
         str('Turn on verbose output from pytest. Same as specifying `-v` in '
             '`args`.')),
        (str('plugins='), str('p'),
         str('Plugins to enable when running pytest.  Same as specifying `-p` in '
             '`args`.')),
        (str('pastebin='), str('b'),
         str("Enable pytest pastebin output. Either 'all' or 'failed'.")),
        (str('args='), str('a'),
         str('Additional arguments to be passed to pytest')),
        (str('remote-data'), str('R'), 'Run tests that download remote data'),
        (str('pep8'), str('8'),
         str('Enable PEP8 checking and disable regular tests. '
             'Same as specifying `--pep8 -k pep8` in `args`. Requires the '
             'pytest-pep8 plugin.')),
        (str('pdb'), str('d'),
         str('Turn on PDB post-mortem analysis for failing tests. '
             'Same as specifying `--pdb` in `args`.')),
        (str('coverage'), str('c'),
         str('Create a coverage report. Requires the pytest-cov '
             'plugin is installed')),
        (str('open-files'), str('o'), 'Fail if any tests leave files open'),
        (str('parallel='), str('n'),
         str('Run the tests in parallel on the specified '
             'number of CPUs.  If parallel is negative, it will use the all '
             'the cores on the machine.  Requires the `pytest-xdist` plugin '
             'is installed.'))

    ]

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

    def finalize_options(self):
        # Normally we would validate the options here, but that's handled in
        # run_tests
        pass

    def run(self):
        self.reinitialize_command('build', inplace=False)
        self.run_command('build')
        build_cmd = self.get_finalized_command('build')
        new_path = os.path.abspath(build_cmd.build_lib)

        # Copy the build to a temporary directory for the purposes of testing
        # - this avoids creating pyc and __pycache__ directories inside the
        # build directory
        tmp_dir = tempfile.mkdtemp(prefix='astropy-test-')
        testing_path = os.path.join(tmp_dir, os.path.basename(new_path))
        shutil.copytree(new_path, testing_path)
        shutil.copy('setup.cfg', testing_path)

        try:

            # Run the tests in a subprocess--this is necessary since new extension
            # modules may have appeared, and this is the easiest way to set up a
            # new environment

            # We need to set a flag in the child's environment so that
            # unnecessary code is not imported before py.test can start
            # up, otherwise the coverage results will be artifically low.
            if sys.version_info[0] >= 3:
                set_flag = "import builtins; builtins._ASTROPY_TEST_ = True"
            else:
                set_flag = "import __builtin__; __builtin__._ASTROPY_TEST_ = True"

            cmd = ('{0}; import {1.package_name}, sys; sys.exit('
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
                   'coverage={1.coverage!r}, '
                   'open_files={1.open_files!r}, '
                   'parallel={1.parallel!r}))')
            cmd = cmd.format(set_flag, self)

            # override the config locations to not make a new directory nor use
            # existing cache or config
            os.environ['XDG_CONFIG_HOME'] = tempfile.mkdtemp('astropy_config')
            os.environ['XDG_CACHE_HOME'] = tempfile.mkdtemp('astropy_cache')
            os.mkdir(os.path.join(os.environ['XDG_CONFIG_HOME'], 'astropy'))
            os.mkdir(os.path.join(os.environ['XDG_CACHE_HOME'], 'astropy'))

            try:
                retcode = subprocess.call([sys.executable, '-c', cmd],
                                          cwd=testing_path, close_fds=False)
            finally:
                # kill the temporary dirs
                shutil.rmtree(os.environ['XDG_CONFIG_HOME'])
                shutil.rmtree(os.environ['XDG_CACHE_HOME'])

            if self.coverage and retcode == 0:

                # Copy the htmlcov from build/lib.../htmlcov to a more
                # obvious place
                if os.path.exists('htmlcov'):
                    shutil.rmtree('htmlcov')
                shutil.copytree(os.path.join(testing_path, 'htmlcov'), 'htmlcov')

                # The coverage report includes the full path to the temporary
                # directory, so we replace all the paths with the true source
                # path. This means that the coverage line-by-line report will
                # only be correct for Python 2 code (since the Python 3 code
                # will be different in the build directory from the source
                # directory as long as 2to3 is needed). Therefore we only do
                # this fix for Python 2.x.

                if six.PY2:
                    import coverage
                    d = coverage.CoverageData()
                    d.read_file(os.path.join(testing_path, '.coverage'))
                    for key in d.lines.keys():
                        new_path = os.path.relpath(os.path.realpath(key),
                                                   os.path.realpath(testing_path))
                        new_path = os.path.abspath(new_path)
                        d.lines[new_path] = d.lines.pop(key)
                    d.write()
                else:
                    shutil.copy2(os.path.join(testing_path, '.coverage'), '.coverage')

        finally:

            # Remove temporary directory
            shutil.rmtree(tmp_dir)

        raise SystemExit(retcode)


class raises(object):
    """
    A decorator to mark that a test should raise a given exception.
    Use as follows::

        @raises(ZeroDivisionError)
        def test_foo():
            x = 1/0
    """
    # pep-8 naming exception -- this is a decorator class
    def __init__(self, exc):
        self._exc = exc

    def __call__(self, func):
        @functools.wraps(func)
        def run_raises_test(*args, **kwargs):
            pytest.raises(self._exc, func, *args, **kwargs)
        return run_raises_test


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
        for module in list(six.itervalues(sys.modules)):
            if hasattr(module, '__warningregistry__'):
                del module.__warningregistry__
        super(catch_warnings, self).__init__(record=True)
        self.classes = classes

    def __enter__(self):
        warning_list = super(catch_warnings, self).__enter__()
        warnings.resetwarnings()
        if len(self.classes) == 0:
            warnings.simplefilter('always')
        else:
            warnings.simplefilter('ignore')
            for cls in self.classes:
                warnings.simplefilter('always', cls)
        return warning_list


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
    from .. import UNICODE_OUTPUT
    from ..extern import six

    UNICODE_OUTPUT.set(False)

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

    UNICODE_OUTPUT.set(True)

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

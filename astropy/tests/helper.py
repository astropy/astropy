# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module prvoides the tools used to internally run the astropy test suite
from the installed astropy.  It makes use of the `pytest` testing framework.
"""

import shlex
import sys
import base64
import zlib
import functools
import os
import subprocess
import shutil

from distutils.core import Command

from .. import test

#If we are in setup.py, we don't want to import astropy.config
if __builtins__.get('_ASTROPY_SETUP_'):
    if os.environ.get('ASTROPY_USE_SYSTEM_PYTEST'):
        USE_SYSTEM_PYTEST = lambda: True
    else:
        USE_SYSTEM_PYTEST = lambda: False
else:
    from ..config import ConfigurationItem

    USE_SYSTEM_PYTEST = ConfigurationItem('use_system_pytest', False,
                                          'Set to True to load system pytest.  '
                                          'This item will *not* be obeyed if '
                                          'using setup.py.  In that case the '
                                          'environment variable '
                                          'ASTROPY_USE_SYSTEM_TEST must be '
                                          'used',
                                          'boolean', 'astropy.tests.helper')

if USE_SYSTEM_PYTEST():
    import pytest

else:
    from ..extern import pytest as extern_pytest

    if sys.version_info >= (3, 0):
        exec("def do_exec_def(co, loc): exec(co, loc)\n")
        extern_pytest.do_exec = do_exec_def

        import pickle
        unpacked_sources = extern_pytest.sources.encode("ascii")
        unpacked_sources = pickle.loads(
            zlib.decompress(base64.decodebytes(unpacked_sources)))
    else:
        exec("def do_exec_def(co, loc): exec co in loc\n")
        extern_pytest.do_exec = do_exec_def

        import cPickle as pickle
        unpacked_sources = pickle.loads(
            zlib.decompress(base64.decodestring(extern_pytest.sources)))

    importer = extern_pytest.DictImporter(unpacked_sources)
    sys.meta_path.append(importer)

    pytest = importer.load_module('pytest')


# pytest marker to mark tests which get data from the web
remote_data = pytest.mark.remote_data


class TestRunner(object):
    def __init__(self, base_path):
        self.base_path = base_path

    def run_tests(self, package=None, test_path=None, args=None, plugins=None,
                  verbose=False, pastebin=None, remote_data=False, pep8=False,
                  pdb=False, coverage=False):
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
                from .. import config
                all_args += (
                    ' --cov-report html --cov astropy '
                    '--cov-config {0}'.format(
                        config.get_data_filename('coveragerc')))

        all_args = shlex.split(
            all_args, posix=not sys.platform.startswith('win'))

        result = pytest.main(args=all_args, plugins=plugins)
    run_tests.__doc__ = test.__doc__


class astropy_test(Command, object):
    user_options = [
        ('package=', 'P',
         "The name of a specific package to test, e.g. 'io.fits' or 'utils'.  "
         "If nothing is specified all default Astropy tests are run."),
        ('test-path=', 't', 'Specify a test location by path. Must be '
         'specified absolutely or relative to the current directory. '
         'May be a single file or directory.'),
        ('verbose-results', 'V',
         'Turn on verbose output from pytest. Same as specifying `-v` in '
         '`args`.'),
        ('plugins=', 'p',
         'Plugins to enable when running pytest.  Same as specifying `-p` in '
         '`args`.'),
        ('pastebin=', 'b',
         "Enable pytest pastebin output. Either 'all' or 'failed'."),
        ('args=', 'a', 'Additional arguments to be passed to pytest'),
        ('remote-data', 'R', 'Run tests that download remote data'),
        ('pep8', '8', 'Enable PEP8 checking and disable regular tests. '
         'Same as specifying `--pep8 -k pep8` in `args`. Requires the '
         'pytest-pep8 plugin.'),
        ('pdb', 'd', 'Turn on PDB post-mortem analysis for failing tests. '
         'Same as specifying `--pdb` in `args`.'),
        ('coverage', 'c', 'Create a coverage report. Requires the pytest-cov '
         'plugin is installed')
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

    def finalize_options(self):
        # Normally we would validate the options here, but that's handled in
        # run_tests
        pass

    def run(self):
        self.reinitialize_command('build', inplace=False)
        self.run_command('build')
        build_cmd = self.get_finalized_command('build')
        new_path = os.path.abspath(build_cmd.build_lib)

        # Run the tests in a subprocess--this is necessary since new extension
        # modules may have appeared, and this is the easiest way to set up a
        # new environment
        cmd = ('import {0}, sys; sys.exit({0}.test({1!r}, {2!r}, ' +
               '{3!r}, {4!r}, {5!r}, {6!r}, {7!r}, {8!r}, {9!r}, {10!r}))')
        cmd = cmd.format(self.package_name,
                         self.package, self.test_path, self.args,
                         self.plugins, self.verbose_results, self.pastebin,
                         self.remote_data, self.pep8, self.pdb, self.coverage)

        retcode = subprocess.call([sys.executable, '-c', cmd],
                                  cwd=new_path, close_fds=False)

        if self.coverage and retcode == 0:
            # Copy the htmlcov from build/lib.../htmlcov to a more
            # obvious place
            if os.path.exists('htmlcov'):
                shutil.rmtree('htmlcov')
            shutil.copytree(os.path.join(new_path, 'htmlcov'), 'htmlcov')

        raise SystemExit(retcode)


class raises:
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

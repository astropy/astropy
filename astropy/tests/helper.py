# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module prvoides the tools used to internally run the astropy test suite
from the installed astropy.  It makes use of the `pytest` testing framework.
"""

import sys
import base64
import zlib
import functools
import os
import subprocess

from distutils.core import Command

try:
    import pytest

    # Check that a recent py.test version is being used
    from distutils import version
    if version.LooseVersion(pytest.__version__) < \
       version.LooseVersion('2.2.0'):
        class VersionError(Exception):
            pass
        raise VersionError("py.test 2.2.0 or later is required")

except ImportError:
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
                  verbose=False, pastebin=None, remote_data=False, pep8=False):
        """
        Run Astropy tests using py.test. A proper set of arguments is
        constructed and passed to `pytest.main`.

        Parameters
        ----------
        package : str, optional
            The name of a specific package to test, e.g. 'io.fits' or 'utils'.
            If nothing is specified all default Astropy tests are run.

        test_path : str, optional
            Specify location to test by path. May be a single file or
            directory. Must be specified absolutely or relative to the
            calling directory.

        args : str, optional
            Additional arguments to be passed to `pytest.main` in the `args`
            keyword argument.

        plugins : list, optional
            Plugins to be passed to `pytest.main` in the `plugins` keyword
            argument.

        verbose : bool, optional
            Convenience option to turn on verbose output from py.test. Passing
            True is the same as specifying `-v` in `args`.

        pastebin : {'failed','all',None}, optional
            Convenience option for turning on py.test pastebin output. Set to
            'failed' to upload info for failed tests, or 'all' to upload info
            for all tests.

        remote_data : bool, optional
            Controls whether to run tests marked with @remote_data. These
            tests use online data and are not run by default. Set to True to
            run these tests.

        pep8 : bool, optional
            Turn on PEP8 checking via the pytest-pep8 plugin and disable normal
            tests. Same as specifying `--pep8 -k pep8` in `args`.

        See Also
        --------
        pytest.main : py.test function wrapped by `run_tests`.

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

        return pytest.main(args=all_args, plugins=plugins)


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
         'pytest-pep8 plugin.')
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

    def finalize_options(self):
        # Normally we would validate the options here, but that's handled in
        # run_tests
        pass

    def run(self):
        self.reinitialize_command('build_py', inplace=False)
        self.run_command('build_py')
        if sys.version_info[0] >= 3:
            build_py_cmd = self.get_finalized_command('build_py')
            new_path = os.path.abspath(build_py_cmd.build_lib)
        else:
            new_path = os.getcwd()

        self.reinitialize_command('build_ext', inplace=True)
        self.run_command('build_ext')

        # Run the tests in a subprocess--this is necessary since new extension
        # modules may have appeared, and this is the easiest way to set up a
        # new environment
        cmd = ('import {0}, sys; sys.exit({0}.test({1!r}, {2!r}, ' +
               '{3!r}, {4!r}, {5!r}, {6!r}, {7!r}, {8!r}))')
        cmd = cmd.format(self.package_name,
                         self.package, self.test_path, self.args,
                         self.plugins, self.verbose_results, self.pastebin,
                         self.remote_data, self.pep8)

        raise SystemExit(subprocess.call([sys.executable, '-c', cmd],
                                         cwd=new_path))


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

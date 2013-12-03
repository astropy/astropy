"""Implements the Astropy TestRunner which is a thin wrapper around py.test."""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import atexit
import multiprocessing
import os
import pkgutil
import shlex
import shutil
import subprocess
import sys
import tempfile
import warnings
import zipfile
import zipimport

from ..config.paths import set_temp_config, set_temp_cache
from ..extern import six
from ..utils import wraps, find_current_module
from ..utils.exceptions import AstropyUserWarning


class TestRunner(object):
    def __init__(self, base_path):
        # base_path should be the root package's directory
        # testing_path is assumed to be the root directory that the package
        # itself is in (such as the source checkout)
        self.base_path = os.path.abspath(base_path)
        self.testing_path = os.path.dirname(self.base_path)
        self.package_name = os.path.basename(self.base_path)

    def run_tests(self, package=None, test_path=None, args=None, plugins=None,
                  verbose=False, pastebin=None, remote_data=False, pep8=False,
                  pdb=False, coverage=False, open_files=False, parallel=0,
                  docs_path=None, skip_docs=False, repeat=None):
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

        pdb : bool, optional
            Turn on PDB post-mortem analysis for failing tests. Same as
            specifying `--pdb` in `args`.

        open_files : bool, optional
            Fail when any tests leave files open.  Off by default, because
            this adds extra run time to the test suite.  Requires the
            ``psutil`` package.

        parallel : int, optional
            When provided, run the tests in parallel on the specified
            number of CPUs.  If parallel is negative, it will use the all
            the cores on the machine.  Requires the `pytest-xdist` plugin.

        docs_path : str, optional
            The path to the documentation .rst files.

        skip_docs : bool, optional
            When `True`, skips running the doctests in the .rst files.

        repeat : int, optional
            If set, specifies how many times each test should be run. This is
            useful for diagnosing sporadic failures.

        See Also
        --------
        pytest.main : py.test function wrapped by `run_tests`.
        """

        # Don't import pytest until it's actually needed to run the tests
        from .helper import pytest

        if coverage:
            warnings.warn(
                "The coverage option is ignored on run_tests, since it "
                "can not be made to work in that context.  Use "
                "'python setup.py test --coverage' instead.",
                AstropyUserWarning)

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
                import pytest_pep8  # pylint: disable=W0611
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
                import psutil  # pylint: disable=W0611
            except ImportError:
                raise SystemError(
                    "open file detection requested, but psutil package "
                    "is not installed.")

            all_args.append('--open-files')

            print("Checking for unclosed files")

        if parallel != 0:
            try:
                import xdist  # pylint: disable=W0611
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
        astropy_config = tempfile.mkdtemp('astropy_config')
        astropy_cache = tempfile.mkdtemp('astropy_cache')

        # This prevents cyclical import problems that make it
        # impossible to test packages that define Table types on their
        # own.
        from ..table import Table  # pylint: disable=W0611

        # Have to use nested with statements for cross-Python support
        # Note, using these context managers here is superfluous if the
        # config_dir or cache_dir options to py.test are in use, but it's
        # also harmless to nest the contexts
        with set_temp_config(astropy_config, delete=True):
            with set_temp_cache(astropy_cache, delete=True):
                return pytest.main(args=all_args, plugins=plugins)

    @wraps(run_tests)
    def run_tests_in_subprocess(self, **kwargs):
        """
        Launches a new Python interpreter in a subprocess and uses that
        subprocess to import the package and run ``package.test()``.  This is
        generally required for py.test to run the tests outside of the
        package's normal import path, for example from a temp directory.

        Arguments are the same as `TestRunner.run_tests`.
        """

        # Construct this modules testing command
        cmd = self._generate_testing_command(**kwargs)

        # Run the tests in a subprocess--this is necessary since
        # new extension modules may have appeared, and this is the
        # easiest way to set up a new environment
        # (Note: Previously this used the -B option to disable generation
        # of .pyc files to work around an issue on Python 3.0 - 3.2 with
        # atomicity of .pyc creation; however we no longer support those
        # Python versions so the workaround has been dropped)
        return subprocess.call([sys.executable, '-c', cmd],
                               cwd=self.testing_path, close_fds=False)

    @classmethod
    def make_test_runner_in(cls, path):
        """
        Constructs a `TestRunner` to run in the given path, and returns a
        ``test()`` function which takes the same arguments as
        `TestRunner.run_tests`.

        The returned ``test()`` function will be defined in the module this
        was called from.  This is used to implement the ``astropy.test()``
        function (or the equivalent for affiliated packages).
        """

        # Check to see if the package was imported from a zip file; if so we
        # need to extract it to a temporary directory in order for py.test to
        # find and run the tests
        importer = pkgutil.get_importer(path)
        if isinstance(importer, zipimport.zipimporter):
            pkg_name = os.path.basename(path)
            test_tempdir = tempfile.mkdtemp(prefix='{0}-test'.format(pkg_name))
            # If the tests result in a segfault or something like that,
            # this will leave junk in the temp dir, but not much we can do
            # about that
            atexit.register(lambda p: os.path.exists(p) and shutil.rmtree(p),
                            test_tempdir)
            with zipfile.ZipFile(importer.archive) as archive:
                archive.extractall(test_tempdir)

            runner = TestRunner(os.path.join(test_tempdir, pkg_name))
            test_meth = runner.run_tests_in_subprocess
        else:
            runner = cls(path)
            test_meth = runner.run_tests

        @wraps(test_meth, ('__doc__',), exclude_args=('self',))
        def test(*args, **kwargs):
            return test_meth(*args, **kwargs)

        module = find_current_module(2)
        if module is not None:
            test.__module__ = module.__name__

        # A somewhat unusual hack, but delete the attached __wrapped__
        # attribute--although this is normally used to tell if the function
        # was wrapped with wraps, on some version of Python this is also
        # used to determine the signature to display in help() which is
        # not useful in this case.  We don't really care in this case if the
        # function was wrapped either
        if hasattr(test, '__wrapped__'):
            del test.__wrapped__

        return test

    def _generate_testing_command(self, **kwargs):
        """
        Build a Python script to run the tests.
        """

        cmd_pre = ''  # Commands to run before the test function
        cmd_post = ''  # Commands to run after the test function

        if kwargs.get('coverage'):
            pre, post = self._generate_coverage_commands(**kwargs)
            cmd_pre += pre
            cmd_post += post

        if six.PY3:
            set_flag = "import builtins; builtins._ASTROPY_TEST_ = True"
        else:
            set_flag = "import __builtin__; __builtin__._ASTROPY_TEST_ = True"

        test_kwargs = ','.join('{0}={1!r}'.format(key, value)
                               for key, value in six.iteritems(kwargs))

        cmd = ('{cmd_pre}{0}; import {1}, sys; result = ('
               '{1}.test({2})); {cmd_post}sys.exit(result)')

        return cmd.format(set_flag, self.package_name, test_kwargs,
                          cmd_pre=cmd_pre, cmd_post=cmd_post)

    def _generate_coverage_commands(self, **kwargs):
        """
        This method creates the post and pre commands if coverage is to be
        generated
        """
        if kwargs.get('parallel'):
            raise ValueError(
                "--coverage can not be used with --parallel")

        try:
            import coverage  # pylint: disable=W0611
        except ImportError:
            raise ImportError(
                "--coverage requires that the coverage package is "
                "installed.")

        # Don't use get_pkg_data_filename here, because it
        # requires importing astropy.config and thus screwing
        # up coverage results for those packages.
        coveragerc = os.path.join(
            self.base_path, 'tests', 'coveragerc')

        # We create a coveragerc that is specific to the version
        # of Python we're running, so that we can mark branches
        # as being specifically for Python 2 or Python 3
        with open(coveragerc, 'r') as fd:
            coveragerc_content = fd.read()
        if six.PY3:
            ignore_python_version = '2'
        else:
            ignore_python_version = '3'
        coveragerc_content = coveragerc_content.replace(
            "{ignore_python_version}", ignore_python_version).replace(
                "{packagename}", self.package_name)
        tmp_coveragerc = os.path.join(self.base_path, 'coveragerc')
        with open(tmp_coveragerc, 'wb') as tmp:
            tmp.write(coveragerc_content.encode('utf-8'))

        cmd_pre = (
            'import coverage; '
            'cov = coverage.coverage(data_file="{0}", config_file="{1}"); '
            'cov.start(); '.format(
                os.path.abspath(".coverage"), tmp_coveragerc))
        cmd_post = (
            'cov.stop(); '
            'from astropy.tests.helper import _save_coverage; '
            '_save_coverage(cov, result, "{0}", "{1}"); '.format(
                os.path.abspath('.'), self.base_path))

        return cmd_pre, cmd_post

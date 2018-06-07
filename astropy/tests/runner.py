"""Implements the Astropy TestRunner which is a thin wrapper around py.test."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import inspect
import os
import copy
import shlex
import sys
import tempfile
import warnings
from collections import OrderedDict

from ..config.paths import set_temp_config, set_temp_cache
from ..extern import six
from ..utils import wraps, find_current_module
from ..utils.exceptions import AstropyWarning, AstropyDeprecationWarning

__all__ = ['TestRunner', 'TestRunnerBase', 'keyword']


class keyword(object):
    """
    A decorator to mark a method as keyword argument for the ``TestRunner``.

    Parameters
    ----------
    default_value : `object`
        The default value for the keyword argument. (Default: `None`)

    priority : `int`
        keyword argument methods are executed in order of descending priority.
    """

    def __init__(self, default_value=None, priority=0):
        self.default_value = default_value
        self.priority = priority

    def __call__(self, f):
        def keyword(*args, **kwargs):
            return f(*args, **kwargs)

        keyword._default_value = self.default_value
        keyword._priority = self.priority
        # Set __doc__ explicitly here rather than using wraps because we want
        # to keep the function name as keyword so we can inspect it later.
        keyword.__doc__ = f.__doc__

        return keyword


class TestRunnerBase(object):
    """
    The base class for the TestRunner.

    A test runner can be constructed by creating a subclass of this class and
    defining 'keyword' methods. These are methods that have the
    `~astropy.tests.runner.keyword` decorator, these methods are used to
    construct allowed keyword arguments to the
    `~astropy.tests.runner.TestRunnerBase.run_tests` method as a way to allow
    customization of individual keyword arguments (and associated logic)
    without having to re-implement the whole
    `~astropy.tests.runner.TestRunnerBase.run_tests` method.

    Examples
    --------

    A simple keyword method::

        class MyRunner(TestRunnerBase):

            @keyword('default_value'):
            def spam(self, spam, kwargs):
                \"\"\"
                spam : `str`
                    The parameter description for the run_tests docstring.
                \"\"\"
                # Return value must be a list with a CLI parameter for pytest.
                return ['--spam={}'.format(spam)]
    """

    def __init__(self, base_path):
        self.base_path = os.path.abspath(base_path)

    def __new__(cls, *args, **kwargs):
        # Before constructing the class parse all the methods that have been
        # decorated with ``keyword``.

        # The objective of this method is to construct a default set of keyword
        # arguments to the ``run_tests`` method. It does this by inspecting the
        # methods of the class for functions with the name ``keyword`` which is
        # the name of the decorator wrapping function. Once it has created this
        # dictionary, it also formats the docstring of ``run_tests`` to be
        # comprised of the docstrings for the ``keyword`` methods.

        # To add a keyword argument to the ``run_tests`` method, define a new
        # method decorated with ``@keyword`` and with the ``self, name, kwargs``
        # signature.
        # Get all 'function' members as the wrapped methods are functions
        if six.PY2:
            functions = inspect.getmembers(cls, predicate=inspect.ismethod)
        else:
            functions = inspect.getmembers(cls, predicate=inspect.isfunction)

        # Filter out anything that's not got the name 'keyword'
        keywords = filter(lambda func: func[1].__name__ == 'keyword', functions)
        # Sort all keywords based on the priority flag.
        sorted_keywords = sorted(keywords, key=lambda x: x[1]._priority, reverse=True)

        cls.keywords = OrderedDict()
        doc_keywords = ""
        for name, func in sorted_keywords:
            # Here we test if the function has been overloaded to return
            # NotImplemented which is the way to disable arguments on
            # subclasses. If it has been disabled we need to remove it from the
            # default keywords dict. We do it in the try except block because
            # we do not have access to an instance of the class, so this is
            # going to error unless the method is just doing `return
            # NotImplemented`.
            try:
                # Second argument is False, as it is normally a bool.
                # The other two are placeholders for objects.
                if func(None, False, None) is NotImplemented:
                    continue
            except Exception:
                pass

            # Construct the default kwargs dict and docstring
            cls.keywords[name] = func._default_value
            if func.__doc__:
                doc_keywords += ' '*8
                doc_keywords += func.__doc__.strip()
                doc_keywords += '\n\n'

        if six.PY2:
            cls.run_tests.__func__.__doc__ = cls.RUN_TESTS_DOCSTRING.format(keywords=doc_keywords)
        else:
            cls.run_tests.__doc__ = cls.RUN_TESTS_DOCSTRING.format(keywords=doc_keywords)

        return super(TestRunnerBase, cls).__new__(cls)

    def _generate_args(self, **kwargs):
        # Update default values with passed kwargs
        # but don't modify the defaults
        keywords = copy.deepcopy(self.keywords)
        keywords.update(kwargs)
        # Iterate through the keywords (in order of priority)
        args = []
        for keyword in keywords.keys():
            func = getattr(self, keyword)
            result = func(keywords[keyword], keywords)

            # Allow disabaling of options in a subclass
            if result is NotImplemented:
                raise TypeError("run_tests() got an unexpected keyword argument {}".format(keyword))

            # keyword methods must return a list
            if not isinstance(result, list):
                raise TypeError("{} keyword method must return a list".format(keyword))

            args += result

        if six.PY2:
            args = [x.encode('utf-8') for x in args]

        return args

    RUN_TESTS_DOCSTRING = \
        """
        Run the tests for the package.

        This method builds arguments for and then calls ``pytest.main``.

        Parameters
        ----------
        {keywords}

        """

    def run_tests(self, **kwargs):
        # The docstring for this method is defined as a class variable.
        # This allows it to be built for each subclass in __new__.

        # Don't import pytest until it's actually needed to run the tests
        import pytest

        # Raise error for undefined kwargs
        allowed_kwargs = set(self.keywords.keys())
        passed_kwargs = set(kwargs.keys())
        if not passed_kwargs.issubset(allowed_kwargs):
            wrong_kwargs = list(passed_kwargs.difference(allowed_kwargs))
            raise TypeError("run_tests() got an unexpected keyword argument {}".format(wrong_kwargs[0]))

        args = self._generate_args(**kwargs)

        # override the config locations to not make a new directory nor use
        # existing cache or config
        astropy_config = tempfile.mkdtemp('astropy_config')
        astropy_cache = tempfile.mkdtemp('astropy_cache')

        # Have to use nested with statements for cross-Python support
        # Note, using these context managers here is superfluous if the
        # config_dir or cache_dir options to py.test are in use, but it's
        # also harmless to nest the contexts
        with set_temp_config(astropy_config, delete=True):
            with set_temp_cache(astropy_cache, delete=True):
                return pytest.main(args=args, plugins=self.keywords['plugins'])

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

        runner = cls(path)

        @wraps(runner.run_tests, ('__doc__',), exclude_args=('self',))
        def test(**kwargs):
            return runner.run_tests(**kwargs)

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


class TestRunner(TestRunnerBase):
    """
    A test runner for astropy tests
    """

    # Increase priority so this warning is displayed first.
    @keyword(priority=1000)
    def coverage(self, coverage, kwargs):
        if coverage:
            warnings.warn(
                "The coverage option is ignored on run_tests, since it "
                "can not be made to work in that context.  Use "
                "'python setup.py test --coverage' instead.",
                AstropyWarning)

        return []

    # test_path depends on self.package_path so make sure this runs before
    # test_path.
    @keyword(priority=1)
    def package(self, package, kwargs):
        """
        package : str, optional
            The name of a specific package to test, e.g. 'io.fits' or 'utils'.
            If nothing is specified all default Astropy tests are run.
        """
        if package is None:
            self.package_path = self.base_path
        else:
            self.package_path = os.path.join(self.base_path,
                                        package.replace('.', os.path.sep))

            if not os.path.isdir(self.package_path):
                raise ValueError('Package not found: {0}'.format(package))

        if not kwargs['test_path']:
            return [self.package_path]

        return []

    @keyword()
    def test_path(self, test_path, kwargs):
        """
        test_path : str, optional
            Specify location to test by path. May be a single file or
            directory. Must be specified absolutely or relative to the
            calling directory.
        """
        all_args = []
        # Ensure that the package kwarg has been run.
        self.package(kwargs['package'], kwargs)
        if test_path:
            base, ext = os.path.splitext(test_path)

            if ext in ('.rst', ''):
                if kwargs['docs_path'] is None:
                    # This shouldn't happen from "python setup.py test"
                    raise ValueError(
                        "Can not test .rst files without a docs_path "
                        "specified.")

                abs_docs_path = os.path.abspath(kwargs['docs_path'])
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

            return all_args + [test_path]

        return []

    @keyword()
    def args(self, args, kwargs):
        """
        args : str, optional
            Additional arguments to be passed to ``pytest.main`` in the ``args``
            keyword argument.
        """
        if args:
            return shlex.split(args, posix=not sys.platform.startswith('win'))

        return []

    @keyword()
    def plugins(self, plugins, kwargs):
        """
        plugins : list, optional
            Plugins to be passed to ``pytest.main`` in the ``plugins`` keyword
            argument.
        """
        return []

    @keyword()
    def verbose(self, verbose, kwargs):
        """
        verbose : bool, optional
            Convenience option to turn on verbose output from py.test. Passing
            True is the same as specifying ``-v`` in ``args``.
        """
        if verbose:
            return ['-v']

        return []

    @keyword()
    def pastebin(self, pastebin, kwargs):
        """
        pastebin : ('failed', 'all', None), optional
            Convenience option for turning on py.test pastebin output. Set to
            'failed' to upload info for failed tests, or 'all' to upload info
            for all tests.
        """
        if pastebin is not None:
            if pastebin in ['failed', 'all']:
                return ['--pastebin={0}'.format(pastebin)]
            else:
                raise ValueError("pastebin should be 'failed' or 'all'")

        return []

    @keyword(default_value='none')
    def remote_data(self, remote_data, kwargs):
        """
        remote_data : {'none', 'astropy', 'any'}, optional
            Controls whether to run tests marked with @remote_data. This can be
            set to run no tests with remote data (``none``), only ones that use
            data from http://data.astropy.org (``astropy``), or all tests that
            use remote data (``any``). The default is ``none``.
        """

        if remote_data is True:
            remote_data = 'any'
        elif remote_data is False:
            remote_data = 'none'
        elif remote_data not in ('none', 'astropy', 'any'):
            warnings.warn("The remote_data option should be one of "
                          "none/astropy/any (found {0}). For backward-compatibility, "
                          "assuming 'any', but you should change the option to be "
                          "one of the supported ones to avoid issues in "
                          "future.".format(remote_data),
                          AstropyDeprecationWarning)
            remote_data = 'any'

        return ['--remote-data={0}'.format(remote_data)]

    @keyword()
    def pep8(self, pep8, kwargs):
        """
        pep8 : bool, optional
            Turn on PEP8 checking via the pytest-pep8 plugin and disable normal
            tests. Same as specifying ``--pep8 -k pep8`` in ``args``.
        """
        if pep8:
            try:
                import pytest_pep8  # pylint: disable=W0611
            except ImportError:
                raise ImportError('PEP8 checking requires pytest-pep8 plugin: '
                                  'http://pypi.python.org/pypi/pytest-pep8')
            else:
                return ['--pep8', '-k', 'pep8']

        return []

    @keyword()
    def pdb(self, pdb, kwargs):
        """
        pdb : bool, optional
            Turn on PDB post-mortem analysis for failing tests. Same as
            specifying ``--pdb`` in ``args``.
        """
        if pdb:
            return ['--pdb']
        return []

    @keyword()
    def open_files(self, open_files, kwargs):
        """
        open_files : bool, optional
            Fail when any tests leave files open.  Off by default, because
            this adds extra run time to the test suite.  Requires the
            ``psutil`` package.
        """
        if open_files:
            if kwargs['parallel'] != 0:
                raise SystemError(
                    "open file detection may not be used in conjunction with "
                    "parallel testing.")

            try:
                import psutil  # pylint: disable=W0611
            except ImportError:
                raise SystemError(
                    "open file detection requested, but psutil package "
                    "is not installed.")

            return ['--open-files']

            print("Checking for unclosed files")

        return []

    @keyword(0)
    def parallel(self, parallel, kwargs):
        """
        parallel : int or 'auto', optional
            When provided, run the tests in parallel on the specified
            number of CPUs.  If parallel is ``'auto'``, it will use the all
            the cores on the machine.  Requires the ``pytest-xdist`` plugin.
        """
        if parallel != 0:
            try:
                from xdist import plugin # noqa
            except ImportError:
                raise SystemError(
                    "running tests in parallel requires the pytest-xdist package")

            return ['-n', six.text_type(parallel)]

        return []

    @keyword()
    def docs_path(self, docs_path, kwargs):
        """
        docs_path : str, optional
            The path to the documentation .rst files.
        """
        if docs_path is not None and not kwargs['skip_docs']:
            if kwargs['package'] is not None:
                docs_path = os.path.join(
                    docs_path, kwargs['package'].replace('.', os.path.sep))
            if not os.path.exists(docs_path):
                warnings.warn(
                    "Can not test .rst docs, since docs path "
                    "({0}) does not exist.".format(docs_path))
                docs_path = None
        if docs_path and not kwargs['skip_docs'] and not kwargs['test_path']:
            return [docs_path, '--doctest-rst']

        return []

    @keyword()
    def skip_docs(self, skip_docs, kwargs):
        """
        skip_docs : `bool`, optional
            When `True`, skips running the doctests in the .rst files.
        """
        # Skip docs is a bool used by docs_path only.
        return []

    @keyword()
    def repeat(self, repeat, kwargs):
        """
        repeat : `int`, optional
            If set, specifies how many times each test should be run. This is
            useful for diagnosing sporadic failures.
        """
        if repeat:
            return ['--repeat={0}'.format(repeat)]

        return []

    # Override run_tests for astropy-specific fixes
    def run_tests(self, **kwargs):

        # This prevents cyclical import problems that make it
        # impossible to test packages that define Table types on their
        # own.
        from ..table import Table  # pylint: disable=W0611

        return super(TestRunner, self).run_tests(**kwargs)

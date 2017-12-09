# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This plugin provides advanced doctest support and enables the testing of .rst
files.
"""
from ...bundled import six

import doctest
import fnmatch
import imp
import os
import re
import sys

import pytest

from .output_checker import OutputChecker, FIX


# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.
def pytest_addoption(parser):

    parser.addoption("--doctest-plus", action="store_true",
                     help="enable running doctests with additional "
                     "features not found in the normal doctest "
                     "plugin")

    parser.addoption("--doctest-rst", action="store_true",
                     help="enable running doctests in .rst documentation")

    parser.addini("doctest_plus", "enable running doctests with additional "
                  "features not found in the normal doctest plugin")

    parser.addini("doctest_norecursedirs",
                  "like the norecursedirs option but applies only to doctest "
                  "collection", type="args", default=())

    parser.addini("doctest_rst",
                  "Run the doctests in the rst documentation",
                  default=False)


# We monkey-patch in our replacement doctest OutputChecker.  Not
# great, but there isn't really an API to replace the checker when
# using doctest.testfile, unfortunately.
doctest.OutputChecker = OutputChecker

REMOTE_DATA = doctest.register_optionflag('REMOTE_DATA')


def pytest_configure(config):

    doctest_plugin = config.pluginmanager.getplugin('doctest')
    if (doctest_plugin is None or config.option.doctestmodules or not
            (config.getini('doctest_plus') or config.option.doctest_plus)):
        return

    # These are the default doctest options we use for everything.
    # There shouldn't be any need to manually put them in doctests
    # themselves.
    opts = (doctest.ELLIPSIS |
            doctest.NORMALIZE_WHITESPACE |
            FIX)

    class DocTestModulePlus(doctest_plugin.DoctestModule):
        # pytest 2.4.0 defines "collect".  Prior to that, it defined
        # "runtest".  The "collect" approach is better, because we can
        # skip modules altogether that have no doctests.  However, we
        # need to continue to override "runtest" so that the built-in
        # behavior (which doesn't do whitespace normalization or
        # handling __doctest_skip__) doesn't happen.
        def collect(self):
            # When running directly from pytest we need to make sure that we
            # don't accidentally import setup.py!
            if self.fspath.basename == "setup.py":
                return
            elif self.fspath.basename == "conftest.py":
                try:
                    module = self.config._conftest.importconftest(self.fspath)
                except AttributeError:  # pytest >= 2.8.0
                    module = self.config.pluginmanager._importconftest(self.fspath)
            else:
                try:
                    module = self.fspath.pyimport()
                    # Just ignore searching modules that can't be imported when
                    # collecting doctests
                except ImportError:
                    return

            # uses internal doctest module parsing mechanism
            finder = DocTestFinderPlus()
            runner = doctest.DebugRunner(
                verbose=False, optionflags=opts, checker=OutputChecker())
            for test in finder.find(module):
                if test.examples:  # skip empty doctests
                    if config.getoption('remote_data', 'none') != 'any':
                        for example in test.examples:
                            if example.options.get(REMOTE_DATA):
                                example.options[doctest.SKIP] = True

                    yield doctest_plugin.DoctestItem(
                        test.name, self, runner, test)

    class DocTestTextfilePlus(doctest_plugin.DoctestItem, pytest.Module):

        # Some pytest plugins such as hypothesis try and access the 'obj'
        # attribute, and by default this returns an error for this class
        # so we override it here to avoid any issues.
        def obj(self):
            pass

        def runtest(self):
            # satisfy `FixtureRequest` constructor...
            self.funcargs = {}
            fixture_request = doctest_plugin._setup_fixtures(self)

            failed, tot = doctest.testfile(
                str(self.fspath), module_relative=False,
                optionflags=opts, parser=DocTestParserPlus(),
                extraglobs=dict(getfixture=fixture_request.getfuncargvalue),
                raise_on_error=True, verbose=False, encoding='utf-8')

        def reportinfo(self):
            """
            Overwrite pytest's ``DoctestItem`` because
            ``DocTestTextfilePlus`` does not have a ``dtest`` attribute
            which is used by pytest>=3.2.0 to return the location of the
            tests.

            For details see `pytest-dev/pytest#2651
            <https://github.com/pytest-dev/pytest/pull/2651>`_.
            """
            return self.fspath, None, "[doctest] %s" % self.name

    class DocTestParserPlus(doctest.DocTestParser):
        """
        An extension to the builtin DocTestParser that handles the
        special directives for skipping tests.

        The directives are:

           - ``.. doctest-skip::``: Skip the next doctest chunk.

           - ``.. doctest-requires:: module1, module2``: Skip the next
             doctest chunk if the given modules/packages are not
             installed.

           - ``.. doctest-skip-all``: Skip all subsequent doctests.
        """

        def parse(self, s, name=None):
            result = doctest.DocTestParser.parse(self, s, name=name)

            # result is a sequence of alternating text chunks and
            # doctest.Example objects.  We need to look in the text
            # chunks for the special directives that help us determine
            # whether the following examples should be skipped.

            required = []
            skip_next = False
            skip_all = False

            for entry in result:
                if isinstance(entry, six.string_types) and entry:
                    required = []
                    skip_next = False
                    lines = entry.strip().splitlines()

                    if '.. doctest-skip-all' in (x.strip() for x in lines):
                        skip_all = True
                        continue

                    if not len(lines):
                        continue

                    last_line = lines[-1]
                    match = re.match(
                        r'\.\.\s+doctest-skip\s*::(\s+.*)?', last_line)
                    if match:
                        marker = match.group(1)
                        if (marker is None or
                                (marker.strip() == 'win32' and
                                 sys.platform == 'win32')):
                            skip_next = True
                            continue

                    match = re.match(
                        r'\.\.\s+doctest-requires\s*::\s+(.*)',
                        last_line)
                    if match:
                        required = re.split(r'\s*,?\s*', match.group(1))
                elif isinstance(entry, doctest.Example):
                    if (skip_all or skip_next or
                        not DocTestFinderPlus.check_required_modules(required)):
                        entry.options[doctest.SKIP] = True

                    if (config.getoption('remote_data', 'none') != 'any' and
                        entry.options.get(REMOTE_DATA)):
                        entry.options[doctest.SKIP] = True

            return result

    config.pluginmanager.register(
        DoctestPlus(DocTestModulePlus, DocTestTextfilePlus,
                    config.getini('doctest_rst') or config.option.doctest_rst),
        'doctestplus')

    # Remove the doctest_plugin, or we'll end up testing the .rst files twice.
    config.pluginmanager.unregister(doctest_plugin)


class DoctestPlus(object):
    def __init__(self, doctest_module_item_cls, doctest_textfile_item_cls,
                 run_rst_doctests):
        """
        doctest_module_item_cls should be a class inheriting
        `pytest.doctest.DoctestItem` and `pytest.File`.  This class handles
        running of a single doctest found in a Python module.  This is passed
        in as an argument because the actual class to be used may not be
        available at import time, depending on whether or not the doctest
        plugin for py.test is available.
        """
        self._doctest_module_item_cls = doctest_module_item_cls
        self._doctest_textfile_item_cls = doctest_textfile_item_cls
        self._run_rst_doctests = run_rst_doctests

        # Directories to ignore when adding doctests
        self._ignore_paths = []

    def pytest_ignore_collect(self, path, config):
        """Skip paths that match any of the doctest_norecursedirs patterns."""

        for pattern in config.getini("doctest_norecursedirs"):
            if path.check(fnmatch=pattern):
                # Apparently pytest_ignore_collect causes files not to be
                # collected by any test runner; for DoctestPlus we only want to
                # avoid creating doctest nodes for them
                self._ignore_paths.append(path)
                break

        return False

    def pytest_collect_file(self, path, parent):
        """Implements an enhanced version of the doctest module from py.test
        (specifically, as enabled by the --doctest-modules option) which
        supports skipping all doctests in a specific docstring by way of a
        special ``__doctest_skip__`` module-level variable.  It can also skip
        tests that have special requirements by way of
        ``__doctest_requires__``.

        ``__doctest_skip__`` should be a list of functions, classes, or class
        methods whose docstrings should be ignored when collecting doctests.

        This also supports wildcard patterns.  For example, to run doctests in
        a class's docstring, but skip all doctests in its modules use, at the
        module level::

            __doctest_skip__ = ['ClassName.*']

        You may also use the string ``'.'`` in ``__doctest_skip__`` to refer
        to the module itself, in case its module-level docstring contains
        doctests.

        ``__doctest_requires__`` should be a dictionary mapping wildcard
        patterns (in the same format as ``__doctest_skip__``) to a list of one
        or more modules that should be *importable* in order for the tests to
        run.  For example, if some tests require the scipy module to work they
        will be skipped unless ``import scipy`` is possible.  It is also
        possible to use a tuple of wildcard patterns as a key in this dict::

            __doctest_requires__ = {('func1', 'func2'): ['scipy']}

        """

        for ignore_path in self._ignore_paths:
            if ignore_path.common(path) == ignore_path:
                return None

        if path.ext == '.py':
            if path.basename == 'conf.py':
                return None

            # Don't override the built-in doctest plugin
            return self._doctest_module_item_cls(path, parent)
        elif self._run_rst_doctests and path.ext == '.rst':
            # Ignore generated .rst files
            parts = str(path).split(os.path.sep)

            # Don't test files that start with a _
            if path.basename.startswith('_'):
                return None

            # Don't test files in directories that start with a '_' if those
            # directories are inside docs. Note that we *should* allow for
            # example /tmp/_q/docs/file.rst but not /tmp/docs/_build/file.rst
            # If we don't find 'docs' in the path, we should just skip this
            # check to be safe. We also want to skip any api sub-directory
            # of docs.
            if 'docs' in parts:
                # We index from the end on the off chance that the temporary
                # directory includes 'docs' in the path, e.g.
                # /tmp/docs/371j/docs/index.rst You laugh, but who knows! :)
                # Also, it turns out lists don't have an rindex method. Huh??!!
                docs_index = len(parts) - 1 - parts[::-1].index('docs')
                if any(x.startswith('_') or x == 'api' for x in parts[docs_index:]):
                    return None

            # TODO: Get better names on these items when they are
            # displayed in py.test output
            return self._doctest_textfile_item_cls(path, parent)


class DocTestFinderPlus(doctest.DocTestFinder):
    """Extension to the default `doctest.DoctestFinder` that supports
    ``__doctest_skip__`` magic.  See `pytest_collect_file` for more details.
    """

    # Caches the results of import attempts
    _import_cache = {}

    @classmethod
    def check_required_modules(cls, mods):
        for mod in mods:
            if mod in cls._import_cache:
                if not cls._import_cache[mod]:
                    return False
            try:
                imp.find_module(mod)
            except ImportError:
                cls._import_cache[mod] = False
                return False
            else:
                cls._import_cache[mod] = True
        return True

    def find(self, obj, name=None, module=None, globs=None,
             extraglobs=None):
        tests = doctest.DocTestFinder.find(self, obj, name, module, globs,
                                           extraglobs)
        if (hasattr(obj, '__doctest_skip__') or
                hasattr(obj, '__doctest_requires__')):
            if name is None and hasattr(obj, '__name__'):
                name = obj.__name__
            else:
                raise ValueError("DocTestFinder.find: name must be given "
                                 "when obj.__name__ doesn't exist: {!r}".format(
                        (type(obj),)))

            def test_filter(test):
                for pat in getattr(obj, '__doctest_skip__', []):
                    if pat == '*':
                        return False
                    elif pat == '.' and test.name == name:
                        return False
                    elif fnmatch.fnmatch(test.name, '.'.join((name, pat))):
                        return False

                reqs = getattr(obj, '__doctest_requires__', {})
                for pats, mods in six.iteritems(reqs):
                    if not isinstance(pats, tuple):
                        pats = (pats,)
                    for pat in pats:
                        if not fnmatch.fnmatch(test.name,
                                               '.'.join((name, pat))):
                            continue
                        if not self.check_required_modules(mods):
                            return False
                return True

            tests = list(filter(test_filter, tests))

        return tests

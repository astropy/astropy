# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
These plugins modify the behavior of py.test and are meant to be imported
into conftest.py in the root directory.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import __future__

from ..extern import six
from ..extern.six.moves import filter

import ast
import doctest
import fnmatch
import imp
import io
import locale
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import types

from .helper import (
    pytest, treat_deprecations_as_exceptions, enable_deprecations_as_exceptions)
from .disable_internet import turn_off_internet, turn_on_internet
from .output_checker import AstropyOutputChecker, FIX, FLOAT_CMP

# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.


def pytest_addoption(parser):
    parser.addoption("--remote-data", action="store_true",
                     help="run tests with online data")
    parser.addoption("--open-files", action="store_true",
                     help="fail if any test leaves files open")

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
doctest.OutputChecker = AstropyOutputChecker


REMOTE_DATA = doctest.register_optionflag('REMOTE_DATA')


def pytest_configure(config):
    treat_deprecations_as_exceptions()

    # Monkeypatch to deny access to remote resources unless explicitly told
    # otherwise
    if not config.getoption('remote_data'):
        turn_off_internet(verbose=config.option.verbose)

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
            if self.fspath.basename == "conftest.py":
                module = self.config._conftest.importconftest(self.fspath)
            else:
                try:
                    module = self.fspath.pyimport()
                    # Just ignore searching modules that can't be imported when
                    # collecting doctests
                except ImportError:
                    raise StopIteration

            # uses internal doctest module parsing mechanism
            finder = DocTestFinderPlus()
            runner = doctest.DebugRunner(verbose=False, optionflags=opts)
            for test in finder.find(module):
                if test.examples:  # skip empty doctests
                    if not config.getvalue("remote_data"):
                        for example in test.examples:
                            if example.options.get(REMOTE_DATA):
                                example.options[doctest.SKIP] = True

                    yield doctest_plugin.DoctestItem(
                        test.name, self, runner, test)

        # This is for py.test prior to 2.4.0
        def runtest(self):
            return

    class DocTestTextfilePlus(doctest_plugin.DoctestTextfile):
        def runtest(self):
            # satisfy `FixtureRequest` constructor...
            self.funcargs = {}
            self._fixtureinfo = doctest_plugin.FuncFixtureInfo((), [], {})
            fixture_request = doctest_plugin.FixtureRequest(self)
            failed, tot = doctest.testfile(
                str(self.fspath), module_relative=False,
                optionflags=opts, parser=DocTestParserPlus(),
                extraglobs=dict(getfixture=fixture_request.getfuncargvalue),
                raise_on_error=True, verbose=False, encoding='utf-8')

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

                    if (not config.getvalue('remote_data') and
                        entry.options.get(REMOTE_DATA)):
                        entry.options[doctest.SKIP] = True

            return result

    config.pluginmanager.register(
        DoctestPlus(DocTestModulePlus, DocTestTextfilePlus,
                    config.getini('doctest_rst') or config.option.doctest_rst),
        'doctestplus')

    # Remove the doctest_plugin, or we'll end up testing the .rst
    # files twice.
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

        if run_rst_doctests and six.PY3:
            self._run_rst_doctests = False

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
            parts = bytes(path).split(os.path.sep)
            if (path.basename.startswith(b'_') or
                any(x.startswith(b'_') for x in parts) or
                any(x == b'api' for x in parts)):
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
                return cls._import_cache[mod]
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
                                 "when obj.__name__ doesn't exist: %r" %
                                 (type(obj),))

            def test_filter(test):
                for pat in getattr(obj, '__doctest_skip__', []):
                    if pat == '*':
                        return False
                    elif pat == '.' and test.name == name:
                        return False
                    elif fnmatch.fnmatch(test.name, '.'.join((name, pat))):
                        return False

                reqs = getattr(obj, '__doctest_requires__', {})
                for pats, mods in list(six.iteritems(reqs)):
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


# Open file detection.
#
# This works by calling out to lsof to get the list of open files held
# by the process both before and after the test.  If something is
# still open after the test that wasn't open before the test, an
# AssertionError is raised.
#
# This is not thread-safe.  We're not currently running our tests
# multi-threaded, but that is worth noting.

SUPPORTS_OPEN_FILE_DETECTION = (
    sys.platform in ('linux', 'linux2', 'darwin'))


def _get_open_file_list():
    fsencoding = sys.getfilesystemencoding()

    sproc = subprocess.Popen(
        ['lsof -F0 -n -p {0}'.format(os.getpid())],
        shell=True, stdout=subprocess.PIPE)
    output = sproc.communicate()[0].strip()
    files = []
    for line in output.split(b'\n'):
        columns = line.split(b'\0')
        mapping = {}
        for column in columns:
            if len(column) >= 2:
                mapping[column[0:1]] = column[1:]

        if (mapping.get(b'f') and
            mapping.get(b'a', b' ') != b' ' and
                mapping.get(b't') == b'REG'):
            # Ignore extension modules -- they may be imported by a
            # test but are never again closed by the runtime.  That's
            # ok.
            for suffix, mode, filetype in imp.get_suffixes():
                if mapping[b'n'].decode(fsencoding).endswith(suffix):
                    break
            else:
                files.append(mapping[b'n'])

    return set(files)


def pytest_runtest_setup(item):
    # Store a list of the currently opened files so we can compare
    # against them when the test is done.
    if SUPPORTS_OPEN_FILE_DETECTION and item.config.getvalue('open_files'):
        item.open_files = _get_open_file_list()

    if ('remote_data' in item.keywords and
            not item.config.getvalue("remote_data")):
        pytest.skip("need --remote-data option to run")


if SUPPORTS_OPEN_FILE_DETECTION:
    def pytest_runtest_teardown(item, nextitem):
        # a "skipped" test will not have been called with
        # pytest_runtest_setup, so therefore won't have an
        # "open_files" member
        if (not item.config.getvalue('open_files') or
                not hasattr(item, 'open_files')):
            return

        start_open_files = item.open_files
        del item.open_files

        open_files = _get_open_file_list()

        # This works in tandem with the test_open_file_detection test to
        # ensure that it creates one extra open file.
        if item.name == 'test_open_file_detection':
            assert len(start_open_files) + 1 == len(open_files)
            return

        not_closed = set()
        for filename in open_files:
            # astropy.log files are allowed to continue to exist
            # between test runs
            if os.path.basename(filename) == 'astropy.log':
                continue

            if filename not in start_open_files:
                not_closed.add(filename)

        if len(not_closed):
            msg = ['File(s) not closed:']
            for name in not_closed:
                msg.append('  {0}'.format(
                    name.decode(sys.getfilesystemencoding())))
            raise AssertionError('\n'.join(msg))


def pytest_report_header(config):
    from .. import __version__

    stdoutencoding = getattr(sys.stdout, 'encoding') or 'ascii'

    s = "\nRunning tests with Astropy version {0}.\n".format(__version__)
    if six.PY2:
        args = [x.decode('utf-8') for x in config.args]
    elif six.PY3:
        args = config.args
    s += "Running tests in {0}.\n\n".format(" ".join(args))

    from platform import platform
    plat = platform()
    if isinstance(plat, bytes):
        plat = plat.decode(stdoutencoding, 'replace')
    s += "Platform: {0}\n\n".format(plat)
    s += "Executable: {0}\n\n".format(sys.executable)
    s += "Full Python Version: \n{0}\n\n".format(sys.version)

    s += "encodings: sys: {0}, locale: {1}, filesystem: {2}".format(
        sys.getdefaultencoding(),
        locale.getpreferredencoding(),
        sys.getfilesystemencoding())
    if sys.version_info < (3, 3, 0):
        s += ", unicode bits: {0}".format(
            int(math.log(sys.maxunicode, 2)))
    s += '\n'

    s += "byteorder: {0}\n".format(sys.byteorder)
    s += "float info: dig: {0.dig}, mant_dig: {0.dig}\n\n".format(
        sys.float_info)

    import numpy
    s += "Numpy: {0}\n".format(numpy.__version__)

    try:
        import scipy
        s += "Scipy: {0}\n".format(scipy.__version__)
    except:
        s += "Scipy: not available\n"

    try:
        import matplotlib
        s += "Matplotlib: {0}\n".format(matplotlib.__version__)
    except:
        s += "Matplotlib: not available\n"

    try:
        import h5py.version
        s += "h5py: {0}\n".format(h5py.version.version)
    except:
        s += "h5py: not available\n"

    special_opts = ["remote_data", "pep8"]
    opts = []
    for op in special_opts:
        if getattr(config.option, op, None):
            opts.append(op)
    if opts:
        s += "Using Astropy options: {0}.\n".format(" ".join(opts))

    if six.PY3 and (config.getini('doctest_rst') or config.option.doctest_rst):
        s += "Running doctests in .rst files is not supported on Python 3.x\n"

    if not six.PY3:
        s = s.encode(stdoutencoding, 'replace')

    return s


def pytest_pycollect_makemodule(path, parent):
    # This is where we set up testing both with and without
    # from __future__ import unicode_literals

    # On Python 3, just do the regular thing that py.test does
    if six.PY3:
        return pytest.Module(path, parent)
    elif six.PY2:
        return Pair(path, parent)


class Pair(pytest.File):
    """
    This class treats a given test .py file as a pair of .py files
    where one has __future__ unicode_literals and the other does not.
    """
    def collect(self):
        # First, just do the regular import of the module to make
        # sure it's sane and valid.  This block is copied directly
        # from py.test
        try:
            mod = self.fspath.pyimport(ensuresyspath=True)
        except SyntaxError:
            import py
            excinfo = py.code.ExceptionInfo()
            raise self.CollectError(excinfo.getrepr(style="short"))
        except self.fspath.ImportMismatchError:
            e = sys.exc_info()[1]
            raise self.CollectError(
                "import file mismatch:\n"
                "imported module %r has this __file__ attribute:\n"
                "  %s\n"
                "which is not the same as the test file we want to collect:\n"
                "  %s\n"
                "HINT: remove __pycache__ / .pyc files and/or use a "
                "unique basename for your test file modules"
                % e.args
            )

        # Now get the file's content.
        with io.open(six.text_type(self.fspath), 'rb') as fd:
            content = fd.read()

        # If the file contains the special marker, only test it both ways.
        if b'TEST_UNICODE_LITERALS' in content:
            # Return the file in both unicode_literal-enabled and disabled forms
            return [
                UnicodeLiteralsModule(mod.__name__, content, self.fspath, self),
                NoUnicodeLiteralsModule(mod.__name__, content, self.fspath, self)
            ]
        else:
            return [pytest.Module(self.fspath, self)]


_RE_FUTURE_IMPORTS = re.compile(br'from __future__ import ((\(.*?\))|([^\n]+))',
                                flags=re.DOTALL)


class ModifiedModule(pytest.Module):
    def __init__(self, mod_name, content, path, parent):
        self.mod_name = mod_name
        self.content = content
        super(ModifiedModule, self).__init__(path, parent)

    def _importtestmodule(self):
        # We have to remove the __future__ statements *before* parsing
        # with compile, otherwise the flags are ignored.
        content = re.sub(_RE_FUTURE_IMPORTS, b'', self.content)

        new_mod = types.ModuleType(self.mod_name)
        new_mod.__file__ = six.text_type(self.fspath)

        if hasattr(self, '_transform_ast'):
            # ast.parse doesn't let us hand-select the __future__
            # statements, but built-in compile, with the PyCF_ONLY_AST
            # flag does.
            tree = compile(
                content, six.text_type(self.fspath), 'exec',
                self.flags | ast.PyCF_ONLY_AST, True)
            tree = self._transform_ast(tree)
            # Now that we've transformed the tree, recompile it
            code = compile(
                tree, six.text_type(self.fspath), 'exec')
        else:
            # If we don't need to transform the AST, we can skip
            # parsing/compiling in two steps
            code = compile(
                content, six.text_type(self.fspath), 'exec',
                self.flags, True)

        pwd = os.getcwd()
        try:
            os.chdir(os.path.dirname(six.text_type(self.fspath)))
            six.exec_(code, new_mod.__dict__)
        finally:
            os.chdir(pwd)
        self.config.pluginmanager.consider_module(new_mod)
        return new_mod


class UnicodeLiteralsModule(ModifiedModule):
    flags = (
        __future__.absolute_import.compiler_flag |
        __future__.division.compiler_flag |
        __future__.print_function.compiler_flag |
        __future__.unicode_literals.compiler_flag
    )


class NoUnicodeLiteralsModule(ModifiedModule):
    flags = (
        __future__.absolute_import.compiler_flag |
        __future__.division.compiler_flag |
        __future__.print_function.compiler_flag
    )

    def _transform_ast(self, tree):
        # When unicode_literals is disabled, we still need to convert any
        # byte string containing non-ascii characters into a Unicode string.
        # If it doesn't decode as utf-8, we assume it's some other kind
        # of byte string and just ultimately leave it alone.

        # Note that once we drop support for Python 3.2, we should be
        # able to remove this transformation and just put explicit u''
        # prefixes in the test source code.

        class NonAsciiLiteral(ast.NodeTransformer):
            def visit_Str(self, node):
                s = node.s
                if isinstance(s, bytes):
                    try:
                        s.decode('ascii')
                    except UnicodeDecodeError:
                        try:
                            s = s.decode('utf-8')
                        except UnicodeDecodeError:
                            pass
                        else:
                            return ast.copy_location(ast.Str(s=s), node)
                return node
        return NonAsciiLiteral().visit(tree)

def pytest_unconfigure():
    """
    Cleanup post-testing
    """
    # restore internet connectivity (only lost if remote_data=False and
    # turn_off_internet previously called)
    # this is harmless / does nothing if socket connections were never disabled
    turn_on_internet()

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
These plugins modify the behavior of py.test and are meant to be imported
into conftest.py in the root directory.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import __future__

from ..extern import six
from ..extern.six.moves import builtins

import ast
import datetime
import io
import locale
import math
import os
import re
import sys
import types
from pkgutil import find_loader
from collections import OrderedDict

import pytest

from ..config.paths import set_temp_config, set_temp_cache
from .helper import treat_deprecations_as_exceptions, ignore_warnings
from .helper import enable_deprecations_as_exceptions  # pylint: disable=W0611
from ..utils.argparse import writeable_directory
from ..utils.introspection import resolve_name

try:
    import importlib.machinery as importlib_machinery
except ImportError:  # Python 2.7
    importlib_machinery = None

pytest_plugins = ['astropy.tests.pytest_repeat']

_PLUGINS_PREFIX = 'astropy.extern.plugins'
for plugin in ['pytest_doctestplus', 'pytest_openfiles', 'pytest_remotedata']:
    if find_loader(plugin) is None:
        pytest_plugins.append('{}.{}.plugin'.format(_PLUGINS_PREFIX, plugin))

# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.

# This makes sure that this module is not collected when running the test
# suite. This is necessary in order to get the test suite to run without errors
# using pytest>=3.7
if getattr(builtins, '_pytest_running', False):
    pytest.skip()


def pytest_addoption(parser):

    parser.addoption("--astropy-config-dir", nargs='?', type=writeable_directory,
                     help="specify directory for storing and retrieving the "
                          "Astropy configuration during tests (default is "
                          "to use a temporary directory created by the test "
                          "runner); be aware that using an Astropy config "
                          "file other than the default can cause some tests "
                          "to fail unexpectedly")

    parser.addoption("--astropy-cache-dir", nargs='?', type=writeable_directory,
                     help="specify directory for storing and retrieving the "
                          "Astropy cache during tests (default is "
                          "to use a temporary directory created by the test "
                          "runner)")
    parser.addini("astropy_config_dir",
                  "specify directory for storing and retrieving the "
                  "Astropy configuration during tests (default is "
                  "to use a temporary directory created by the test "
                  "runner); be aware that using an Astropy config "
                  "file other than the default can cause some tests "
                  "to fail unexpectedly", default=None)

    parser.addini("astropy_cache_dir",
                  "specify directory for storing and retrieving the "
                  "Astropy cache during tests (default is "
                  "to use a temporary directory created by the test "
                  "runner)", default=None)


def pytest_configure(config):
    treat_deprecations_as_exceptions()

def pytest_runtest_setup(item):
    config_dir = item.config.getini('astropy_config_dir')
    cache_dir = item.config.getini('astropy_cache_dir')

    # Command-line options can override, however
    config_dir = item.config.getoption('astropy_config_dir') or config_dir
    cache_dir = item.config.getoption('astropy_cache_dir') or cache_dir

    # We can't really use context managers directly in py.test (although
    # py.test 2.7 adds the capability), so this may look a bit hacky
    if config_dir:
        item.set_temp_config = set_temp_config(config_dir)
        item.set_temp_config.__enter__()
    if cache_dir:
        item.set_temp_cache = set_temp_cache(cache_dir)
        item.set_temp_cache.__enter__()



def pytest_runtest_teardown(item, nextitem):
    if hasattr(item, 'set_temp_cache'):
        item.set_temp_cache.__exit__()
    if hasattr(item, 'set_temp_config'):
        item.set_temp_config.__exit__()


PYTEST_HEADER_MODULES = OrderedDict([('Numpy', 'numpy'),
                                     ('Scipy', 'scipy'),
                                     ('Matplotlib', 'matplotlib'),
                                     ('h5py', 'h5py'),
                                     ('Pandas', 'pandas')])

# This always returns with Astropy's version
from .. import __version__

TESTED_VERSIONS = OrderedDict([('Astropy', __version__)])


def pytest_report_header(config):

    try:
        stdoutencoding = sys.stdout.encoding or 'ascii'
    except AttributeError:
        stdoutencoding = 'ascii'

    if six.PY2:
        args = [x.decode('utf-8') for x in config.args]
    else:
        args = config.args

    # TESTED_VERSIONS can contain the affiliated package version, too
    if len(TESTED_VERSIONS) > 1:
        for pkg, version in TESTED_VERSIONS.items():
            if pkg not in ['Astropy', 'astropy_helpers']:
                s = "\nRunning tests with {0} version {1}.\n".format(
                    pkg, version)
    else:
        s = "\nRunning tests with Astropy version {0}.\n".format(
            TESTED_VERSIONS['Astropy'])

    # Per https://github.com/astropy/astropy/pull/4204, strip the rootdir from
    # each directory argument
    if hasattr(config, 'rootdir'):
        rootdir = str(config.rootdir)
        if not rootdir.endswith(os.sep):
            rootdir += os.sep

        dirs = [arg[len(rootdir):] if arg.startswith(rootdir) else arg
                for arg in args]
    else:
        dirs = args

    s += "Running tests in {0}.\n\n".format(" ".join(dirs))

    s += "Date: {0}\n\n".format(datetime.datetime.now().isoformat()[:19])

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

    for module_display, module_name in six.iteritems(PYTEST_HEADER_MODULES):
        try:
            with ignore_warnings(DeprecationWarning):
                module = resolve_name(module_name)
        except ImportError:
            s += "{0}: not available\n".format(module_display)
        else:
            try:
                version = module.__version__
            except AttributeError:
                version = 'unknown (no __version__ attribute)'
            s += "{0}: {1}\n".format(module_display, version)

    # Helpers version
    if 'astropy_helpers' in TESTED_VERSIONS:
        astropy_helpers_version = TESTED_VERSIONS['astropy_helpers']
    else:
        try:
            from ..version import astropy_helpers_version
        except ImportError:
            astropy_helpers_version = None

    if astropy_helpers_version:
        s += "astropy_helpers: {0}\n".format(astropy_helpers_version)

    special_opts = ["remote_data", "pep8"]
    opts = []
    for op in special_opts:
        op_value = getattr(config.option, op, None)
        if op_value:
            if isinstance(op_value, six.string_types):
                op = ': '.join((op, op_value))
            opts.append(op)
    if opts:
        s += "Using Astropy options: {0}.\n".format(", ".join(opts))

    if six.PY2:
        s = s.encode(stdoutencoding, 'replace')

    return s


def pytest_pycollect_makemodule(path, parent):
    # This is where we set up testing both with and without
    # from __future__ import unicode_literals

    # On Python 3, just do the regular thing that py.test does
    if six.PY2:
        return Pair(path, parent)
    else:
        return pytest.Module(path, parent)


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
                "imported module {!r} has this __file__ attribute:\n"
                "  {}\n"
                "which is not the same as the test file we want to collect:\n"
                "  {}\n"
                "HINT: remove __pycache__ / .pyc files and/or use a "
                "unique basename for your test file modules".format(e.args))

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
        content = re.sub(_RE_FUTURE_IMPORTS, b'\n', self.content)

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


def pytest_terminal_summary(terminalreporter):
    """Output a warning to IPython users in case any tests failed."""

    try:
        get_ipython()
    except NameError:
        return

    if not terminalreporter.stats.get('failed'):
        # Only issue the warning when there are actually failures
        return

    terminalreporter.ensure_newline()
    terminalreporter.write_line(
        'Some tests are known to fail when run from the IPython prompt; '
        'especially, but not limited to tests involving logging and warning '
        'handling.  Unless you are certain as to the cause of the failure, '
        'please check that the failure occurs outside IPython as well.  See '
        'http://docs.astropy.org/en/stable/known_issues.html#failing-logging-'
        'tests-when-running-the-tests-in-ipython for more information.',
        yellow=True, bold=True)

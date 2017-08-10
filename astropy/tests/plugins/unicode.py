# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This plugin provide support for unicode literals in Python 2.7.  
"""
#TODO: this plugin is necessary only for as long as we support Python 2
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import __future__


import ast
import io
import os
import re
import sys
import types

import pytest

from ...extern import six


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

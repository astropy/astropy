# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib
import secrets
from textwrap import dedent

import pytest

from astropy.utils.parsing import _TAB_HEADER


def _docstring_canary():
    """Docstring that's here just to check for -OO."""


@pytest.mark.skipif(not _docstring_canary.__doc__, reason="Test cannot be run with -OO")
def test_generate_parser(tmp_path, monkeypatch):
    # Write Python code into the temporary directory, so that the
    # generated tables will also go into the temporary directory.
    # We use a unique suffix so that the test can be run multiple times
    # without weirdness due to module caching.
    suffix = secrets.token_hex(16)
    lexer_file = tmp_path / f"test_parsing_lexer_{suffix}.py"
    lexer_file.write_text(
        dedent(
            rf"""
        from astropy.utils.parsing import lex

        def make_lexer():
            tokens = ('NUMBER', 'PLUS')
            t_PLUS = r'\+'

            def t_NUMBER(t):
                r'\d+'
                t.value = int(t.value)
                return t

            return lex('test_parsing_lextab_{suffix}', 'test_parsing_lexer_{suffix}')
        """
        )
    )
    parser_file = tmp_path / f"test_parsing_parser_{suffix}.py"
    parser_file.write_text(
        dedent(
            rf"""
        from astropy.utils.parsing import yacc

        def make_parser():
            tokens = ('NUMBER', 'PLUS')

            def p_expression_number(p):
                'expression : NUMBER'
                p[0] = p[1]

            def p_expression_plus(p):
                'expression : expression PLUS NUMBER'
                p[0] = p[1] + p[3]

            return yacc('test_parsing_parsetab_{suffix}', 'test_parsing_parser_{suffix}')
        """
        )
    )

    monkeypatch.syspath_prepend(tmp_path)

    lexer_mod = importlib.import_module(f"test_parsing_lexer_{suffix}")
    lexer = lexer_mod.make_lexer()
    parser_mod = importlib.import_module(f"test_parsing_parser_{suffix}")
    parser = parser_mod.make_parser()
    result = parser.parse("1+2+3", lexer=lexer)
    assert result == 6

    lextab = (tmp_path / f"test_parsing_lextab_{suffix}.py").read_text()
    assert lextab.startswith(_TAB_HEADER.format(package=f"test_parsing_lexer_{suffix}"))
    parsetab = (tmp_path / f"test_parsing_parsetab_{suffix}.py").read_text()
    assert parsetab.startswith(
        _TAB_HEADER.format(package=f"test_parsing_parser_{suffix}")
    )

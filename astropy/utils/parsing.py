# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Wrappers for PLY to provide thread safety.
"""

_TAB_HEADER = """# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest {package}
#
# You can then commit the changes to this file.

"""

import contextlib
import functools
import re
import os
import threading


_LOCK = threading.RLock()


def _add_tab_header(filename, package):
    with open(filename, 'r') as f:
        contents = f.read()

    with open(filename, 'w') as f:
        f.write(_TAB_HEADER.format(package=package))
        f.write(contents)


@contextlib.contextmanager
def _patch_get_caller_module_dict(module):
    """Temporarily replace the module's get_caller_module_dict.

    It's wrapped to offset the lookup depth by one level, so that we get the
    function which defines the symbols instead of the lex/yacc wrapper in
    this module.
    """
    original = module.get_caller_module_dict

    @functools.wraps(original)
    def wrapper(levels):
        # Add 2, not 1, because the wrapper itself adds another level
        return original(levels + 2)

    module.get_caller_module_dict = wrapper
    yield
    module.get_caller_module_dict = original


def lex(lextab, package, reflags=int(re.VERBOSE)):
    """Create a lexer from local variables.

    It automatically compiles the lexer in optimized mode, writing to
    `lextab` in the same directory as the calling file.

    This function is thread-safe. The returned lexer is *not* thread-safe, but
    if it is used exclusively with a single parser returned by :func:`yacc`
    then it will be safe.

    The `package` is inserted into the instructions in the generated file that
    explain which package to test to rebuild it.

    It is only intended to work with callbacks defined in a closure, not a class
    or module.
    """
    from astropy.extern.ply import lex

    caller_file = lex.get_caller_module_dict(2)['__file__']
    lextab_filename = os.path.join(os.path.dirname(caller_file), lextab + '.py')
    with _LOCK:
        lextab_exists = os.path.exists(lextab_filename)
        with _patch_get_caller_module_dict(lex):
            lexer = lex.lex(optimize=True, lextab=lextab,
                            outputdir=os.path.dirname(caller_file),
                            reflags=reflags)
        if not lextab_exists:
            _add_tab_header(lextab_filename, package)
        return lexer


class _ThreadSafeParser:
    def __init__(self, parser):
        self.parser = parser
        self._lock = threading.RLock()

    def parse(self, *args, **kwargs):
        with self._lock:
            return self.parser.parse(*args, **kwargs)


def yacc(tabmodule, package):
    """Create a parser from local variables.

    It automatically compiles the parser in optimized mode, writing to
    `tabmodule` in the same directory as the calling file.

    This function is thread-safe, and the returned parser is also thread-safe,
    provided that it does not share a lexer with any other parser.

    The `package` is inserted into the instructions in the generated file that
    explain which package to test to rebuild it.

    It is only intended to work with callbacks defined in a closure, not a class
    or module.
    """
    from astropy.extern.ply import yacc

    caller_file = yacc.get_caller_module_dict(2)['__file__']
    tab_filename = os.path.join(os.path.dirname(caller_file), tabmodule + '.py')
    with _LOCK:
        tab_exists = os.path.exists(tab_filename)
        with _patch_get_caller_module_dict(yacc):
            parser = yacc.yacc(tabmodule=tabmodule,
                               outputdir=os.path.dirname(caller_file),
                               debug=False, optimize=True, write_tables=True)
        if not tab_exists:
            _add_tab_header(tab_filename, package)

    return _ThreadSafeParser(parser)

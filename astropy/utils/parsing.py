# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Wrappers for PLY to provide thread safety.
"""

import contextlib
import functools
import os
import re
import threading

__all__ = ["lex", "ThreadSafeParser", "yacc"]


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
_LOCK = threading.RLock()


def _add_tab_header(filename, package):
    with open(filename) as f:
        contents = f.read()

    with open(filename, "w") as f:
        f.write(_TAB_HEADER.format(package=package))
        f.write(contents)


@contextlib.contextmanager
def _patch_get_caller_module_dict(module):
    """Temporarily replace the module's get_caller_module_dict.

    This is a function inside ``ply.lex`` and ``ply.yacc`` (each has a copy)
    that is used to retrieve the caller's local symbols. Here, we patch the
    function to instead retrieve the grandparent's local symbols to account
    for a wrapper layer.
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
    ``lextab`` in the same directory as the calling file.

    This function is thread-safe. The returned lexer is *not* thread-safe, but
    if it is used exclusively with a single parser returned by :func:`yacc`
    then it will be safe.

    It is only intended to work with lexers defined within the calling
    function, rather than at class or module scope.

    Parameters
    ----------
    lextab : str
        Name for the file to write with the generated tables, if it does not
        already exist (without ``.py`` suffix).
    package : str
        Name of a test package which should be run with pytest to regenerate
        the output file. This is inserted into a comment in the generated
        file.
    reflags : int
        Passed to ``ply.lex``.
    """
    from astropy.extern.ply import lex

    caller_file = lex.get_caller_module_dict(2)["__file__"]
    caller_dir = os.path.dirname(caller_file)
    lextab_filename_py = os.path.join(caller_dir, lextab + ".py")
    lextab_filename_pyc = os.path.join(caller_dir, lextab + ".pyc")
    with _LOCK:
        lextab_exists = os.path.exists(lextab_filename_py) or os.path.exists(
            lextab_filename_pyc
        )
        with _patch_get_caller_module_dict(lex):
            lexer = lex.lex(
                optimize=True,
                lextab=lextab,
                outputdir=caller_dir,
                reflags=reflags,
            )
        if not lextab_exists:
            _add_tab_header(lextab_filename_py, package)
        return lexer


class ThreadSafeParser:
    """Wrap a parser produced by ``ply.yacc.yacc``.

    It provides a :meth:`parse` method that is thread-safe.
    """

    def __init__(self, parser):
        self.parser = parser
        self._lock = threading.RLock()

    def parse(self, *args, **kwargs):
        """Run the wrapped parser, with a lock to ensure serialization."""
        with self._lock:
            return self.parser.parse(*args, **kwargs)


def yacc(tabmodule, package):
    """Create a parser from local variables.

    It automatically compiles the parser in optimized mode, writing to
    ``tabmodule`` in the same directory as the calling file.

    This function is thread-safe, and the returned parser is also thread-safe,
    provided that it does not share a lexer with any other parser.

    It is only intended to work with parsers defined within the calling
    function, rather than at class or module scope.

    Parameters
    ----------
    tabmodule : str
        Name for the file to write with the generated tables, if it does not
        already exist (without ``.py`` suffix).
    package : str
        Name of a test package which should be run with pytest to regenerate
        the output file. This is inserted into a comment in the generated
        file.
    """
    from astropy.extern.ply import yacc

    caller_file = yacc.get_caller_module_dict(2)["__file__"]
    caller_dir = os.path.dirname(caller_file)
    tab_filename_py = os.path.join(caller_dir, tabmodule + ".py")
    tab_filename_pyc = os.path.join(caller_dir, tabmodule + ".pyc")
    with _LOCK:
        tab_exists = os.path.exists(tab_filename_py) or os.path.exists(tab_filename_pyc)
        with _patch_get_caller_module_dict(yacc):
            parser = yacc.yacc(
                tabmodule=tabmodule,
                outputdir=caller_dir,
                debug=False,
                optimize=True,
                write_tables=True,
            )
        if not tab_exists:
            _add_tab_header(tab_filename_py, package)

    return ThreadSafeParser(parser)

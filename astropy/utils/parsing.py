# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Wrappers for PLY to provide thread safety.
"""

from __future__ import annotations

import contextlib
import functools
import re
import threading
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Generator
    from types import ModuleType

    from astropy.extern.ply.lex import Lexer
    from astropy.extern.ply.yacc import LRParser

__all__ = ["ThreadSafeParser", "lex", "yacc"]


_TAB_HEADER = """# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest {package}
#
# You can then commit the changes to this file.

"""
_LOCK = threading.RLock()


@contextlib.contextmanager
def _patch_ply_module(
    module: ModuleType, file: Path, package: str
) -> Generator[None, None, None]:
    """Temporarily replace the module's get_caller_module_dict.

    This is a function inside ``ply.lex`` and ``ply.yacc`` (each has a copy)
    that is used to retrieve the caller's local symbols. Here, we patch the
    function to instead retrieve the grandparent's local symbols to account
    for a wrapper layer.

    Additionally, a custom header is inserted into any files ``ply`` writes.
    """
    original = module.get_caller_module_dict

    @functools.wraps(original)
    def wrapper(levels):
        # Add 2, not 1, because the wrapper itself adds another level
        return original(levels + 2)

    file_exists = file.exists() or file.with_suffix(".pyc").exists()
    module.get_caller_module_dict = wrapper
    yield
    module.get_caller_module_dict = original
    if not file_exists:
        file.write_text(_TAB_HEADER.format(package=package) + file.read_text())


def lex(lextab: str, package: str, reflags: int = int(re.VERBOSE)) -> Lexer:
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

    caller_dir = Path(lex.get_caller_module_dict(2)["__file__"]).parent
    with _LOCK, _patch_ply_module(lex, caller_dir / (lextab + ".py"), package):
        return lex.lex(
            optimize=True, lextab=lextab, outputdir=caller_dir, reflags=reflags
        )


class ThreadSafeParser:
    """Wrap a parser produced by ``ply.yacc.yacc``.

    It provides a :meth:`parse` method that is thread-safe.
    """

    def __init__(self, parser: LRParser) -> None:
        self.parser = parser
        self._lock = threading.RLock()

    def parse(self, *args, **kwargs):
        """Run the wrapped parser, with a lock to ensure serialization."""
        with self._lock:
            return self.parser.parse(*args, **kwargs)


def yacc(tabmodule: str, package: str) -> ThreadSafeParser:
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

    caller_dir = Path(yacc.get_caller_module_dict(2)["__file__"]).parent
    with _LOCK, _patch_ply_module(yacc, caller_dir / (tabmodule + ".py"), package):
        parser = yacc.yacc(
            tabmodule=tabmodule,
            outputdir=caller_dir,
            debug=False,
            optimize=True,
            write_tables=True,
        )
    return ThreadSafeParser(parser)

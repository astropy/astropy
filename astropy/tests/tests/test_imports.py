# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pkgutil
import subprocess
import sys
from pathlib import Path
from textwrap import dedent
from types import ModuleType

import pytest

import astropy


def test_imports():
    """
    This just imports all modules in astropy, making sure they don't have any
    dependencies that sneak through
    """

    def onerror(name):
        # We should raise any legitimate error that occurred, but not
        # any warnings which happen to be caught because of our pytest
        # settings (e.g., DeprecationWarning).
        try:
            raise  # noqa: PLE0704
        except Warning:
            pass

    for imper, nm, ispkg in pkgutil.walk_packages(
        ["astropy"], "astropy.", onerror=onerror
    ):
        imper.find_spec(nm)


def test_toplevel_namespace():
    import astropy

    d = dir(astropy)
    assert "os" not in d
    assert "log" in d
    assert "test" in d
    assert "sys" not in d


def test_toplevel_lazy_imports():
    # Check that subpackages are loaded on demand.
    cmd = dedent("""
    import astropy, sys
    assert 'astropy.units' not in sys.modules
    astropy.units
    assert 'astropy.units' in sys.modules
    """)
    cp = subprocess.check_call([sys.executable, "-c", cmd])
    assert cp == 0


def test_toplevel_attribute_error():
    # Ensure that our __getattr__ does not leak an import error or so.
    with pytest.raises(AttributeError, match="module 'astropy' has no"):
        astropy.nonsense


def test_completeness_toplevel__all__():
    # Implicitly check that all items in __all__ exist.
    all_items = {getattr(astropy, attr) for attr in astropy.__all__}
    # Verify that the list of modules in __all__ is complete.
    module_names = {
        item.__name__.partition(".")[2]
        for item in all_items
        if isinstance(item, ModuleType)
    }
    module_dirs = {
        f.name
        for f in Path(astropy.__file__).parent.glob("[a-z]*")
        if f.is_dir() and f.name != "extern"
    }
    assert module_names == module_dirs

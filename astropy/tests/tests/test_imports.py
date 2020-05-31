# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pkgutil


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
            raise
        except Warning:
            pass

    for imper, nm, ispkg in pkgutil.walk_packages(['astropy'], 'astropy.',
                                                  onerror=onerror):
        imper.find_module(nm)


def test_toplevel_namespace():
    import astropy
    d = dir(astropy)
    assert 'os' not in d
    assert 'log' in d
    assert 'test' in d
    assert 'sys' not in d

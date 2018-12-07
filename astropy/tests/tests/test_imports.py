# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pkgutil
import os
import types


def test_imports():
    """
    This just imports all modules in astropy, making sure they don't have any
    dependencies that sneak through
    """

    from astropy.utils import find_current_module

    pkgornm = find_current_module(1).__name__.split('.')[0]

    if isinstance(pkgornm, str):
        package = pkgutil.get_loader(pkgornm).load_module(pkgornm)
    elif (isinstance(pkgornm, types.ModuleType) and
            '__init__' in pkgornm.__file__):
        package = pkgornm
    else:
        msg = 'test_imports is not determining a valid package/package name'
        raise TypeError(msg)

    if hasattr(package, '__path__'):
        pkgpath = package.__path__
    elif hasattr(package, '__file__'):
        pkgpath = os.path.split(package.__file__)[0]
    else:
        raise AttributeError('package to generate config items for does not '
                             'have __file__ or __path__')

    prefix = package.__name__ + '.'

    def onerror(name):
        # A legitimate error occurred in a module that wasn't excluded
        raise

    for imper, nm, ispkg in pkgutil.walk_packages(pkgpath, prefix,
                                                  onerror=onerror):
        imper.find_module(nm)


def test_toplevel_namespace():
    import astropy
    d = dir(astropy)
    assert 'os' not in d
    assert 'log' in d
    assert 'test' in d
    assert 'sys' not in d

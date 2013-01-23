# Licensed under a 3-clause BSD style license - see LICENSE.rst

def test_imports():
    """
    This just imports all modules in astropy, making sure they don't have any
    dependencies that sneak through
    """

    from os.path import split
    from types import ModuleType
    from pkgutil import get_loader, walk_packages

    from ...utils import find_current_module

    pkgornm = find_current_module(1).__name__.split('.')[0]

    if isinstance(pkgornm, basestring):
        package = get_loader(pkgornm).load_module(pkgornm)
    elif isinstance(pkgornm, ModuleType) and '__init__' in pkgornm.__file__:
        package = pkgornm
    else:
        msg = 'test_imports is not determining a valid package/package name'
        raise TypeError(msg)

    if hasattr(package, '__path__'):
        pkgpath = package.__path__
    elif hasattr(package, '__file__'):
        pkgpath = split(package.__file__)[0]
    else:
        raise AttributeError('package to generate config items for does not '
                             'have __file__ or __path__')

    for imper, nm, ispkg in walk_packages(pkgpath, package.__name__ + '.'):
        imper.find_module(nm)

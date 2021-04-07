# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Checks for optional dependencies using lazy import from
`PEP 562 <https://www.python.org/dev/peps/pep-0562/>`_.
"""
import importlib

# First, the top-level packages:
# TODO: This list is a duplicate of the dependencies in setup.cfg "all", but
# some of the package names are different from the pip-install name (e.g.,
# beautifulsoup4 -> bs4).
_optional_deps = ['bleach', 'bottleneck', 'bs4', 'bz2', 'h5py', 'html5lib',
                  'IPython', 'jplephem', 'lxml', 'matplotlib', 'mpmath',
                  'pandas', 'PIL', 'pkg_resources', 'pytz', 'scipy', 'skyfield',
                  'sortedcontainers', 'lzma', 'yaml']
_deps = {k: k for k in _optional_deps}

# Any subpackages that have different import behavior:
_deps['plt'] = 'matplotlib.pyplot'

__all__ = [f"HAS_{pkg}" for pkg in _deps]


def __getattr__(name):
    if name in __all__:
        module_name = name[4:]

        try:
            importlib.import_module(_deps[module_name])
        except (ImportError, ModuleNotFoundError):
            return False
        return True

    # I don't think we can ever get here, because if a name doesn't exist, it
    # should fail to import from this module. But adding this for safety!
    raise AttributeError(f"Module {__name__!r} has no attribute {name!r}. To "
                         "add support for checking if a package name is "
                         "available for import, add the package name")

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Checks for optional dependencies using lazy import from
`PEP 562 <https://www.python.org/dev/peps/pep-0562/>`_.
"""
import importlib


__all__ = [  # noqa
    'HAS_bleach', 'HAS_bottleneck', 'HAS_bs4', 'HAS_bz2', 'HAS_h5py',
    'HAS_html5lib', 'HAS_IPython', 'HAS_jplephem', 'HAS_lxml', 'HAS_matplotlib',
    'HAS_mpmath', 'HAS_pandas', 'HAS_PIL', 'HAS_pkg_resources', 'HAS_pytz',
    'HAS_scipy', 'HAS_skyfield', 'HAS_sortedcontainers', 'HAS_lzma', 'HAS_yaml'
]


def __getattr__(name):
    if name in __all__:
        module_name = name[4:]

        try:
            importlib.import_module(module_name)
        except (ImportError, ModuleNotFoundError):
            return False
        return True

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

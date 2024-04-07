# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Checks for optional dependencies using lazy import from
`PEP 562 <https://www.python.org/dev/peps/pep-0562/>`_.
"""

from importlib.util import find_spec

# First, the top-level packages:
# TODO: This list is a duplicate of the dependencies in pyproject.toml "all", but
# some of the package names are different from the pip-install name (e.g.,
# beautifulsoup4 -> bs4).
_optional_deps = [
    "asdf_astropy",
    "bleach",
    "bottleneck",
    "bs4",
    "bz2",
    "dask",
    "fsspec",
    "h5py",
    "html5lib",
    "IPython",
    "jplephem",
    "lxml",
    "matplotlib",
    "mpmath",
    "pandas",
    "PIL",
    "pytz",
    "s3fs",
    "scipy",
    "skyfield",
    "sortedcontainers",
    "lzma",
    "pyarrow",
    "pytest_mpl",
    "array_api_strict",
]
_deps = {k.upper(): k for k in _optional_deps}

# Any subpackages that have different import behavior:
_deps["PLT"] = "matplotlib"

__all__ = [f"HAS_{pkg}" for pkg in _deps]


def __getattr__(name):
    if name in __all__:
        return find_spec(_deps[name.removeprefix("HAS_")]) is not None

    raise AttributeError(f"Module {__name__!r} has no attribute {name!r}.")

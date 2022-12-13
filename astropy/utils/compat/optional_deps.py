# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Checks for optional dependencies using lazy import from
`PEP 562 <https://www.python.org/dev/peps/pep-0562/>`_.
"""
import importlib
import warnings

# First, the top-level packages:
# TODO: This list is a duplicate of the dependencies in setup.cfg "all", but
# some of the package names are different from the pip-install name (e.g.,
# beautifulsoup4 -> bs4).
_optional_deps = [
    "asdf",
    "asdf_astropy",
    "bleach",
    "bottleneck",
    "bs4",
    "bz2",
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
]
_formerly_optional_deps = ["yaml"]  # for backward compatibility
_deps = {k.upper(): k for k in _optional_deps + _formerly_optional_deps}

# Any subpackages that have different import behavior:
_deps["PLT"] = "matplotlib.pyplot"

__all__ = [f"HAS_{pkg}" for pkg in _deps]


def __getattr__(name):
    if name in __all__:
        module_name = name[4:]

        if module_name == "YAML":
            from astropy.utils.exceptions import AstropyDeprecationWarning

            warnings.warn(
                "PyYaml is now a strict dependency. HAS_YAML is deprecated as "
                "of v5.0 and will be removed in a subsequent version.",
                category=AstropyDeprecationWarning,
            )

        try:
            importlib.import_module(_deps[module_name])
        except (ImportError, ModuleNotFoundError):
            return False
        return True

    raise AttributeError(f"Module {__name__!r} has no attribute {name!r}.")

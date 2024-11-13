# Licensed under a 3-clause BSD style license - see LICENSE.rst


from pathlib import Path

from numpy import get_include as get_numpy_include
from setuptools import Extension

ROOT = Path(__file__).parent.relative_to(Path.cwd())


def get_extensions():
    sources = [ROOT / "_np_utils.pyx", ROOT / "_column_mixins.pyx"]
    include_dirs = [get_numpy_include()]

    exts = [
        Extension(
            name=f"astropy.table.{source.stem}",
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
            sources=[str(s) for s in sources],
            include_dirs=include_dirs,
        )
        for source in sources
    ]

    return exts

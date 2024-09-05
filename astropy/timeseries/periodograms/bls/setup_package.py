# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pathlib import Path

from numpy import get_include as get_numpy_include
from setuptools import Extension

BLS_ROOT = Path(__file__).parent.resolve().relative_to(Path.cwd())


def get_extensions():
    ext = Extension(
        "astropy.timeseries.periodograms.bls._impl",
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        sources=[
            str(BLS_ROOT / "bls.c"),
            str(BLS_ROOT / "_impl.pyx"),
        ],
        include_dirs=[get_numpy_include()],
    )
    return [ext]

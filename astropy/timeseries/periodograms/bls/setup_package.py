# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from numpy import get_include as get_numpy_include
from setuptools import Extension

BLS_ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    ext = Extension(
        "astropy.timeseries.periodograms.bls._impl",
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        sources=[
            os.path.join(BLS_ROOT, "bls.c"),
            os.path.join(BLS_ROOT, "_impl.pyx"),
        ],
        include_dirs=[get_numpy_include()],
    )
    return [ext]

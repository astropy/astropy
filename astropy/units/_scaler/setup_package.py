# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Copied from astropy/convolution/setup_package.py

import os

from numpy import get_include as get_numpy_include
from setuptools import Extension

C_SCALER_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_SCALER_PKGDIR, filename) for filename in ["src/scaler.c"]]
INCLUDE_FILES = [
    os.path.join(C_SCALER_PKGDIR, filename)
    for filename in ["src/scaler_limited_api_workarounds.h"]
]


def get_extensions():
    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations
    _scaler_ext = Extension(
        name="astropy.units._scaler.scaler",
        sources=SRC_FILES,
        depends=INCLUDE_FILES,
        include_dirs=[get_numpy_include()],
        language="c",
    )

    return [_scaler_ext]

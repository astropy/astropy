# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Copied from astropy/convolution/setup_package.py

from pathlib import Path

from numpy import get_include as get_numpy_include
from setuptools import Extension

C_TIME_PKGDIR = Path(__file__).parent.resolve().relative_to(Path.cwd())

SRC_FILES = [str(C_TIME_PKGDIR / "src" / "parse_times.c")]


def get_extensions():
    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations
    _time_ext = Extension(
        name="astropy.time._parse_times",
        sources=SRC_FILES,
        include_dirs=[get_numpy_include()],
        language="c",
    )

    return [_time_ext]

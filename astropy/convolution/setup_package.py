# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys

import numpy
from setuptools import Extension

C_CONVOLVE_PKGDIR = os.path.relpath(os.path.dirname(__file__))

extra_compile_args = ["-UNDEBUG"]
if not sys.platform.startswith("win"):
    extra_compile_args.append("-fPIC")


def get_extensions():
    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations
    sources = [
        os.path.join(C_CONVOLVE_PKGDIR, "_convolve.pyx"),
        os.path.join(C_CONVOLVE_PKGDIR, "src", "convolve.c"),
    ]
    _convolve_ext = Extension(
        name="astropy.convolution._convolve",
        extra_compile_args=extra_compile_args,
        include_dirs=[numpy.get_include()],
        sources=sources,
    )
    return [_convolve_ext]

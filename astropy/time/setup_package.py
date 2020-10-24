# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Copied from astropy/convolution/setup_package.py

import os
from setuptools import Extension

import numpy

C_TIME_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_TIME_PKGDIR, filename)
             for filename in ['src/parse_times.c']]


def get_extensions():
    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations
    _time_ext = Extension(name='astropy.time._parse_times',
                          sources=SRC_FILES,
                          include_dirs=[numpy.get_include()],
                          language='c')

    return [_time_ext]

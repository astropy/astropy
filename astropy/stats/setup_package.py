# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
from setuptools import Extension

import numpy

C_STATS_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_STATS_PKGDIR, '_fast_sigma_clipping.c')]


def get_extensions():
    _sigma_clip_ext = Extension(name='astropy.stats._fast_sigma_clipping', sources=SRC_FILES,
                                include_dirs=[numpy.get_include()],
                                language='c')

    return [_sigma_clip_ext]

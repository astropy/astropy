# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from setuptools import Extension

import numpy

SRCDIR = os.path.join(os.path.relpath(os.path.dirname(__file__)), 'src')

SRCFILES = ['wirth_select.c', 'compute_bounds.c', 'fast_sigma_clip.c']

SRCFILES = [os.path.join(SRCDIR, srcfile) for srcfile in SRCFILES]


def get_extensions():
    _sigma_clip_ext = Extension(name='astropy.stats._fast_sigma_clip', sources=SRCFILES,
                                include_dirs=[numpy.get_include()],
                                language='c')

    return [_sigma_clip_ext]

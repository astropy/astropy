# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import numpy
from setuptools import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))
SRCFILES = ["wirth_select.c", "compute_bounds.c", "fast_sigma_clip.c"]
SRCFILES = [os.path.join(ROOT, "src", srcfile) for srcfile in SRCFILES]


def get_extensions():
    _sigma_clip_ext = Extension(
        name="astropy.stats._fast_sigma_clip",
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        sources=SRCFILES,
        include_dirs=[numpy.get_include()],
        language="c",
    )
    _stats_ext = Extension(
        name="astropy.stats._stats",
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        sources=[os.path.join(ROOT, "_stats.pyx")],
        include_dirs=[numpy.get_include()],
    )

    return [_sigma_clip_ext, _stats_ext]

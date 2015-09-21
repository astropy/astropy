# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import

import os
from distutils.extension import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = [os.path.join(ROOT, "_np_utils.pyx")]
    include_dirs = ['numpy']
    libraries = []

    table_ext = Extension(
        name="astropy.table._np_utils",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,)

    return [table_ext]


def requires_2to3():
    return False

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from distutils.extension import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = [os.path.join(ROOT, "cyjoin.pyx")]
    include_dirs = ['numpy']
    libraries = []

    table_ext = Extension(
        name="astropy.table.cyjoin",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,)

    return [table_ext]

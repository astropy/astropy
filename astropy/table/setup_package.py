# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from distutils.extension import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = ["_np_utils.pyx", "_column_mixins.pyx"]
    include_dirs = ['numpy']

    exts = [
        Extension(name='astropy.table.' + os.path.splitext(source)[0],
                  sources=[os.path.join(ROOT, source)],
                  include_dirs=include_dirs)
        for source in sources
    ]

    return exts

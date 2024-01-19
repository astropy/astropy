# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from numpy import get_include as get_numpy_include
from setuptools import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = ["_np_utils.pyx", "_column_mixins.pyx"]
    include_dirs = [get_numpy_include()]

    exts = [
        Extension(
            name="astropy.table." + os.path.splitext(source)[0],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
            sources=[os.path.join(ROOT, source)],
            include_dirs=include_dirs,
        )
        for source in sources
    ]

    return exts

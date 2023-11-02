# Licensed under a 3-clause BSD style license

import os

import numpy
from setuptools import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = [
        os.path.join(ROOT, "cparser.pyx"),
        os.path.join(ROOT, "src", "tokenizer.c"),
    ]
    ascii_ext = Extension(
        name="astropy.io.ascii.cparser",
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        include_dirs=[numpy.get_include()],
        sources=sources,
    )
    return [ascii_ext]

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from os.path import join

import numpy
from setuptools import Extension

BLS_ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    ext = Extension(
        "astropy.timeseries.periodograms.bls._impl",
        sources=[
            join(BLS_ROOT, "bls.c"),
            join(BLS_ROOT, "_impl.pyx"),
        ],
        include_dirs=[numpy.get_include()],
    )
    return [ext]

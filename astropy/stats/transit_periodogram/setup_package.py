# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from os.path import join

from distutils.core import Extension


TRANSIT_PERIODOGRAM_ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    ext = Extension(
        "astropy.stats.transit_periodogram._impl",
        sources=[
            join(TRANSIT_PERIODOGRAM_ROOT, "transit_periodogram.c"),
            join(TRANSIT_PERIODOGRAM_ROOT, "_impl.pyx"),
        ],
        include_dirs=["numpy"],
    )
    return [ext]

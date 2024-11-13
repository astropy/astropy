# Licensed under a 3-clause BSD style license - see LICENSE.rst

from os.path import dirname, join, relpath

from setuptools import Extension

ASTROPY_UTILS_ROOT = dirname(__file__)


def get_extensions():
    return [
        Extension(
            "astropy.utils._compiler",
            [relpath(join(ASTROPY_UTILS_ROOT, "src", "compiler.c"))],
        )
    ]

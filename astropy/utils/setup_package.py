# Licensed under a 3-clause BSD style license - see LICENSE.rst

from setuptools import Extension
from os.path import dirname, join, abspath

ASTROPY_UTILS_ROOT = dirname(__file__)


def get_extensions():
    return [
        Extension('astropy.utils._compiler',
                  [abspath(join(ASTROPY_UTILS_ROOT, 'src', 'compiler.c'))])
    ]

# Licensed under a 3-clause BSD style license - see LICENSE.rst

from distutils.core import Extension
from os.path import dirname, join, relpath

ASTROPY_UTILS_ROOT = dirname(__file__)


def get_extensions():
    return [
        Extension('astropy.utils._compiler',
                  [relpath(join(ASTROPY_UTILS_ROOT, 'src', 'compiler.c'))])
    ]

from distutils.core import Extension
from os.path import dirname, join, relpath


def get_extensions():
    ROOT = dirname(__file__)

    return [
        Extension('astropy.utils._compiler',
                  [relpath(join(ROOT, 'src', 'compiler.c'))])
        ]

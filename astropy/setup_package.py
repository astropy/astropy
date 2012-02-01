from distutils.core import Extension
from os.path import dirname, join

def get_extensions():
    ROOT = dirname(__file__)

    return [
        Extension('astropy._compiler',
                  [join(ROOT, 'src', 'compiler.c')])
        ]

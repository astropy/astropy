from distutils.core import Extension
from os.path import dirname, join, relpath


def get_extensions():
    ROOT = dirname(__file__)

    return [
        Extension('astropy.utils._compiler',
                  [relpath(join(ROOT, 'src', 'compiler.c'))])
        ]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.utils.tests': ['data/*.dat']
        }

import os

from distutils.core import Extension
from glob import glob

from ... import setup_helpers

def get_extensions():
    return [
        Extension('astropy.io.fits.compression',
                  glob(os.path.join(os.path.dirname(__file__), 'src/*.c')),
                  include_dirs=[setup_helpers.get_numpy_include_path()],
                  extra_compile_args=['-Wno-unused-function',
                                      '-Wno-strict-prototypes'])
    ]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.io.fits.tests': ['data/*.fits']}

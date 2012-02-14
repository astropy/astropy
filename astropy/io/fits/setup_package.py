# Licensed under a 3-clause BSD style license - see PYFITS.rst

import os

from distutils.core import Extension
from glob import glob

from ... import setup_helpers


def get_extensions():
    if setup_helpers.get_compiler_option() != 'msvc':
        extra_compile_args = ['-Wno-unused-function',
                              '-Wno-strict-prototypes']
    else:
        extra_compile_args = []
    return [
        Extension(
            'astropy.io.fits.compression',
            [os.path.relpath(x) for x in
             glob(os.path.join(os.path.dirname(__file__), 'src/*.c'))],
            include_dirs=[setup_helpers.get_numpy_include_path()],
            extra_compile_args=extra_compile_args)
    ]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.io.fits.tests': ['data/*.fits']}


def get_legacy_alias():
    return setup_helpers.add_legacy_alias('pyfits', 'astropy.io.fits')

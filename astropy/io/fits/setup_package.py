# Licensed under a 3-clause BSD style license - see PYFITS.rst

import os

from distutils.core import Extension
from glob import glob

from astropy import setup_helpers


def get_extensions():
    if setup_helpers.get_compiler_option() != 'msvc':
        # All of these switches are to silence warnings from compiling CFITSIO
        extra_compile_args = ['-Wno-unused-variable', '-Wno-parentheses',
                              '-Wno-uninitialized', '-Wno-format',
                              '-Wno-strict-prototypes', '-Wno-unused',
                              '-Wno-comments', '-Wno-switch']
    else:
        extra_compile_args = []

    cfitsio_path = os.path.join('cextern', 'cfitsio')
    cfitsio_files = glob(os.path.join('cextern', 'cfitsio', '*.c'))
    compression_module_files = [os.path.relpath(fname) for fname in
                                glob(os.path.join(os.path.dirname(__file__),
                                                  'src', '*.c'))]

    return [
        Extension(
            'astropy.io.fits.compression',
            cfitsio_files + compression_module_files,
            # 'numpy' will be replaced with the proper path to the numpy
            # includes
            include_dirs=[cfitsio_path, 'numpy'],
            extra_compile_args=extra_compile_args)
    ]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.io.fits.tests': [os.path.join('data', '*.fits')]}


def get_legacy_alias():
    return setup_helpers.add_legacy_alias(
        'pyfits', 'astropy.io.fits', '3.2.dev', {'__svn_revision__': '1927'})

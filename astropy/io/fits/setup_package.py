# Licensed under a 3-clause BSD style license - see PYFITS.rst

import os

from distutils.core import Extension
from glob import glob

from astropy import setup_helpers


def get_extensions():
    # 'numpy' will be replaced with the proper path to the numpy includes
    include_dirs = ['numpy']
    library_dirs = []
    libraries = []
    source_files = [os.path.relpath(fname) for fname in
                    glob(os.path.join(os.path.dirname(__file__),
                                      'src', '*.c'))]
    extra_compile_args = []

    if not setup_helpers.use_system_library('cfitsio'):
        if setup_helpers.get_compiler_option() == 'msvc':
            # These come from the CFITSIO vcc makefile
            extra_compile_args = [
                    '/D', '"WIN32"',
                    '/D', '"_WINDOWS"',
                    '/D', '"_MBCS"',
                    '/D', '"_USRDLL"',
                    '/D', '"_CRT_SECURE_NO_DEPRECATE"']
        else:
            # All of these switches are to silence warnings from compiling CFITSIO
            extra_compile_args = ['-Wno-unused-variable', '-Wno-parentheses',
                                  '-Wno-uninitialized', '-Wno-format',
                                  '-Wno-strict-prototypes', '-Wno-unused',
                                  '-Wno-comments', '-Wno-switch']

        cfitsio_path = os.path.join('cextern', 'cfitsio')
        cfitsio_files = glob(os.path.join(cfitsio_path, '*.c'))
        include_dirs.append(cfitsio_path)
        source_files.extend(cfitsio_files)
    else:
        setup_helpers.pkg_config(['cfitsio'], ['cfitsio'], include_dirs,
                                 library_dirs, libraries)

    return [
        Extension(
            'astropy.io.fits.compression',
            source_files,
            include_dirs=include_dirs,
            extra_compile_args=extra_compile_args,
            libraries=libraries,
            library_dirs=library_dirs)
    ]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.io.fits.tests': [os.path.join('data', '*.fits')]}


def get_external_libraries():
    return ['cfitsio']


def get_legacy_alias():
    return setup_helpers.add_legacy_alias(
        'pyfits', 'astropy.io.fits', '3.2.dev', {'__svn_revision__': '1927'})

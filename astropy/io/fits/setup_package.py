# Licensed under a 3-clause BSD style license - see PYFITS.rst
from __future__ import absolute_import

import os

from distutils.core import Extension
from glob import glob

from astropy_helpers import setup_helpers
from astropy_helpers.distutils_helpers import get_distutils_build_option


def _get_compression_extension():
    # 'numpy' will be replaced with the proper path to the numpy includes
    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append('numpy')
    cfg['sources'].append(os.path.join(os.path.dirname(__file__), 'src',
                                       'compressionmodule.c'))

    if not setup_helpers.use_system_library('cfitsio'):
        if setup_helpers.get_compiler_option() == 'msvc':
            # These come from the CFITSIO vcc makefile
            cfg['extra_compile_args'].extend([
                    '/D', '"WIN32"',
                    '/D', '"_WINDOWS"',
                    '/D', '"_MBCS"',
                    '/D', '"_USRDLL"',
                    '/D', '"_CRT_SECURE_NO_DEPRECATE"'])
        else:
            cfg['extra_compile_args'].extend([
                '-Wno-declaration-after-statement'
            ])

            if not get_distutils_build_option('debug'):
                # All of these switches are to silence warnings from compiling
                # CFITSIO
                cfg['extra_compile_args'].extend([
                    '-Wno-unused-variable', '-Wno-parentheses',
                    '-Wno-uninitialized', '-Wno-format',
                    '-Wno-strict-prototypes', '-Wno-unused', '-Wno-comments',
                    '-Wno-switch', '-Wno-strict-aliasing', '-Wno-return-type',
                    '-Wno-address'
                ])

        cfitsio_path = os.path.join('cextern', 'cfitsio')
        cfitsio_files = glob(os.path.join(cfitsio_path, '*.c'))
        cfg['include_dirs'].append(cfitsio_path)
        cfg['sources'].extend(cfitsio_files)
    else:
        cfg.update(setup_helpers.pkg_config(['cfitsio'], ['cfitsio']))

    return Extension('astropy.io.fits.compression', **cfg)


def get_extensions():
    return [_get_compression_extension()]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.io.fits.tests': [os.path.join('data', '*.fits')]}


def get_external_libraries():
    return ['cfitsio']


def requires_2to3():
    return False

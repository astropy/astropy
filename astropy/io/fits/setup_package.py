# Licensed under a 3-clause BSD style license - see PYFITS.rst

import os

from distutils.core import Extension
from glob import glob

from astropy import setup_helpers


def get_extensions():
    # 'numpy' will be replaced with the proper path to the numpy includes
    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append('numpy')
    cfg['sources'].extend(
        os.path.relpath(fname) for fname in
        glob(os.path.join(os.path.dirname(__file__), 'src', '*.c')))

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
            # All of these switches are to silence warnings from compiling CFITSIO
            cfg['extra_compile_args'].extend([
                '-Wno-unused-variable', '-Wno-parentheses',
                '-Wno-uninitialized', '-Wno-format', '-Wno-strict-prototypes',
                '-Wno-unused', '-Wno-comments', '-Wno-switch'])

        cfitsio_path = os.path.join('cextern', 'cfitsio')
        cfitsio_files = glob(os.path.join(cfitsio_path, '*.c'))
        cfg['include_dirs'].append(cfitsio_path)
        cfg['sources'].extend(cfitsio_files)
    else:
        cfg.update(setup_helpers.pkg_config(['cfitsio'], ['cfitsio']))

    return [Extension('astropy.io.fits.compression', **cfg)]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.io.fits.tests': [os.path.join('data', '*.fits')]}


def get_external_libraries():
    return ['cfitsio']


def get_legacy_alias():
    return setup_helpers.add_legacy_alias(
        'pyfits', 'astropy.io.fits', '3.2.dev', {'__svn_revision__': '1927'})

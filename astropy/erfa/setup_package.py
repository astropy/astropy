# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from distutils.extension import Extension

from astropy_helpers import setup_helpers

ERFAPKGDIR = os.path.relpath(os.path.dirname(__file__))


def _generate_erfa_wrapper():
    import astropy.erfa.cython_generator as cython_generator
    srcdir = cython_generator.DEFAULT_ERFA_LOC
    templateloc = cython_generator.DEFAULT_TEMPLATE_LOC
    cython_generator.main(srcdir, 'astropy/erfa/erfa.py', templateloc)


def get_extensions():

    _generate_erfa_wrapper()

    sources = [os.path.join(ERFAPKGDIR, "erfa.pyx")]
    include_dirs = ['numpy']
    libraries = []

    if setup_helpers.use_system_library('erfa'):
        libraries.append('erfa')
    else:
        # get all of the .c files in the cextern/erfa directory
        erfafns = os.listdir(os.path.join(ERFAPKGDIR, '..', '..', 'cextern', 'erfa'))
        sources.extend(['cextern/erfa/'+fn for fn in erfafns if fn.endswith('.c')])

        include_dirs.append('cextern/erfa')

    erfa_ext = Extension(
        name="astropy.erfa._erfa",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language="c",)

    return [erfa_ext]


def get_external_libraries():
    return ['erfa']


def requires_2to3():
    return False


def get_package_data():
    return {'astropy.erfa': ['erfa.pyx.templ', 'erfa.pyx.templ']}

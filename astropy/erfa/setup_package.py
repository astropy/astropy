# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import glob

from distutils.extension import Extension

from astropy_helpers import setup_helpers

ERFAPKGDIR = os.path.relpath(os.path.dirname(__file__))

ERFA_SRC = os.path.abspath(os.path.join(ERFAPKGDIR, '..','..','cextern','erfa'))

SRC_FILES = glob.glob(os.path.join(ERFA_SRC, '*'))
SRC_FILES += [os.path.join(ERFAPKGDIR, filename)
              for filename in ['erfa.py.templ', 'erfa.pyx.templ', 'cython_generator.py']]

GEN_FILES = [os.path.join(ERFAPKGDIR, 'erfa.py'), os.path.join(ERFAPKGDIR, 'erfa.pyx')]


def preprocess_source():

    # Generating the ERFA wrappers should only be done if needed. This also
    # ensures that it is not done for any release tarball since those will
    # include erfa.py and erfa.pyx.
    if all(os.path.exists(filename) for filename in GEN_FILES):

        # Determine modification times
        erfa_mtime = max(os.path.getmtime(filename) for filename in SRC_FILES)
        gen_mtime = min(os.path.getmtime(filename) for filename in GEN_FILES)

        # If generated source is recent enough, don't update
        if gen_mtime > erfa_mtime:
            return

    import imp
    gen = imp.load_source('cython_generator',
                          os.path.join(ERFAPKGDIR, 'cython_generator.py'))

    gen.main(gen.DEFAULT_ERFA_LOC,
             os.path.join(ERFAPKGDIR, 'erfa.py'),
             gen.DEFAULT_TEMPLATE_LOC,
             verbose=False)


def get_extensions():
    sources = [os.path.join(ERFAPKGDIR, "erfa.pyx")]
    include_dirs = ['numpy']
    libraries = []

    if setup_helpers.use_system_library('erfa'):
        libraries.append('erfa')
    else:
        # get all of the .c files in the cextern/erfa directory
        erfafns = os.listdir(ERFA_SRC)
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

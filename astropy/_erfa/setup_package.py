# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import glob

from distutils.extension import Extension
from distutils.dep_util import newer

import numpy
from extension_helpers import import_file

ERFAPKGDIR = os.path.relpath(os.path.dirname(__file__))

ERFA_SRC = os.path.abspath(os.path.join(ERFAPKGDIR, '..', '..',
                                        'cextern', 'erfa'))

SRC_FILES = glob.glob(os.path.join(ERFA_SRC, '*'))
SRC_FILES += [os.path.join(ERFAPKGDIR, filename)
              for filename in ['pav2pv.c', 'pv2pav.c', 'erfa_additions.h',
                               'ufunc.c.templ', 'core.py.templ',
                               'erfa_generator.py']]

GEN_FILES = [os.path.join(ERFAPKGDIR, 'core.py'),
             os.path.join(ERFAPKGDIR, 'ufunc.c')]


def get_extensions():
    gen_files_exist = all(os.path.isfile(fn) for fn in GEN_FILES)
    gen_files_outdated = False
    if os.path.isdir(ERFA_SRC):
        # assume thet 'erfaversion.c' is updated at each release at least
        src = os.path.join(ERFA_SRC, 'erfaversion.c')
        gen_files_outdated = any(newer(src, fn) for fn in GEN_FILES)
    elif not gen_files_exist:
        raise RuntimeError(
            'Missing "liberfa" source files, unable to generate '
            '"erfa/ufunc.c" and "erfa/core.py". '
            'Please check your source tree. '
            'Maybe "git submodule update" could help.')

    if not gen_files_exist or gen_files_outdated:
        print('Run "erfa_generator.py"')
        #cmd = [sys.executable, 'erfa_generator.py', ERFA_SRC, '--quiet']
        #subprocess.run(cmd, check=True)

        gen = import_file(os.path.join(ERFAPKGDIR, 'erfa_generator.py'))
        gen.main(verbose=False)

    sources = [os.path.join(ERFAPKGDIR, fn)
               for fn in ("ufunc.c", "pav2pv.c", "pv2pav.c")]
    include_dirs = [numpy.get_include()]
    libraries = []

    if (int(os.environ.get('ASTROPY_USE_SYSTEM_ERFA', 0)) or
            int(os.environ.get('ASTROPY_USE_SYSTEM_ALL', 0))):
        libraries.append('erfa')
    else:
        # get all of the .c files in the cextern/erfa directory
        erfafns = os.listdir(ERFA_SRC)
        sources.extend(['cextern/erfa/' + fn
                        for fn in erfafns if fn.endswith('.c')])

        include_dirs.append('cextern/erfa')

    erfa_ext = Extension(
        name="astropy._erfa.ufunc",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language="c",)

    return [erfa_ext]

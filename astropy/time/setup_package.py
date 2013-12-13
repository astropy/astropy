# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from distutils.extension import Extension

from astropy_helpers import setup_helpers

TIMEROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = [os.path.join(TIMEROOT, "erfa_time.pyx")]
    include_dirs = ['numpy']
    libraries = []

    if setup_helpers.use_system_library('erfa'):
        libraries.append('erfa')
    else:
        sources.append("cextern/erfa/erfa.c")
        include_dirs.append('cextern/erfa')

    time_ext = Extension(
        name="astropy.time.erfa_time",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language="c",)

    return [time_ext]


def get_external_libraries():
    return ['erfa']

def requires_2to3():
    return False

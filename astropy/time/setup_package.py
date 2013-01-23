# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from distutils.extension import Extension

from astropy import setup_helpers

TIMEROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = [os.path.join(TIMEROOT, "sofa_time.pyx")]
    include_dirs = ['numpy']
    libraries = []

    if setup_helpers.use_system_library('sofa'):
        libraries.append('sofa_c')
    else:
        sources.append("cextern/sofa/sofa.c")
        include_dirs.append('cextern/sofa')

    time_ext = Extension(
        name="astropy.time.sofa_time",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language="c",)

    return [time_ext]

def get_external_libraries():
    return ['sofa']

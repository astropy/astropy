# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from distutils.extension import Extension

C_CONVOLVE_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_CONVOLVE_PKGDIR, filename) for filename in ['src/boundary_none.c', 'src/boundary_padded.c']]


def get_extensions():
    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations
    lib_convolve_none = Extension(name='astropy.convolution.lib_convolve_none',
                                  sources=[os.path.join(C_CONVOLVE_PKGDIR, 'src/boundary_none.c')],
                                  include_dirs=["numpy"],
                                  language='c')
    lib_convolve_padded = Extension(name='astropy.convolution.lib_convolve_padded',
                                    sources=[os.path.join(C_CONVOLVE_PKGDIR, 'src/boundary_padded.c')],
                                    include_dirs=["numpy"],
                                    language='c')

    return [lib_convolve_none, lib_convolve_padded]

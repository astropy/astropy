# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
from distutils.extension import Extension
from astropy_helpers.openmp_helpers import add_openmp_flags_if_available

C_CONVOLVE_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_CONVOLVE_PKGDIR, filename)
              for filename in ['src/boundary_none.c',
                               'src/boundary_padded.c']]

extra_compile_args=['-UNDEBUG']
if not sys.platform.startswith('win'):
    extra_compile_args.append('-fPIC')

def get_extensions():
    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations
    lib_convolve_ext = Extension(name='astropy.convolution.lib_convolve', sources=SRC_FILES,
                 extra_compile_args=extra_compile_args,
                 include_dirs=["numpy"],
                 language='c')

    add_openmp_flags_if_available(lib_convolve_ext)

    return [lib_convolve_ext]

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
from distutils.extension import Extension

C_CONVOLVE_PKGDIR = os.path.relpath(os.path.dirname(__file__))

SRC_FILES = [os.path.join(C_CONVOLVE_PKGDIR, filename)
              for filename in ['src/convolve.c']]

extra_compile_args=['-UNDEBUG']
if not sys.platform.startswith('win'):
    extra_compile_args.append('-fPIC')

def get_extensions():
    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations
    _convolve_ext = Extension(name='astropy.convolution._convolve', sources=SRC_FILES,
                              extra_compile_args=extra_compile_args,
                              include_dirs=["numpy"],
                              language='c')

    return [_convolve_ext]

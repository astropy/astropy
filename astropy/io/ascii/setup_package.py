# Licensed under a 3-clause BSD style license

import os
from distutils.extension import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = [os.path.join(ROOT, 'cparser.pyx'),
               os.path.join(ROOT, 'src', 'tokenizer.c')]
    ascii_ext = Extension(
        name="astropy.io.ascii.cparser",
        include_dirs=["numpy"],
        sources=sources)
    return [ascii_ext]

# Licensed under a 3-clause BSD style license

import os

from setuptools import Extension

SRC_DIR = os.path.relpath(os.path.join(os.path.dirname(__file__), "src"))


def get_extensions():
    return [
        Extension(
            "astropy.io.fits.hdu.compressed._compression",
            sources=[
                os.path.join(SRC_DIR, "compression.c"),
                os.path.join(SRC_DIR, "unquantize.c"),
                os.path.join("cextern", "cfitsio", "lib", "pliocomp.c"),
                os.path.join("cextern", "cfitsio", "lib", "ricecomp.c"),
                os.path.join("cextern", "cfitsio", "lib", "fits_hcompress.c"),
                os.path.join("cextern", "cfitsio", "lib", "fits_hdecompress.c"),
                os.path.join("cextern", "cfitsio", "lib", "quantize.c"),
            ],
            include_dirs=[SRC_DIR],
        )
    ]

# Licensed under a 3-clause BSD style license

import os
import sysconfig

from setuptools import Extension

SRC_DIR = os.path.relpath(os.path.join(os.path.dirname(__file__), "src"))

USE_PY_LIMITED_API = not sysconfig.get_config_var("Py_GIL_DISABLED")


def get_extensions():
    kwargs = {}
    if USE_PY_LIMITED_API:
        kwargs["py_limited_api"] = True
        kwargs['define_macros'] = [("Py_LIMITED_API", "0x030B0000")]

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
            **kwargs
        )
    ]

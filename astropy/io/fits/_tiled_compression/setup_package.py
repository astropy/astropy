# Licensed under a 3-clause BSD style license

import os
from collections import defaultdict

from setuptools import Extension

from extension_helpers import pkg_config

SRC_DIR = os.path.join(os.path.dirname(__file__), "src")


def get_extensions():

    cfg = defaultdict(list)
    cfg["sources"].extend(
        [
            os.path.join(SRC_DIR, "compression.c"),
            os.path.join(SRC_DIR, "unquantize.c"),
        ]
    )

    if int(os.environ.get("ASTROPY_USE_SYSTEM_CFITSIO", 0)) or int(
        os.environ.get("ASTROPY_USE_SYSTEM_ALL", 0)
    ):
        for k, v in pkg_config(["cfitsio"], ["cfitsio"]).items():
            cfg[k].extend(v)
    else:
        cfg["sources"].extend(
            [
                os.path.join("cextern", "cfitsio", "lib", "pliocomp.c"),
                os.path.join("cextern", "cfitsio", "lib", "ricecomp.c"),
                os.path.join("cextern", "cfitsio", "lib", "fits_hcompress.c"),
                os.path.join("cextern", "cfitsio", "lib", "fits_hdecompress.c"),
                os.path.join("cextern", "cfitsio", "lib", "quantize.c"),
            ]
        )

    cfg["include_dirs"].append(SRC_DIR)

    return [Extension("astropy.io.fits._tiled_compression._compression", **cfg)]

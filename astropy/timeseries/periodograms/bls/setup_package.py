# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sysconfig
from pathlib import Path

from numpy import get_include as get_numpy_include
from setuptools import Extension

BLS_ROOT = Path(__file__).parent.resolve().relative_to(Path.cwd())

USE_PY_LIMITED_API = not sysconfig.get_config_var("Py_GIL_DISABLED")


def get_extensions():
    kwargs = {}
    if USE_PY_LIMITED_API:
        kwargs["py_limited_api"] = True
        kwargs["define_macros"] = [("Py_LIMITED_API", "0x030B0000")]

    ext = Extension(
        "astropy.timeseries.periodograms.bls._impl",
        sources=[
            str(BLS_ROOT / "bls.c"),
            str(BLS_ROOT / "_impl.pyx"),
        ],
        include_dirs=[get_numpy_include()],
        **kwargs,
    )
    return [ext]

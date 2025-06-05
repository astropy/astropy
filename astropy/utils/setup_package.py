# Licensed under a 3-clause BSD style license - see LICENSE.rst

from os.path import dirname, join, relpath
import sysconfig

from setuptools import Extension

ASTROPY_UTILS_ROOT = dirname(__file__)
USE_PY_LIMITED_API = not sysconfig.get_config_var("Py_GIL_DISABLED")


def get_extensions():
    kwargs = {}
    if USE_PY_LIMITED_API:
        kwargs["py_limited_api"] = True
        kwargs['define_macros'] = [("Py_LIMITED_API", "0x030B0000")]

    return [
        Extension(
            "astropy.utils._compiler",
            [relpath(join(ASTROPY_UTILS_ROOT, "src", "compiler.c"))],
            **kwargs
        )
    ]

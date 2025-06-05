# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Copied from astropy/convolution/setup_package.py

from pathlib import Path

from numpy import get_include as get_numpy_include
from setuptools import Extension
import sysconfig

C_TIME_PKGDIR = Path(__file__).parent.resolve().relative_to(Path.cwd())

SRC_FILES = [str(C_TIME_PKGDIR / "src" / "parse_times.c")]

USE_PY_LIMITED_API = not sysconfig.get_config_var("Py_GIL_DISABLED")


def get_extensions():
    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations

    kwargs = {}
    if USE_PY_LIMITED_API:
        kwargs["py_limited_api"] = True
        kwargs['define_macros'] = [("Py_LIMITED_API", "0x030B0000")]

    _time_ext = Extension(
        name="astropy.time._parse_times",
        sources=SRC_FILES,
        include_dirs=[get_numpy_include()],
        language="c",
        **kwargs
    )
    return [_time_ext]

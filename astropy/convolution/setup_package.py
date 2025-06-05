# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
from pathlib import Path
import sysconfig

from numpy import get_include as get_numpy_include
from setuptools import Extension

C_CONVOLVE_PKGDIR = Path(__file__).parent.resolve().relative_to(Path.cwd())

extra_compile_args = ["-UNDEBUG"]
if not sys.platform.startswith("win"):
    extra_compile_args.append("-fPIC")

USE_PY_LIMITED_API = not sysconfig.get_config_var("Py_GIL_DISABLED")


def get_extensions():
    kwargs = {}
    if USE_PY_LIMITED_API:
        kwargs["py_limited_api"] = True
        kwargs['define_macros'] = [("Py_LIMITED_API", "0x030B0000")]

    # Add '-Rpass-missed=.*' to ``extra_compile_args`` when compiling with clang
    # to report missed optimizations
    sources = [
        str(C_CONVOLVE_PKGDIR / "_convolve.pyx"),
        str(C_CONVOLVE_PKGDIR / "src" / "convolve.c"),
    ]
    _convolve_ext = Extension(
        name="astropy.convolution._convolve",
        extra_compile_args=extra_compile_args,
        include_dirs=[get_numpy_include()],
        sources=sources,
        **kwargs
    )
    return [_convolve_ext]

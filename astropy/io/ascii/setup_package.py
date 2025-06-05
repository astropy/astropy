# Licensed under a 3-clause BSD style license

from pathlib import Path
import sysconfig

from numpy import get_include as get_numpy_include
from setuptools import Extension

ROOT = Path(__file__).parent.relative_to(Path.cwd())
USE_PY_LIMITED_API = not sysconfig.get_config_var("Py_GIL_DISABLED")


def get_extensions():
    kwargs = {}
    if USE_PY_LIMITED_API:
        kwargs["py_limited_api"] = True
        kwargs['define_macros'] = [("Py_LIMITED_API", "0x030B0000")]

    sources = [
        str(ROOT / "cparser.pyx"),
        str(ROOT / "src" / "tokenizer.c"),
    ]
    ascii_ext = Extension(
        name="astropy.io.ascii.cparser",
        include_dirs=[get_numpy_include()],
        sources=sources,
        **kwargs
    )
    return [ascii_ext]

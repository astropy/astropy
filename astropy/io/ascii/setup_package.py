# Licensed under a 3-clause BSD style license

from pathlib import Path

from numpy import get_include as get_numpy_include
from setuptools import Extension

ROOT = Path(__file__).parent.relative_to(Path.cwd())


def get_extensions():
    sources = [
        str(ROOT / "cparser.pyx"),
        str(ROOT / "src" / "tokenizer.c"),
    ]
    ascii_ext = Extension(
        name="astropy.io.ascii.cparser",
        include_dirs=[get_numpy_include()],
        sources=sources,
    )
    return [ascii_ext]

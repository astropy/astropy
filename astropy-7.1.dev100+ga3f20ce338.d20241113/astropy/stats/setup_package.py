# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

from pathlib import Path

from numpy import get_include as get_numpy_include
from setuptools import Extension

ROOT = Path(__file__).parent.resolve().relative_to(Path.cwd())
SRCFILES = ["wirth_select.c", "compute_bounds.c", "fast_sigma_clip.c"]
SRCFILES = [str(ROOT / "src" / srcfile) for srcfile in SRCFILES]


def get_extensions() -> list[Extension, Extension]:
    _sigma_clip_ext = Extension(
        name="astropy.stats._fast_sigma_clip",
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        sources=SRCFILES,
        include_dirs=[get_numpy_include()],
        language="c",
    )
    _stats_ext = Extension(
        name="astropy.stats._stats",
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        sources=[str(ROOT / "_stats.pyx")],
        include_dirs=[get_numpy_include()],
    )

    return [_sigma_clip_ext, _stats_ext]

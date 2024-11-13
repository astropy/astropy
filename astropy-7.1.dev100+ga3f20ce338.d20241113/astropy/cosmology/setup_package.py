# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
from os.path import relpath
from pathlib import Path

from setuptools import Extension

ASTROPY_COSMOLOGY_ROOT = Path(__file__).parent


if sys.platform.startswith("win"):
    # on windows, -Werror (and possibly -Wall too) isn't recognized
    extra_compile_args = []
else:
    # be extra careful with this extension as it calls PyErr_WarnEx
    # with a formatted message whose size cannot be determined at compile time,
    # which is never done within the standard library
    extra_compile_args = ["-Werror", "-Wall"]


def get_extensions():
    return [
        Extension(
            "astropy.cosmology._signature_deprecations",
            [relpath(Path(ASTROPY_COSMOLOGY_ROOT, "_signature_deprecations.c"))],
            extra_compile_args=extra_compile_args,
        ),
    ]

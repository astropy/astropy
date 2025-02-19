#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the pyproject.toml file.

from setuptools import setup

from extension_helpers import get_extensions

ext_modules = get_extensions()

# Specify the minimum version for the Numpy C-API
for ext in ext_modules:
    if ext.include_dirs and "numpy" in ext.include_dirs[0]:
        ext.define_macros.append(("NPY_TARGET_VERSION", "NPY_1_23_API_VERSION"))
        ext.define_macros.append(("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"))

setup(ext_modules=ext_modules)

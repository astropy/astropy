#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the pyproject.toml file.

import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.install import install

from extension_helpers import get_extensions

ext_modules = get_extensions()


def install_stubs(build_lib, output_dir):
    sys.path.insert(0, build_lib)
    try:
        from astropy.units.scripts.typestubs_for_units import (  # noqa: PLC0415
            create_stubs,
        )

        create_stubs(str(Path(output_dir) / "astropy" / "units"), log_level="WARN")
    finally:
        # Undo the path modification.
        sys.path.pop(0)


class InstallWithStubs(install):
    """Post-installation command for installation mode."""

    def run(self):
        super().run()
        install_stubs(self.build_lib, self.root)


# Specify the minimum version for the Numpy C-API
for ext in ext_modules:
    if ext.include_dirs and "numpy" in ext.include_dirs[0]:
        ext.define_macros.append(("NPY_TARGET_VERSION", "NPY_1_24_API_VERSION"))
        ext.define_macros.append(("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"))

setup(
    ext_modules=ext_modules,
    cmdclass={
        "install": InstallWithStubs,
    },
)

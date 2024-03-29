#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the pyproject.toml file.

import sys

# First provide helpful messages if contributors try and run legacy commands
# for tests or docs.

TEST_HELP = """
Note: running tests is no longer done using 'python setup.py test'. Instead
you will need to run:

    tox -e test

If you don't already have tox installed, you can install it with:

    pip install tox

If you only want to run part of the test suite, you can also use pytest
directly with::

    pip install -e .[test]
    pytest

For more information, see:

  https://docs.astropy.org/en/latest/development/testguide.html#running-tests
"""

if "test" in sys.argv:
    print(TEST_HELP)
    sys.exit(1)

DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py build_docs'. Instead you will need to run:

    tox -e build_docs

If you don't already have tox installed, you can install it with:

    pip install tox

You can also build the documentation with Sphinx directly using::

    pip install -e .[docs]
    cd docs
    make html

For more information, see:

  https://docs.astropy.org/en/latest/install.html#builddocs
"""

if "build_docs" in sys.argv or "build_sphinx" in sys.argv:
    print(DOCS_HELP)
    sys.exit(1)


# Only import these if the above checks are okay
# to avoid masking the real problem with import error.
from setuptools import setup  # noqa: E402

from extension_helpers import get_extensions  # noqa: E402

ext_modules = get_extensions()

# Specify the minimum version for the Numpy C-API
for ext in ext_modules:
    if ext.include_dirs and "numpy" in ext.include_dirs[0]:
        ext.define_macros.append(("NPY_TARGET_VERSION", "NPY_1_23_API_VERSION"))

setup(ext_modules=ext_modules)

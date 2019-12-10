#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys

from setuptools import setup
from extension_helpers import get_extensions

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

    pip install -e .
    pytest

For more information, see:

  http://docs.astropy.org/en/latest/development/testguide.html#running-tests
"""

if 'test' in sys.argv:
    print(TEST_HELP)
    sys.exit(1)

DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py build_docs'. Instead you will need to run:

    tox -e build_docs

If you don't already have tox installed, you can install it with:

    pip install tox

For more information, see:

  http://docs.astropy.org/en/latest/install.html#builddocs
"""

if 'build_docs' in sys.argv or 'build_sphinx' in sys.argv:
    print(DOCS_HELP)
    sys.exit(1)

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the setup.cfg file. Here we mainly set up
# setup_requires and install_requires since these are determined
# programmatically.

setup(use_scm_version={'write_to': os.path.join('astropy', 'version.py')}, ext_modules=get_extensions())

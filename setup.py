#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the setup.cfg file. Here we mainly set up
# setup_requires and install_requires since these are determined
# programmatically.

import os
import builtins

# Because we have a pyproject.toml file, the isolated build environment
# doesn't allow the ah_bootstrap file to be imported unless the current
# directory is added to the Python path
import sys
sys.path.append(os.path.abspath("."))
import ah_bootstrap  # noqa

from astropy_helpers.distutils_helpers import is_distutils_display_option
from astropy_helpers.setup_helpers import setup

from setuptools.config import read_configuration

# We set up the following variable because we then use this in astropy/__init__.py
# to make sure that we aren't importing astropy during the setup process (we used
# to do this)
builtins._ASTROPY_CORE_SETUP_ = True

if is_distutils_display_option():
    # Avoid installing setup_requires dependencies if the user just
    # queries for information
    setup_requires = []
else:
    setup_requires = read_configuration('setup.cfg')['options']['setup_requires']
    # Make sure we have the packages needed for building astropy, but do not
    # require them when installing from an sdist as the c files are included.
    if not os.path.exists(os.path.join(os.path.dirname(__file__), 'PKG-INFO')):
        setup_requires.extend(['cython>=0.29.13', 'jinja2>=2.7'])

setup(setup_requires=setup_requires)

#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import builtins
from setuptools import setup

# We set up the following variable because we then use this in astropy/__init__.py
# to make sure that we aren't importing astropy during the setup process (we used
# to do this)
builtins._ASTROPY_CORE_SETUP_ = True

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the setup.cfg file. Here we mainly set up
# setup_requires and install_requires since these are determined
# programmatically.

from extension_helpers import get_extensions
setup(use_scm_version={'write_to': os.path.join('astropy', 'version.py')}, ext_modules=get_extensions())

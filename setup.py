#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys

from distutils.version import LooseVersion

# We require setuptools 30.3.0 or later for the configuration in setup.cfg to
# work properly.
import setuptools
if LooseVersion(setuptools.__version__) < LooseVersion('30.3.0'):
    sys.stderr.write("ERROR: Astropy requires setuptools 30.3.0 or later "
                     "(found {0})".format(setuptools.__version__))
    sys.exit(1)

from setuptools.config import read_configuration
conf = read_configuration('setup.cfg')

# This is the same check as astropy/__init__.py but this one has to
# happen before importing ah_bootstrap
minimum_python_version = conf['options']['python_requires'].replace('>', '').replace('<', '').replace('=', '')
if sys.version_info < tuple((int(val) for val in minimum_python_version.split('.'))):
    sys.stderr.write("ERROR: Astropy requires Python {} or later\n".format(
        minimum_python_version))
    sys.exit(1)

import os
import glob

import ah_bootstrap
from setuptools import setup

from astropy_helpers.setup_helpers import (
    register_commands, get_package_info, get_debug_option)
from astropy_helpers.distutils_helpers import is_distutils_display_option
from astropy_helpers.git_helpers import get_git_devstr
from astropy_helpers.version_helpers import generate_version_py

import astropy

NAME = conf['metadata']['name']

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = conf['metadata']['version']

# Indicates if this version is a release version
RELEASE = 'dev' not in VERSION

if not RELEASE:
    VERSION += get_git_devstr(False)

# Populate the dict of setup command overrides; this should be done before
# invoking any other functionality from distutils since it can potentially
# modify distutils' behavior.
cmdclassd = register_commands(NAME, VERSION, RELEASE)

# Freeze build information in version.py
generate_version_py(NAME, VERSION, RELEASE, get_debug_option(NAME),
                    uses_git=not RELEASE)

# Get configuration information from all of the various subpackages.
# See the docstring for setup_helpers.update_package_files for more
# details.
package_info = get_package_info()

# Add the project-global data
package_info['package_data'].setdefault('astropy', []).append('data/*')

min_numpy_version = 'numpy>=' + astropy.__minimum_numpy_version__

if is_distutils_display_option():
    # Avoid installing setup_requires dependencies if the user just
    # queries for information
    setup_requires = []
else:
    setup_requires = [min_numpy_version]

    # Make sure we have the packages needed for building astropy, but do not require them
    # when installing from an sdist as the c files are included there.
    if not os.path.exists(os.path.join(os.path.dirname(__file__), 'PKG-INFO')):
        setup_requires.extend(['cython>=0.21', 'jinja2>=2.7'])

install_requires = [min_numpy_version]

# Avoid installing setup_requires dependencies if the user just
# queries for information
if is_distutils_display_option():
    setup_requires = []

setup(version=VERSION,
      setup_requires=setup_requires,
      install_requires=install_requires,
      long_description=astropy.__doc__,
      cmdclass=cmdclassd,
      **package_info
)

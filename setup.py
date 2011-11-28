#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Use "distribute" - the setuptools fork that supports python 3.
from distribute_setup import use_setuptools
use_setuptools()

import os
import glob
from setuptools import setup, find_packages

import astropy
from astropy import setup_helpers
from astropy.tests.helper import astropy_test
from astropy.version_helper import get_git_devstr, generate_version_py

version = '0.0dev'

# Indicates if this version is a release version
release = 'dev' not in version

# Adjust the compiler in case the default on this platform is to use a
# broken one.
setup_helpers.adjust_compiler()

# Indicate that we are in building mode
setup_helpers.set_build_mode()

if not release:
    version += get_git_devstr(False)
generate_version_py('astropy', version, release,
                    setup_helpers.get_debug_option())

# Use the find_packages tool to locate all packages and modules
packagenames = find_packages()

# Treat everything in scripts except README.rst as a script to be installed
scripts = glob.glob('scripts/*')
scripts.remove('scripts/README.rst')

# Check that Numpy is installed.
# NOTE: We cannot use setuptools/distribute/packaging to handle this
# dependency for us, since some of the subpackages need to be able to
# access numpy at build time, and they are configured before
# setuptools has a chance to check and resolve the dependency.
setup_helpers.check_numpy()

# This dictionary stores the command classes used in setup below
cmdclassd = {'test': astropy_test}

# Additional C extensions that are not Cython-based should be added here.
extensions = []

# A dictionary to keep track of all package data to install
package_data = {'astropy': ['data/*']}

# A dictionary to keep track of extra packagedir mappings
package_dirs = {}

# Extra files to install - distutils calls them "data_files", but this
# shouldn't be used for data files - rather any other files that should be
# installed in a special place
data_files = []

# Update extensions, package_data, and data_files from any sub-packages that
# define their own extension modules and package data.  See the docstring for
# setup_helpers.update_package_files for more details.
setup_helpers.update_package_files('astropy', extensions, package_data,
                                   data_files, packagenames, package_dirs)

if setup_helpers.HAVE_CYTHON and not release:
    from Cython.Distutils import build_ext
    # Builds Cython->C if in dev mode and Cython is present
    cmdclassd['build_ext'] = build_ext

if setup_helpers.AstropyBuildSphinx is not None:
    cmdclassd['build_sphinx'] = setup_helpers.AstropyBuildSphinx


setup(name='astropy',
      version=version,
      description='Community-developed python astronomy tools',
      packages=packagenames,
      package_data=package_data,
      package_dir=package_dirs,
      data_files=data_files,
      ext_modules=extensions,
      scripts=scripts,
      requires=['numpy'],  # scipy not required, but strongly recommended
      install_requires=['numpy'],
      provides=['astropy'],
      author='The Astropy Team',
      author_email='astropy.team@gmail.com',
      license='BSD',
      url='http://astropy.org',
      long_description=astropy.__doc__,
      cmdclass=cmdclassd,
      zip_safe=False,
      use_2to3=True
      )

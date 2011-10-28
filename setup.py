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
from astropy.version_helper import _get_git_devstr, _generate_version_py


VERSION = '0.0dev'
RELEASE = not VERSION.endswith('dev')

if not RELEASE:
    VERSION += _get_git_devstr(False)
_generate_version_py(VERSION, RELEASE, setup_helpers.get_debug_option())

# Use the find_packages tool to locate all packages and modules other than
# those that are in tests/
packages = find_packages(exclude=['tests'])

# Treat everything in scripts except README.rst as a script to be installed
scripts = glob.glob('scripts/*')
scripts.remove('scripts/README.rst')

# Check that Numpy is installed.
# NOTE: We can not use setuptools/distribute/packaging to handle this
# dependency for us, since some of the subpackages need to be able to
# access numpy at build time, and they are configured before
# setuptools has a chance to check and resolve the dependency.
setup_helpers.check_numpy()

# This dictionary stores the command classes used in setup below
cmdclassd = {'test': astropy_test}

# A dictionary to keep track of all package data to install
package_data = {'astropy': ['data/*']}

# Additional C extensions that are not Cython-based should be added here.
extensions = []

# Extra data files
data_files = []

# For each of the setup_package.py modules, extract any information
# that is needed to install them.
for package in setup_helpers.iter_setup_packages():
    if hasattr(package, 'get_extensions'):
        extensions.extend(package.get_extensions())
    if hasattr(package, 'get_package_data'):
        package_data.update(package.get_package_data())
    if hasattr(package, 'get_data_files'):
        data_files.extend(package.get_data_files())

extensions.extend(setup_helpers.get_cython_extensions())

if setup_helpers.HAVE_CYTHON and not RELEASE:
    from Cython.Distutils import build_ext
    cmdclassd['build_ext'] = build_ext

# Implement a version of build_sphinx that automatically creates the
# docs/_build dir - this is needed because github won't create the _build dir
# because it has no tracked files

try:

    from sphinx.setup_command import BuildDoc

    class astropy_build_sphinx(BuildDoc):

        def finalize_options(self):

            from distutils.cmd import DistutilsOptionError

            if self.build_dir is not None:
                if os.path.isfile(self.build_dir):
                    raise DistutilsOptionError('Attempted to build_sphinx ' + \
                                               'into a file ' + self.build_dir)
                self.mkpath(self.build_dir)

            return BuildDoc.finalize_options(self)

    cmdclassd['build_sphinx'] = astropy_build_sphinx

except ImportError:  # Sphinx not present
    pass


setup(name='astropy',
      version=VERSION,
      description='Community-developed python astronomy tools',
      packages=packages,
      package_data=package_data,
      ext_modules=extensions,
      scripts=scripts,
      requires=['numpy', 'scipy'],
      install_requires=['numpy'],
      provides=['astropy'],
      author='The Astropy Team',
      author_email='astropy.team@gmail.com',
      #license = '', #TODO: decide on a license
      url='http://astropy.org',
      long_description=astropy.__doc__,
      cmdclass=cmdclassd,
      zip_safe=False,
      data_files=data_files
      )

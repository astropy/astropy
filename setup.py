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

version = '0.0dev'
#indicates if this version is a release version
release = 'dev' not in version

# Adjust the compiler in case the default on this platform is to use a
# broken one.
setup_helpers.adjust_compiler()
print 'here'
if not release:
    version += _get_git_devstr(False)
_generate_version_py(version, release, setup_helpers.get_debug_option())

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

# A dictionary to keep track of all package data to install
package_data = {'astropy': ['data/*']}

# Additional C extensions that are not Cython-based should be added here.
extensions = []

# Extra data files
data_files = []

# For each of the setup_package.py modules, extract any information
# that is needed to install them.
for setuppkg in setup_helpers.iter_setup_packages():
    if hasattr(setuppkg, 'get_extensions'):
        extensions.extend(setuppkg.get_extensions())
    
    if hasattr(setuppkg, 'get_package_data'):
        #TBD: decide if this should be removed in favor of the data loading mechanism in config.data
        package_data.update(setuppkg.get_package_data())
    if hasattr(setuppkg, 'get_data_files'):
        data_files.extend(setuppkg.get_data_files())

extensions.extend(setup_helpers.get_cython_extensions())

if setup_helpers.HAVE_CYTHON and not release:
    from Cython.Distutils import build_ext
    cmdclassd['build_ext'] = build_ext

# Implement a version of build_sphinx that automatically creates the
# docs/_build dir - this is needed because github won't create the _build dir
# because it has no tracked files

try:

    from sphinx.setup_command import BuildDoc

    class AstropyBuildSphinx(BuildDoc):
        """
        This class 
        """
        def finalize_options(self):
            from distutils.cmd import DistutilsOptionError

            if self.build_dir is not None:
                if os.path.isfile(self.build_dir):
                    raise DistutilsOptionError('Attempted to build_sphinx ' + \
                                               'into a file ' + self.build_dir)
                self.mkpath(self.build_dir)

            return BuildDoc.finalize_options(self)

    cmdclassd['build_sphinx'] = AstropyBuildSphinx

except ImportError:  # Sphinx not present
    pass


setup(name='astropy',
      version=version,
      description='Community-developed python astronomy tools',
      packages=packagenames,
      package_data=package_data,
      data_files=data_files,
      ext_modules=extensions,
      scripts=scripts,
      requires=['numpy'], #scipy not required, but strongly recommended
      install_requires=['numpy'],
      provides=['astropy'],
      author='The Astropy Team',
      author_email='astropy.team@gmail.com',
      license = 'BSD',
      url='http://astropy.org',
      long_description=astropy.__doc__,
      cmdclass=cmdclassd,
      zip_safe=False
      )

#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Use "distribute" - the setuptools fork that supports python 3.
from distribute_setup import use_setuptools
use_setuptools()

import os
import glob
from setuptools import setup, find_packages
from distutils import log

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
from numpy import get_include as get_numpy_include
numpy_includes = get_numpy_include()


# This dictionary stores the command classes used in setup below
cmdclassd = {'test': astropy_test}


# Additional C extensions that are not Cython-based should be added here.
extensions = []

# A dictionary to keep track of all package data to install
package_data = {'astropy': ['data/*']}

# Extra files to install - distutils calls them "data_files", but this shouldn't
# be used for data files - rather any other files that should be installed in a
# special place
data_files = []

# For each of the setup_package.py modules, extract any information
# that is needed to install them.
for pkgnm,setuppkg in setup_helpers.iter_setup_packages():
    #get_extensions must include any Cython extensions by their .pyx filename.
    if hasattr(setuppkg, 'get_extensions'):
        extensions.extend(setuppkg.get_extensions())
    
    if hasattr(setuppkg, 'get_package_data'):
        #TBD: decide if this should be removed in favor of the data loading mechanism in config.data
        package_data.update(setuppkg.get_package_data())
    if hasattr(setuppkg, 'get_data_files'):
        data_files.extend(setuppkg.get_data_files())

#locate any .pyx files not already specified, and add their extensions in. 
#The default include dirs include numpy to facilitate numerical work.
extensions.extend(setup_helpers.get_cython_extensions('astropy',extensions,
                                                      [numpy_includes]))

#now remove extensions that have the special name 'skip_cython', as they exist
#only to indicate that the cython extensions shouldn't be built
for i, ext in reversed(list(enumerate(extensions))):
        if ext.name=='skip_cython':
            del extensions[i]

if setup_helpers.HAVE_CYTHON and not release:
    from Cython.Distutils import build_ext
    #builds Cython->C if in dev mode and Cython is present
    cmdclassd['build_ext'] = build_ext
else:
    
    #otherwise, replace .pyx with C-equivalents, unless c files are missing
    todel = []
    for i,ext in enumerate(extensions):
        for j,s in enumerate(ext.sources):
            if i not in todel and s.endswith('.pyx'):
                cfn = s[:-4]+'.c'
                if os.path.isfile(cfn):
                    ext.sources[j] = cfn
                else:
                    msg = 'Could not find c-file {0} for {1}, skipping extension {2}'
                    log.warn(msg.format(cfn,s,ext.name))
                    todel.append(i)
                    
    
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

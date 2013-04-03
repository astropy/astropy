#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import imp
try:
    # This incantation forces distribute to be used (over setuptools) if it is
    # available on the path; otherwise distribute will be downloaded.
    import pkg_resources
    distribute = pkg_resources.get_distribution('distribute')
    if pkg_resources.get_distribution('setuptools') != distribute:
        sys.path.insert(1, distribute.location)
        distribute.activate()
        imp.reload(pkg_resources)
except:  # There are several types of exceptions that can occur here
    from distribute_setup import use_setuptools
    use_setuptools()

import glob
import os
from setuptools import setup, find_packages

#A dirty hack to get around some early import/configurations ambiguities
if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
builtins._ASTROPY_SETUP_ = True

import astropy
from astropy.setup_helpers import (register_commands, adjust_compiler,
                                   filter_packages, update_package_files,
                                   get_debug_option)
from astropy.version_helpers import get_git_devstr, generate_version_py

NAME = 'astropy'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.2.1'

# Indicates if this version is a release version
RELEASE = 'dev' not in VERSION

if not RELEASE:
    VERSION += get_git_devstr(False)

DOWNLOAD_BASE_URL = 'http://pypi.python.org/packages/source/a/astropy'

# Populate the dict of setup command overrides; this should be done before
# invoking any other functionality from distutils since it can potentially
# modify distutils' behavior.
cmdclassd = register_commands(NAME, VERSION, RELEASE)

# Adjust the compiler in case the default on this platform is to use a
# broken one.
adjust_compiler(NAME)

# Freeze build information in version.py
generate_version_py(NAME, VERSION, RELEASE, get_debug_option())

# Use the find_packages tool to locate all packages and modules
packagenames = filter_packages(find_packages())

# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
           if os.path.basename(fname) != 'README.rst']

# Additional C extensions that are not Cython-based should be added here.
extensions = []

# A dictionary to keep track of all package data to install
package_data = {'astropy': ['data/*']}

# A dictionary to keep track of extra packagedir mappings
package_dirs = {}

# Update extensions, package_data, packagenames and package_dirs from
# any sub-packages that define their own extension modules and package
# data.  See the docstring for setup_helpers.update_package_files for
# more details.
update_package_files(NAME, extensions, package_data, packagenames,
                     package_dirs)

# Currently the only entry points installed by Astropy are hooks to
# zest.releaser for doing Astropy's releases
entry_points = {}
for hook in [('prereleaser', 'middle'), ('releaser', 'middle'),
             ('postreleaser', 'before'), ('postreleaser', 'middle')]:
    hook_ep = 'zest.releaser.' + '.'.join(hook)
    hook_name = 'astropy.release.' + '.'.join(hook)
    hook_func = 'astropy.utils.release:' + '_'.join(hook)
    entry_points[hook_ep] = ['%s = %s' % (hook_name, hook_func)]


setup(name=NAME,
      version=VERSION,
      description='Community-developed python astronomy tools',
      packages=packagenames,
      package_data=package_data,
      package_dir=package_dirs,
      ext_modules=extensions,
      scripts=scripts,
      requires=['numpy'],  # scipy not required, but strongly recommended
      install_requires=['numpy'],
      provides=[NAME],
      author='The Astropy Developers',
      author_email='astropy.team@gmail.com',
      license='BSD',
      url='http://astropy.org',
      long_description=astropy.__doc__,
      download_url='%s/astropy-%s.tar.gz' % (DOWNLOAD_BASE_URL, VERSION),
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: C',
          'Programming Language :: Cython',
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: Implementation :: CPython',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics'
      ],
      cmdclass=cmdclassd,
      zip_safe=False,
      use_2to3=True,
      entry_points=entry_points
)

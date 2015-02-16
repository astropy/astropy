#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import glob
import os
import sys

import ah_bootstrap
from setuptools import setup

#A dirty hack to get around some early import/configurations ambiguities
if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
builtins._ASTROPY_SETUP_ = True

import astropy
from astropy_helpers.setup_helpers import (
    register_commands, adjust_compiler, get_package_info, get_debug_option,
    is_distutils_display_option)
from astropy_helpers.git_helpers import get_git_devstr
from astropy_helpers.version_helpers import generate_version_py

NAME = 'astropy'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.4.6.dev'

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
generate_version_py(NAME, VERSION, RELEASE, get_debug_option(NAME),
                    uses_git=not RELEASE)

# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
           if os.path.basename(fname) != 'README.rst']

# Get configuration information from all of the various subpackages.
# See the docstring for setup_helpers.update_package_files for more
# details.
package_info = get_package_info()

# Add the project-global data
package_info['package_data'].setdefault('astropy', []).append('data/*')

# Currently the only entry points installed by Astropy are hooks to
# zest.releaser for doing Astropy's releases
entry_points = {}
for hook in [('prereleaser', 'middle'), ('releaser', 'middle'),
             ('postreleaser', 'before'), ('postreleaser', 'middle')]:
    hook_ep = 'zest.releaser.' + '.'.join(hook)
    hook_name = 'astropy.release.' + '.'.join(hook)
    hook_func = 'astropy.utils.release:' + '_'.join(hook)
    entry_points[hook_ep] = ['%s = %s' % (hook_name, hook_func)]


setup_requires = ['numpy>=' + astropy.__minimum_numpy_version__]
install_requires = ['numpy>=' + astropy.__minimum_numpy_version__]
# Avoid installing setup_requires dependencies if the user just
# queries for information
if is_distutils_display_option():
    setup_requires = []


setup(name=NAME,
      version=VERSION,
      description='Community-developed python astronomy tools',
      scripts=scripts,
      requires=['numpy'],  # scipy not required, but strongly recommended
      setup_requires=setup_requires,
      install_requires=install_requires,
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
      entry_points=entry_points,
      **package_info
)

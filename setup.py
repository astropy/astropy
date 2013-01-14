#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys

# Use "distribute" - the setuptools fork that supports python 3. We allow
# users to specify to use the system installation if needed, mainly for
# package managers.
if '--use-system-distribute' in sys.argv:
    import setuptools
    sys.argv.remove('--use-system-distribute')
else:
    from distribute_setup import use_setuptools
    use_setuptools()
    import setuptools

from distutils.command import sdist

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
from astropy import setup_helpers
from astropy.version_helper import get_git_devstr, generate_version_py

#version should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
version = '0.2.dev'

# Indicates if this version is a release version
release = 'dev' not in version

download_base_url = 'http://pypi.python.org/packages/source/a/astropy'

# Adjust the compiler in case the default on this platform is to use a
# broken one.
setup_helpers.adjust_compiler()

if not release:
    version += get_git_devstr(False)
generate_version_py('astropy', version, release,
                    setup_helpers.get_debug_option())

# Use the find_packages tool to locate all packages and modules
packagenames = find_packages()
packagenames = setup_helpers.filter_packages(packagenames)

# Treat everything in scripts except README.rst as a script to be installed
scripts = glob.glob(os.path.join('scripts', '*'))
scripts.remove(os.path.join('scripts', 'README.rst'))

# This dictionary stores the command classes used in setup below
cmdclassd = {'test': setup_helpers.setup_test_command('astropy'),

             # Use distutils' sdist because it respects package_data.
             # setuptools/distributes sdist requires duplication of
             # information in MANIFEST.in
             'sdist': sdist.sdist,

             # Use a custom build command which understands additional
             # commandline arguments
             'build': setup_helpers.AstropyBuild,

             # Use a custom install command which understands additional
             # commandline arguments
             'install': setup_helpers.AstropyInstall,

             'register': setup_helpers.AstropyRegister
             }

try:
    import bdist_mpkg
except ImportError:
    pass
else:
    # Use a custom command to build a dmg (on MacOS X)
    cmdclassd['bdist_dmg'] = setup_helpers.bdist_dmg

if setup_helpers.should_build_with_cython(release):
    from Cython.Distutils import build_ext
    # Builds Cython->C if in dev mode and Cython is present
    cmdclassd['build_ext'] = setup_helpers.wrap_build_ext(build_ext)
else:
    cmdclassd['build_ext'] = setup_helpers.wrap_build_ext()

if setup_helpers.HAVE_SPHINX:
    cmdclassd['build_sphinx'] = setup_helpers.AstropyBuildSphinx

# Set our custom command class mapping in setup_helpers, so that
# setup_helpers.get_distutils_option will use the custom classes.
setup_helpers.cmdclassd = cmdclassd

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
setup_helpers.update_package_files('astropy', extensions, package_data,
                                   packagenames, package_dirs)

# Currently the only entry points installed by Astropy are hooks to
# zest.releaser for doing Astropy's releases
entry_points = {}
for hook in [('releaser', 'middle'), ('postreleaser', 'before')]:
    hook_ep = 'zest.releaser.' + '.'.join(hook)
    hook_name = 'astropy.release.' + '.'.join(hook)
    hook_func = 'astropy.utils.release:' + '_'.join(hook)
    entry_points[hook_ep] = ['%s = %s' % (hook_name, hook_func)]


setup(name='astropy',
      version=version,
      description='Community-developed python astronomy tools',
      packages=packagenames,
      package_data=package_data,
      package_dir=package_dirs,
      ext_modules=extensions,
      scripts=scripts,
      requires=['numpy'],  # scipy not required, but strongly recommended
      install_requires=['numpy'],
      provides=['astropy'],
      author='The Astropy Developers',
      author_email='astropy.team@gmail.com',
      license='BSD',
      url='http://astropy.org',
      long_description=astropy.__doc__,
      download_url='%s/astropy-%s.tar.gz' % (download_base_url, version),
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

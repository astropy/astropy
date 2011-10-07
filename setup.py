#!/usr/bin/env python

# Use "distribute" - the setuptools fork that supports python 3.
from distribute_setup import use_setuptools
use_setuptools()

import os
import glob
from setuptools import setup, find_packages, Extension
from distutils.command.build_py import build_py as du_build_py
from warnings import warn

import astropy
from astropy.version import version as version_string
from astropy.version import release

# Use the find_packages tool to locate all packages and modules other than
# those that are in tests/
packages = find_packages(exclude=['tests'])

# Treat everything in scripts except README.rst as a script to be installed
scripts = glob.glob('scripts/*')
scripts.remove('scripts/README.rst')

# This dictionary stores the command classes used in setup below
cmdclassd = {}

# C extensions that are not Cython-based should be added here.
extensions = []

# Look for Cython files - compile with Cython if it is not a release
# and Cython is installed. Otherwise, use the .c files that live next
# to the Cython files.

try:
    import Cython
    have_cython = True
except ImportError:
    have_cython = False

pyxfiles = []
for  dirpath, dirnames, filenames in os.walk('astropy'):
    modbase = dirpath.replace(os.sep, '.')
    for fn in filenames:
        if fn.endswith('.pyx'):
            fullfn = os.path.join(dirpath, fn)
            extmod = modbase + '.' + fn[:-4]  # Package must match file name
            pyxfiles.append((extmod, fullfn))

if not release and have_cython:

    from Cython.Distutils import build_ext as cython_build_ext
    cmdclassd['build_ext'] = cython_build_ext

    # Add .pyx files
    for extmod, pyxfn in pyxfiles:
        extensions.append(Extension(extmod, [pyxfn]))

else:

    # Add .c files
    for extmod, pyxfn in pyxfiles:
        cfn = pyxfn[:-4] + '.c'
        if os.path.exists(cfn):
            extensions.append(Extension(extmod, [cfn]))
        else:
            warnstr = 'Could not find Cython-generated C extension ' + \
                 '{0} - The {1} module will be skipped.'.format(cfn, extmod)
            warn(warnstr)

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


# This customize build_py command replaces the source version.py with
# one that has a fixed version number.  This is necessary to support version
# numbers pulled from the git repository.

class astropy_build_py(du_build_py):

    def run(self):

        from astropy.version import _get_version_py_str

        res = du_build_py.run(self)

        versfile = os.path.join(self.build_lib, 'astropy', 'version.py')

        print 'Freezing version number to', versfile

        with open(versfile, 'w') as f:  # This overwrites the actual version.py
            f.write(_get_version_py_str())

        return res

cmdclassd['build_py'] = astropy_build_py

setup(name='AstroPy',
      version=version_string,
      description='Community-developed python astronomy tools',
      packages=packages,
      package_data={'astropy': ['data/*']},
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
      zip_safe=False
)

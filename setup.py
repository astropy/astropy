#!/usr/bin/env python
from __future__ import division

#Use "distribute" - the setuptools fork that supports python 3.
from distribute_setup import use_setuptools
use_setuptools()

from setuptools import setup,find_packages
from glob import glob

import astropy
from astropy.version import version as versionstring

NAME = 'Astropy'

SHORT_DESCRIPTION = 'Community-developed python astronomy tools'

#Take the long description from the root package docstring
LONG_DESCRIPTION = astropy.__doc__


#use the find_packages tool to locate all packages and modules other than
#those that are in "tests" 
astropyspkgs = find_packages(exclude=['tests'])

#treat everything in scripts except README.rst as a script to be installed
astropyscripts = glob('scripts/*')
astropyscripts.remove('scripts/README.rst')

#Implement a version of build_sphinx that automatically creates the docs/_build
#dir - this is needed because github won't create the _build dir because it 
#initially has no files
try:
    from sphinx.setup_command import BuildDoc
    
    class apy_build_sphinx(BuildDoc):
        def finalize_options(self):
            from os.path import isfile    
            from distutils.cmd import DistutilsOptionError
            
            if self.build_dir is not None:
                if isfile(self.build_dir):
                    raise DistutilsOptionError('Attempted to build_sphinx into a file '+self.build_dir)
                self.mkpath(self.build_dir)
            return BuildDoc.finalize_options(self)
            
except ImportError: #sphinx not present
    apy_build_sphinx = None
    
cmdclassd = {}
if apy_build_sphinx is not None:
    cmdclassd['build_sphinx'] = apy_build_sphinx
    
setup(name=NAME,
      version=versionstring,
      description=SHORT_DESCRIPTION,
      packages=astropypkgs,
      package_data={'astropy':['data/*']},
      scripts=astropyscripts,
      requires=['numpy','scipy'],
      install_requires=['numpy'],
      provides=['astropy'], 
      author='Astropy Team, Erik Tollerud, Thomas Robitaille, and Perry Greenfield',
      author_email='astropy.team@gmail.com',
      #license = '', #TODO: decide on a license
      url='http://astropy.org',
      long_description=LONG_DESCRIPTION,
      cmdclass = cmdclassd
)

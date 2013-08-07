#!/usr/bin/env python

from distutils.core import setup

setup(name='WCSAxesExperimental',
      version='0.0.0',
      description='Experimental WCSAxes implementatin',
      packages=['wcsaxes'],
      provides=['wcsaxes'],
      requires=['numpy', 'matplotlib', 'astropy'],
      keywords=['Scientific/Engineering'],
      classifiers=[
          "Development Status :: 1 - Planning",
          "Programming Language :: Python",
      ],
      )

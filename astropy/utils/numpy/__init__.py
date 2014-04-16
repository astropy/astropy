# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import, print_function

from distutils import version

from numpy import __version__ as numpy_version
NUMPY_VERSION = version.LooseVersion(numpy_version)

PR4622 = NUMPY_VERSION < version.LooseVersion('9.9.9')

from .lib.stride_tricks import *

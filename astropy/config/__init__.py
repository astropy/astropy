# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains configuration and setup utilities for the
`astropy` project. This includes all functionality related to the
affiliated package index.
"""

from .paths import *
from .configuration import *
from .data import *
from .affiliated import *
from .logging_helper import *
del log  # we remove the log here because it is imported in astropy/__init__.py

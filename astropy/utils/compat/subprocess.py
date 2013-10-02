"""
A replacement wrapper around the subprocess module that adds
check_output (which was only added to Python in 2.7.

Instead of importing subprocess, other modules should use this as follows::

    from astropy.utils.compat import subprocess

This module is safe to import from anywhere within astropy.
"""
from __future__ import absolute_import, print_function

import subprocess

# python2.7 and later provide a check_output method
if not hasattr(subprocess, 'check_output'):
    from ._subprocess_py2 import check_output

from subprocess import *

# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

from distutils.core import setup, Extension

import os
import sys

if __name__ == '__main__':
    sys.path.insert(0, os.getcwd())

    from astropy import wcs

    wcsapi_test_module = Extension(
        str('wcsapi_test'),
        include_dirs=[
            os.path.join(wcs.get_include(), 'astropy_wcs'),
            os.path.join(wcs.get_include(), 'wcslib')
        ],
        # Use the *full* name to the c file, since we can't change the cwd
        # during testing
        sources=[str(os.path.join(os.path.dirname(__file__), 'wcsapi_test.c'))])

    setup(
        name='wcsapi_test',
        ext_modules=[wcsapi_test_module])

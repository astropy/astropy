# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import sys

if __name__ == '__main__':
    astropy_path = sys.argv[-1]
    sys.argv = sys.argv[:-1]
    sys.path.insert(0, astropy_path)

    from astropy import wcs
    from astropy import setup_helpers
    from distutils.core import setup, Extension

    wcsapi_test_module = Extension(
        str('wcsapi_test'),
        include_dirs=[
            setup_helpers.get_numpy_include_path(),
            os.path.join(wcs.get_include(), 'astropy_wcs'),
            os.path.join(wcs.get_include(), 'wcslib')
        ],
        # Use the *full* name to the c file, since we can't change the cwd
        # during testing
        sources=[str(os.path.join(os.path.dirname(__file__), 'wcsapi_test.c'))])

    setup(
        name='wcsapi_test',
        ext_modules=[wcsapi_test_module])

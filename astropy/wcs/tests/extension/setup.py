# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys

if __name__ == '__main__':
    from setuptools import Extension, setup

    import numpy as np

    from astropy import wcs

    if sys.platform == 'win32':
        # These are written into wcsconfig.h, but that file is not
        # used by all parts of wcslib.
        define_macros = [
            ('YY_NO_UNISTD_H', None),
            ('_CRT_SECURE_NO_WARNINGS', None),
            ('_NO_OLDNAMES', None),  # for mingw32
            ('NO_OLDNAMES', None),  # for mingw64
            ('__STDC__', None)  # for MSVC
        ]
    else:
        define_macros = []

    # Below, we include a typical system include in case astropy was
    # installed with ASTROPY_USE_SYSTEM_WCSLIB.  Good enough for the
    # test, but a proper implementation would need to look better at
    # where wcslib might be stalled.
    wcsapi_test_module = Extension(
        'wcsapi_test',
        include_dirs=[
            np.get_include(),
            os.path.join(wcs.get_include(), 'astropy_wcs'),
            os.path.join(wcs.get_include(), 'wcslib'),
            os.path.join('/usr/include', 'wcslib'),
        ],
        # Use the *full* name to the c file, since we can't change the cwd
        # during testing
        sources=[str(os.path.join(os.path.dirname(__file__),
                                  'wcsapi_test.c'))],
        define_macros=define_macros)

    setup(
        name='wcsapi_test',
        ext_modules=[wcsapi_test_module])

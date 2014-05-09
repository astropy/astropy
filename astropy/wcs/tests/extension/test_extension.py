# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import subprocess
import sys


def test_wcsapi_extension(tmpdir):
    # Test that we can build a simple C extension with the astropy.wcs C API

    setup_path = os.path.join(os.path.dirname(__file__), 'setup.py')

    # Build the extension
    subprocess.check_call([
        sys.executable, setup_path,
        'install', '--install-lib={0}'.format(tmpdir)])

    code = """
    import sys
    sys.path.insert(0, "{0}")
    import wcsapi_test
    sys.exit(wcsapi_test.test())
    """

    code = code.strip().replace('\n', '; ').format(str(tmpdir))

    # Import and run the extension
    subprocess.check_call([sys.executable, '-c', code])

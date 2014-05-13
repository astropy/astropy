# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import subprocess
import sys


def test_wcsapi_extension(tmpdir):
    # Test that we can build a simple C extension with the astropy.wcs C API

    setup_path = os.path.dirname(__file__)

    env = os.environ.copy()
    env['PYTHONPATH'] = str(tmpdir) + ':' + env.get('PYTHONPATH', '')

    # Build the extension
    subprocess.check_call(
        [sys.executable, 'setup.py',
         'install', '--install-lib={0}'.format(tmpdir)],
        cwd=setup_path,
        env=env
    )

    code = """
    import sys
    import wcsapi_test
    sys.exit(wcsapi_test.test())
    """

    code = code.strip().replace('\n', '; ')

    # Import and run the extension
    subprocess.check_call(
        [sys.executable, '-c', code],
        env=env)

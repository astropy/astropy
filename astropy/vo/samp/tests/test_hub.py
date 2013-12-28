# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import time
import tempfile

from ....tests.helper import remote_data
from ... import samp


@remote_data
def test_SAMPHubServer():
    """Test that SAMPHub can be instantiated"""
    fileobj, lockfile = tempfile.mkstemp()
    samp.SAMPHubServer(web_profile=False, lockfile=lockfile)
    if os.path.exists(lockfile):
        os.remove(lockfile)


@remote_data
def test_SAMPHubServer_run():
    """Test that SAMPHub can be run"""
    fileobj, lockfile = tempfile.mkstemp()
    hub = samp.SAMPHubServer(web_profile=False, lockfile=lockfile)
    hub.start()
    time.sleep(1)
    hub.stop()
    if os.path.exists(lockfile):
        os.remove(lockfile)

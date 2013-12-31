# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import time
import tempfile

from ....tests.helper import remote_data

from ..hub import SAMPHubServer

# By default, tests should not use the internet.
from ..utils import ALLOW_INTERNET
ALLOW_INTERNET.set(False)


def test_SAMPHubServer():
    """Test that SAMPHub can be instantiated"""
    fileobj, lockfile = tempfile.mkstemp()
    SAMPHubServer(web_profile=False, lockfile=lockfile)
    if os.path.exists(lockfile):
        os.remove(lockfile)


def test_SAMPHubServer_run():
    """Test that SAMPHub can be run"""
    fileobj, lockfile = tempfile.mkstemp()
    hub = SAMPHubServer(web_profile=False, lockfile=lockfile)
    hub.start()
    time.sleep(1)
    hub.stop()
    if os.path.exists(lockfile):
        os.remove(lockfile)

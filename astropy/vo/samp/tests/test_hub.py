# Licensed under a 3-clause BSD style license - see LICENSE.rst

import time

from ..hub import SAMPHubServer

# By default, tests should not use the internet.
from ..utils import ALLOW_INTERNET

def setup_module(module):
    ALLOW_INTERNET.set(False)


def test_SAMPHubServer():
    """Test that SAMPHub can be instantiated"""
    SAMPHubServer(web_profile=False, mode='multiple', pool_size=1)


def test_SAMPHubServer_run():
    """Test that SAMPHub can be run"""
    hub = SAMPHubServer(web_profile=False, mode='multiple', pool_size=1)
    hub.start()
    time.sleep(1)
    hub.stop()

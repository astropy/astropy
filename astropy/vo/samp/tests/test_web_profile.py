"""
Test the web profile using Python classes that have been adapted to act like a
web client. We can only put a single test here because only one hub can run
with the web profile active, and the user might want to run the tests in
parallel.
"""

import time
import threading
import pickle
import tempfile

from .. import SAMPIntegratedClient, SAMPHubServer, SAMP_STATUS_OK
from .web_profile_test_helpers import (AlwaysApproveWebProfileDialog,
                                       SAMPIntegratedWebClient)


def write_output(mtype, params):
    filename = params['verification_file']
    f = open(filename, 'wb')
    pickle.dump(mtype, f)
    pickle.dump(params, f)
    f.close()


def check_output(mtype, params, timeout=None):
    filename = params['verification_file']
    start = time.time()
    while True:
        try:
            f = open(filename, 'rb')
            rec_mtype = pickle.load(f)
            rec_params = pickle.load(f)
            f.close()
            return
        except (IOError, EOFError):
            if timeout is not None and time.time() - start > timeout:
                raise Exception("Timeout while waiting for file: {0}".format(filename))

    assert rec_mtype == mtype
    assert rec_params == params
    

class Receiver(object):

    def __init__(self, client):
        self.client = client

    def receive_notification(self, private_key, sender_id, mtype, params, extra):
        write_output(mtype, params)

    def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
        self.receive_notification(private_key, sender_id, mtype, params, extra)
        self.client.reply(msg_id, {"samp.status": SAMP_STATUS_OK,
                                   "samp.result": {"txt": "test"}})

    
def temporary_filename():
    return tempfile.mkstemp()[1]
    
    
def test_web_profile():

    d = AlwaysApproveWebProfileDialog()
    t = threading.Thread(target=d.poll)
    t.start()

    lockfile = temporary_filename()

    h = SAMPHubServer(web_profile_dialog=d, lockfile=lockfile)
    h.start()

    c1 = SAMPIntegratedClient()
    c1.connect(hub=h)

    c2 = SAMPIntegratedWebClient()
    c2.connect()

    rec2 = Receiver(c1)
    c2.bind_receive_notification('samp.load.votable', rec2.receive_notification)
    c2.bind_receive_call('samp.load.votable', rec2.receive_call)

    # Test Notify

    params = {'verification_file':temporary_filename(),
              'parameter1':'abcde',
              'parameter2':1331}

    c1.notify(c2.get_public_id(), {'samp.mtype':'samp.load.votable', 'samp.params':params})

    check_output('samp.load.votable', params)

    # Test Call

    params = {'verification_file':temporary_filename(),
              'parameter1':'abcde',
              'parameter2':1331}

    c1.call(c2.get_public_id(), 'tag', {'samp.mtype':'samp.load.votable', 'samp.params':params})

    check_output('samp.load.votable', params)

    c1.disconnect()
    c2.disconnect()

    h.stop()
    d.stop()

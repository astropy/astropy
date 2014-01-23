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

from ....extern import six
from ....extern.six.moves.urllib.request import Request, urlopen
from ....utils.data import get_readable_fileobj

from .. import SAMPIntegratedClient, SAMPHubServer, SAMP_STATUS_OK
from .web_profile_test_helpers import (AlwaysApproveWebProfileDialog,
                                       SAMPIntegratedWebClient)
from ..web_profile import CROSS_DOMAIN, CLIENT_ACCESS_POLICY
from ..errors import SAMPHubError


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

    with SAMPHubServer(web_profile_dialog=d, lockfile=lockfile) as h:

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

        check_output('samp.load.votable', params, timeout=60)

        # Test Call

        params = {'verification_file':temporary_filename(),
                  'parameter1':'abcde',
                  'parameter2':1331}

        c1.call(c2.get_public_id(), 'tag', {'samp.mtype':'samp.load.votable', 'samp.params':params})

        check_output('samp.load.votable', params, timeout=60)

        c1.disconnect()
        c2.disconnect()
        d.stop()

        # Now check some additional queries to the server

        with get_readable_fileobj('http://localhost:21012/crossdomain.xml') as f:
            assert f.read() == CROSS_DOMAIN

        with get_readable_fileobj('http://localhost:21012/clientaccesspolicy.xml') as f:
            assert f.read() == CLIENT_ACCESS_POLICY

        # Check headers

        req = Request('http://127.0.0.1:21012/crossdomain.xml')
        req.add_header('Origin', 'test_web_profile')
        resp = urlopen(req)

        if six.PY2:
            assert resp.info().getheader('Access-Control-Allow-Origin') == 'test_web_profile'
            assert resp.info().getheader('Access-Control-Allow-Headers') == 'Content-Type'
            assert resp.info().getheader('Access-Control-Allow-Credentials') == 'true'
        else:
            assert resp.getheader('Access-Control-Allow-Origin') == 'test_web_profile'
            assert resp.getheader('Access-Control-Allow-Headers') == 'Content-Type'
            assert resp.getheader('Access-Control-Allow-Credentials') == 'true'

import time
from . import SAMPIntegratedClient
from ...table import Table


class Receiver(object):
    def __init__(self, client):
        self.client = client
        self.received = False
    def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
        self.params = params
        self.received = True
        self.client.reply(msg_id, {"samp.status": "samp.ok", "samp.result": {}})
    def receive_notification(self, private_key, sender_id, mtype, params, extra):
        self.params = params
        self.received = True


def wait_until_received(object, timeout=None, step=0.1):
    start_time = time.time()
    while not (hasattr(object, 'received') and object.received):
        time.sleep(step)
        if timeout is not None and time.time() - start_time > timeout:
            raise AttributeError("Timeout while waiting for message to be received".format(attribute, object))


def receive_table(timeout=None):

    client = SAMPIntegratedClient()
    client.connect()
    r = Receiver(client)
    client.bind_receive_call("table.load.votable", r.receive_call)
    client.bind_receive_notification("table.load.votable",
                                     r.receive_notification)

    try:
        wait_until_received(r, timeout=timeout)
        t = Table.read(r.params['url'])
    finally:
        client.disconnect()

    return t

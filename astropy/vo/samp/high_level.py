import time
import tempfile

from ...table import Table

from .constants import SAMP_ICON
from . import SAMPIntegratedClient
from ... import __version__
from ...extern.six.moves.urllib.parse import urljoin

__all__ = ['send_table', 'receive_table', 'list_table_clients']


def declare_metadata(client):

    metadata = {"samp.name": "astropy",
                "samp.description.text": "The Astropy Project",
                "samp.icon.url": client.client._xmlrpcAddr + "/samp/icon",
                "samp.documentation.url": "http://docs.astropy.org/en/stable/vo/samp",
                "author.name": "The Astropy Collaboration",
                "home.page": "http://www.astropy.org",
                "astropy.version": __version__
                }

    client.declare_metadata(metadata)


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


def list_table_clients():
    """
    List all SAMP clients that can read in tables
    """
    client = SAMPIntegratedClient()
    client.connect()
    declare_metadata(client)
    table_clients = {}
    for c in client.get_subscribed_clients('table.load.votable'):
        table_clients[c] = client.get_metadata(c)['samp.name']
    client.disconnect()
    return table_clients


def receive_table(timeout=None):

    client = SAMPIntegratedClient()
    client.connect()
    declare_metadata(client)
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


def send_table(table, destination='all', timeout='10'):

    # Write the table out to a temporary file
    output_file = tempfile.NamedTemporaryFile()
    table.write(output_file, format='votable')

    client = SAMPIntegratedClient()
    client.connect()
    declare_metadata(client)

    all_clients = client.get_registered_clients()

    message = {}
    message['samp.mtype'] = "table.load.votable"
    message['samp.params'] = {"url": urljoin('file:', output_file.name),
                              "name": "Table from Astropy"}

    if destination == 'all':
        for c in client.get_subscribed_clients('table.load.votable'):
            client.call_and_wait(c, message, timeout=timeout)
    else:
        client.call_and_wait(destination, message, timeout=timeout)

    client.disconnect()

    output_file.close()

import time
import tempfile

import numpy as np

from ...table import Table
from ...nddata import NDData
from ...io import fits

from .constants import SAMP_ICON
from . import SAMPIntegratedClient
from ... import __version__
from ...extern.six.moves.urllib.parse import urljoin

__all__ = ['send', 'receive', 'list_clients']


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
        self.mtype = mtype
        self.params = params
        self.received = True
        self.client.reply(msg_id, {"samp.status": "samp.ok", "samp.result": {}})
    def receive_notification(self, private_key, sender_id, mtype, params, extra):
        self.mtype = mtype
        self.params = params
        self.received = True


def wait_until_received(object, timeout=None, step=0.1):
    start_time = time.time()
    while not (hasattr(object, 'received') and object.received):
        time.sleep(step)
        if timeout is not None and time.time() - start_time > timeout:
            raise AttributeError("Timeout while waiting for message to be received".format(attribute, object))


def list_clients():
    """
    List all SAMP clients that can read in tables
    """
    client = SAMPIntegratedClient()
    client.connect()
    declare_metadata(client)
    table_clients = {}
    for c in client.get_subscribed_clients('table.load.votable'):
        table_clients[c] = client.get_metadata(c)
    client.disconnect()
    return table_clients


def receive(timeout=None):
    """
    Receive data from SAMP clients.

    Parameters
    ----------
    timeout : int, optional
        How long to wait for before giving up
    """

    client = SAMPIntegratedClient()
    client.connect()
    declare_metadata(client)
    r = Receiver(client)

    for mtype in ['table.load.votable', 'table.load.fits', 'image.load.fits']:
        client.bind_receive_call(mtype, r.receive_call)
        client.bind_receive_notification(mtype, r.receive_notification)

    try:
        wait_until_received(r, timeout=timeout)
        if r.mtype.startswith('table'):
            data = Table.read(r.params['url'])
        else:
            data = fits.open(r.params['url'])
    finally:
        client.disconnect()

    return data


def send(data, name, destination='all', timeout=10):
    """
    Send data to SAMP clients.

    Parameters
    ----------
    data : `~astropy.table.table.Table` or `~astropy.nddata.nddata.NDData` or `~numpy.ndarray` or `~astropy.io.fits.PrimaryHDU` or `~astropy.io.fits.ImageHDU` or `~astropy.io.fits.BinTableHDU` or `~astropy.io.fits.TableHDU`
        The data to send over SAMP
    name : str, optional
        The name of the dataset to use in other SAMP clients
    destination : str, optional
        The client to send the data to. By default, the data is broadcast to
        all SAMP clients. You can find the full list of available clients, use
        the :func:`~astropy.vo.samp.high_level.list_clients` function. As a
        convenience, you can also use ``'ds9'``, ``'topcat'``, and ``aladin'``
        and :func:`~astropy.vo.samp.high_level.send` will try and identify the
        correct client.
    timeout : int, optional
        The timeout for the request.
    """

    message = {}
    output_file = tempfile.NamedTemporaryFile()

    if isinstance(data, Table):

        data.write(output_file, format='votable')
        message['samp.mtype'] = "table.load.votable"

    elif isinstance(data, NDData):

        data.write(output_file, format='fits')
        message['samp.mtype'] = "image.load.fits"

    elif isinstance(data, np.ndarray):

        if data.dtype.fields is None:
            fits.writeto(output_file, data)
            message['samp.mtype'] = "image.load.fits"
        else:
            data = Table(data)
            data.write(output_file, format='votable')
            message['samp.mtype'] = "table.load.votable"

    elif isinstance(data, (ImageHDU, PrimaryHDU)):

        data.writeto(output_file)
        message['samp.mtype'] = "image.load.fits"

    elif isinstance(data, (BinTableHDU, TableHDU)):

        data.writeto(output_file)
        message['samp.mtype'] = "table.load.fits"

    else:

        raise TypeError("Unrecognized data type: {0}".format(type(data)))

    message['samp.params'] = {"url": urljoin('file:', output_file.name),
                              "name": name}

    client = SAMPIntegratedClient()
    client.connect()
    declare_metadata(client)

    if destination == 'all':
        for c in client.get_subscribed_clients(message['samp.mtype']):
            client.call_and_wait(c, message, timeout=str(timeout))
    elif destination in ['ds9', 'topcat', 'aladin']:
        clients = list_clients()
        for client_id in clients:
            name = clients[client_id]['samp.name']
            if destination in name.lower():
                client.call_and_wait(client_id, message, timeout=str(timeout))
    else:
        client.call_and_wait(destination, message, timeout=str(timeout))

    client.disconnect()

    output_file.close()

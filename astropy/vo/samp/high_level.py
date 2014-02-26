import time
import tempfile

import numpy as np

from ...table import Table
from ...nddata import NDData
from ...io import fits

from ...utils.exceptions import TimeoutError
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
            raise TimeoutError("Timeout while waiting for message to be received")


def list_clients():
    """
    List all SAMP clients that can read in tables

    Returns
    -------
    client_table : `~astropy.table.table.Table`
        Table of clients and summary of what they support
    """
    client = SAMPIntegratedClient()
    client.connect()
    declare_metadata(client)
    table_clients = []
    for c in client.get_subscribed_clients('table.load.votable'):
        meta = client.get_metadata(c)
        subs = client.get_subscriptions(c)
        alias = 'none'
        for a in ['ds9', 'topcat', 'aladin']:
            if a in meta['samp.name'].lower():
                alias = a
        supports_table = 'Yes' if any(s.startswith('table.load') for s in subs) else 'No'
        supports_image = 'Yes' if any(s.startswith('image.load') for s in subs) else 'No'
        table_clients.append((c, alias, meta['samp.name'], supports_table, supports_image))
    client.disconnect()
    table_clients = Table(zip(*table_clients),
                          names=['id', 'alias', 'name', 'table', 'image'], masked=True)
    table_clients['alias'].mask = table_clients['alias'] == 'none'
    table_clients.sort('id')
    return table_clients


def receive(timeout=None, hub=None):
    """
    Receive data from SAMP clients.

    Parameters
    ----------
    timeout : int, optional
        How long to wait for before giving up
    hub : `~astropy.vo.samp.hub.SAMPHubServer`, optional
        The hub to receive the data through

    Returns
    -------
    data : `~astropy.table.table.Table` or `~astropy.io.fits.HDUList`
        The data received over SAMP.
    """

    client = SAMPIntegratedClient()
    client.connect(hub=hub)
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


def send(data, name, destination='all', timeout=10, hub=None):
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
    hub : `~astropy.vo.samp.hub.SAMPHubServer`, optional
        The hub to send the data through
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

    elif isinstance(data, (fits.ImageHDU, fits.PrimaryHDU)):

        data.writeto(output_file)
        message['samp.mtype'] = "image.load.fits"

    elif isinstance(data, (fits.BinTableHDU, fits.TableHDU)):

        data.writeto(output_file)
        message['samp.mtype'] = "table.load.fits"

    else:

        raise TypeError("Unrecognized data type: {0}".format(type(data)))

    message['samp.params'] = {"url": urljoin('file:', output_file.name),
                              "name": name}

    client = SAMPIntegratedClient()
    client.connect(hub=hub)
    declare_metadata(client)

    if destination == 'all':
        for c in client.get_subscribed_clients(message['samp.mtype']):
            client.call_and_wait(c, message, timeout=str(timeout))
    elif destination in ['ds9', 'topcat', 'aladin']:
        clients = list_clients()
        for target_client in clients:
            name = target_client['name']
            if destination in name.lower():
                client.call_and_wait(str(target_client['id']), message,
                                     timeout=str(timeout))
    else:
        client.call_and_wait(destination, message, timeout=str(timeout))

    client.disconnect()

    output_file.close()

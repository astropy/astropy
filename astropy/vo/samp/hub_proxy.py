# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy

from ...extern.six.moves import xmlrpc_client as xmlrpc

from .errors import SAMPHubError
from .utils import ServerProxyPool
from .lockfile_helpers import get_main_running_hub

from .constants import SSL_SUPPORT

if SSL_SUPPORT:
    from .ssl_utils import SafeTransport


__all__ = ['SAMPHubProxy']


class SAMPHubProxy(object):
    """
    Proxy class to simplify the client interaction with a SAMP hub (via the
    standard profile).
    """

    def __init__(self):
        self.proxy = None
        self._connected = False

    @property
    def is_connected(self):
        """
        Whether the hub proxy is currently connected to a hub.
        """
        return self._connected

    def connect(self, hub=None, hub_params=None,
                key_file=None, cert_file=None, cert_reqs=0,
                ca_certs=None, ssl_version=None, pool_size=20):
        """
        Connect to the current SAMP Hub.

        Parameters
        ----------
        hub : `~astropy.vo.samp.SAMPHubServer`, optional
            The hub to connect to.

        hub_params : dict, optional
            Optional dictionary containing the lock-file content of the Hub
            with which to connect. This dictionary has the form
            ``{<token-name>: <token-string>, ...}``.

        key_file : str, optional
            The path to a file containing the private key for SSL connections.
            If the certificate file (``cert_file``) contains the private key,
            then ``key_file`` can be omitted.

        cert_file : str, optional
            The path to a file which contains a certificate to be used to
            identify the local side of the secure connection.

        cert_reqs : int, optional
            Whether a certificate is required from the server side of the
            connection, and whether it will be validated if provided. It must
            be one of the three values `ssl.CERT_NONE` (certificates ignored),
            `ssl.CERT_OPTIONAL` (not required, but validated if provided), or
            `ssl.CERT_REQUIRED` (required and validated). If the value of this
            parameter is not `ssl.CERT_NONE`, then the ``ca_certs`` parameter
            must point to a file of CA certificates.

        ca_certs : str, optional
            The path to a file containing a set of concatenated "Certification
            Authority" certificates, which are used to validate the
            certificate passed from the Hub end of the connection.

        ssl_version : int, optional
            Which version of the SSL protocol to use. Typically, the
            server chooses a particular protocol version, and the
            client must adapt to the server's choice. Most of the
            versions are not interoperable with the other versions. If
            not specified, the default SSL version is taken from the
            default in the installed version of the Python standard
            `ssl` library.  See the `ssl` documentation for more
            information.

        pool_size : int, optional
            The number of socket connections opened to communicate with the
            Hub.
        """

        self._connected = False
        self.lockfile = {}

        if hub is not None and hub_params is not None:
            raise ValueError("Cannot specify both hub and hub_params")

        if hub_params is None:

            if hub is not None:
                if not hub.is_running:
                    raise SAMPHubError("Hub is not running")
                else:
                    hub_params = hub.params
            else:
                hub_params = get_main_running_hub()

        try:

            url = hub_params["samp.hub.xmlrpc.url"].replace("\\", "")

            if SSL_SUPPORT and url[0:5] == "https":

                transport = SafeTransport(key_file, cert_file, cert_reqs,
                                          ca_certs, ssl_version)

                self.proxy = ServerProxyPool(pool_size, xmlrpc.ServerProxy,
                                             url, transport=transport,
                                             allow_none=1)

            else:

                self.proxy = ServerProxyPool(pool_size, xmlrpc.ServerProxy,
                                             url, allow_none=1)

            self.ping()

            self.lockfile = copy.deepcopy(hub_params)
            self._connected = True

        except xmlrpc.ProtocolError as p:
            # 401 Unauthorized
            if p.errcode == 401:
                raise SAMPHubError("Unauthorized access. Basic Authentication "
                                   "required or failed.")
            else:
                raise SAMPHubError("Protocol Error {}: {}".format(p.errcode,
                                                                  p.errmsg))

    def disconnect(self):
        """
        Disconnect from the current SAMP Hub.
        """
        self.proxy = None
        self._connected = False
        self.lockfile = {}

    def server_close(self):
        self.proxy.server_close()

    @property
    def _samp_hub(self):
        """
        Property to abstract away the path to the hub, which allows this class
        to be used for other profiles.
        """
        return self.proxy.samp.hub

    def ping(self):
        """
        Proxy to ``ping`` SAMP Hub method (Standard Profile only).
        """
        return self._samp_hub.ping()

    def set_xmlrpc_callback(self, private_key, xmlrpc_addr):
        """
        Proxy to ``setXmlrpcCallback`` SAMP Hub method (Standard Profile only).
        """
        return self._samp_hub.setXmlrpcCallback(private_key, xmlrpc_addr)

    def register(self, secret):
        """
        Proxy to ``register`` SAMP Hub method.
        """
        return self._samp_hub.register(secret)

    def unregister(self, private_key):
        """
        Proxy to ``unregister`` SAMP Hub method.
        """
        return self._samp_hub.unregister(private_key)

    def declare_metadata(self, private_key, metadata):
        """
        Proxy to ``declareMetadata`` SAMP Hub method.
        """
        return self._samp_hub.declareMetadata(private_key, metadata)

    def get_metadata(self, private_key, client_id):
        """
        Proxy to ``getMetadata`` SAMP Hub method.
        """
        return self._samp_hub.getMetadata(private_key, client_id)

    def declare_subscriptions(self, private_key, subscriptions):
        """
        Proxy to ``declareSubscriptions`` SAMP Hub method.
        """
        return self._samp_hub.declareSubscriptions(private_key, subscriptions)

    def get_subscriptions(self, private_key, client_id):
        """
        Proxy to ``getSubscriptions`` SAMP Hub method.
        """
        return self._samp_hub.getSubscriptions(private_key, client_id)

    def get_registered_clients(self, private_key):
        """
        Proxy to ``getRegisteredClients`` SAMP Hub method.
        """
        return self._samp_hub.getRegisteredClients(private_key)

    def get_subscribed_clients(self, private_key, mtype):
        """
        Proxy to ``getSubscribedClients`` SAMP Hub method.
        """
        return self._samp_hub.getSubscribedClients(private_key, mtype)

    def notify(self, private_key, recipient_id, message):
        """
        Proxy to ``notify`` SAMP Hub method.
        """
        return self._samp_hub.notify(private_key, recipient_id, message)

    def notify_all(self, private_key, message):
        """
        Proxy to ``notifyAll`` SAMP Hub method.
        """
        return self._samp_hub.notifyAll(private_key, message)

    def call(self, private_key, recipient_id, msg_tag, message):
        """
        Proxy to ``call`` SAMP Hub method.
        """
        return self._samp_hub.call(private_key, recipient_id, msg_tag, message)

    def call_all(self, private_key, msg_tag, message):
        """
        Proxy to ``callAll`` SAMP Hub method.
        """
        return self._samp_hub.callAll(private_key, msg_tag, message)

    def call_and_wait(self, private_key, recipient_id, message, timeout):
        """
        Proxy to ``callAndWait`` SAMP Hub method.
        """
        return self._samp_hub.callAndWait(private_key, recipient_id, message,
                                          timeout)

    def reply(self, private_key, msg_id, response):
        """
        Proxy to ``reply`` SAMP Hub method.
        """
        return self._samp_hub.reply(private_key, msg_id, response)

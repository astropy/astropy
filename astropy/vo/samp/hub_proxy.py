# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import re
import copy
import sys
import traceback

from ...extern.six import StringIO

from .hub import SAMPHubServer
from .errors import SAMPHubError
from .utils import ServerProxyPool, SafeTransport, xmlrpc

try:
    import ssl
except ImportError:
    SSL_SUPPORT = False
else:
    SSL_SUPPORT = True

__all__ = ['SAMPHubProxy']


class SAMPHubProxy(object):
    """
    Proxy class to simplify the client interaction with a SAMP hub.
    """

    def __init__(self):
        self.proxy = None
        self._connected = False
        self.lockfile = {}

    @property
    def is_connected(self):
        """
        Testing method to verify the proxy connection with a running hub.

        Returns
        -------
        isConnected : bool
            True if the proxy is connected to a Hub, False otherwise
        """
        return self._connected

    @staticmethod
    def get_running_hubs():
        """
        Return a dictionary containing the lock-file contents of all the currently
        running hubs (single and/or multiple mode).

        The dictionary format is:
        `{<lock-file>: {<token-name>: <token-string>, ...}, ...}`

        where `{<lock-file>}` is the lock-file name, `{<token-name>}` and `{<token-string>}`
        are the lock-file tokens (name and content).

        Returns
        -------
        running_hubs : dict
            Lock-file contents of all the currently running hubs.
        """

        hubs = {}
        lockfilename = ""

        # HUB SINGLE INSTANCE MODE

        # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
        if "SAMP_HUB" in os.environ:
            # For the time being I assume just the std profile supported.
            if os.environ["SAMP_HUB"].startswith("std-lockurl:"):
                lockfilename = os.environ["SAMP_HUB"][len("std-lockurl:"):]
        else:
            if "HOME" in os.environ:
                # UNIX
                lockfilename = os.path.join(os.environ["HOME"], ".samp")
            else:
                # Windows
                lockfilename = os.path.join(os.environ["USERPROFILE"], ".samp")

        hub_is_running, lockfiledict = SAMPHubServer.check_running_hub(lockfilename)

        if hub_is_running:
            hubs[lockfilename] = lockfiledict

        # HUB MULTIPLE INSTANCE MODE

        lockfiledir = ""

        if "HOME" in os.environ:
            # UNIX
            lockfiledir = os.path.join(os.environ["HOME"], ".samp-1")
        else:
            # Windows
            lockfiledir = os.path.join(os.environ["USERPROFILE"], ".samp-1")

        if os.path.isdir(lockfiledir):
            for filename in os.listdir(lockfiledir):
                if re.match('samp\\-hub\\-\d+\\-\d+', filename) != None:
                    lockfilename = os.path.join(lockfiledir, filename)
                    hub_is_running, lockfiledict = SAMPHubServer.check_running_hub(lockfilename)
                    if hub_is_running:
                        hubs[lockfilename] = lockfiledict

        return hubs

    def connect(self, hub_params=None, user=None, password=None,
                key_file=None, cert_file=None, cert_reqs=0,
                ca_certs=None, ssl_version=1, pool_size=20):
        """
        Connect to the current SAMP Hub.

        If a SAMP Hub is not running or refuses the connection, then a `SAMPHubError` is raised.

        Parameters
        ----------
        hub_params : dict
            Optional dictionary containing the lock-file content of the Hub with which to connect.
            This dictionary has the form `{<token-name>: <token-string>, ...}`.

        user : str
            In case of Basic Authenticated connections, `user` specifies the user name.

        password : str
            In case of Basic Authenticated connections, `password` specifies the user password.

        key_file : str
            Set the file containing the private key for SSL connections. If the
            certificate file (`certfile`) contains the private key, then `keyfile` can be omitted.

        cert_file : str
            Specify the file which contains a certificate to be used to identify the
            local side of the secure connection.

        cert_reqs : int
            The parameter `cert_reqs` specifies whether a certificate is required
            from the server side of the connection, and whether it will be validated if provided. It
            must be one of the three values `ssl.CERT_NONE` (certificates ignored), `ssl.CERT_OPTIONAL`
            (not required, but validated if provided), or `ssl.CERT_REQUIRED` (required and validated).
            If the value of this parameter is not `ssl.CERT_NONE`, then the `ca_certs` parameter must
            point to a file of CA certificates.

        ca_certs : str
            The `ca_certs` file contains a set of concatenated "Certification Authority"
            certificates, which are used to validate the certificate passed from the server end of the
            connection.

        ssl_version : int
            The `ssl_version` option specifies which version of the SSL protocol to use.
            Typically, the server chooses a particular protocol version, and the client must adapt to the
            server's choice. Most of the versions are not interoperable with the other versions. If not
            specified the default SSL version is `ssl.PROTOCOL_SSLv3`. This version provides the most
            compatibility with other versions server side. Other SSL protocol versions are:
            `ssl.PROTOCOL_SSLv2`, `ssl.PROTOCOL_SSLv23` and `ssl.PROTOCOL_TLSv1`.

        pool_size : int
            The number of socket connections opened to communicate with the Hub.
        """

        self._connected = False
        self.lockfile = {}

        if hub_params == None:
            hubs = SAMPHubProxy.get_running_hubs()
            if len(hubs.keys()) > 0:
                # Use Single instance hub by default
                lockfilename = ""
                # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
                if "SAMP_HUB" in os.environ:
                    # For the time being I assume just the std profile supported.
                    if os.environ["SAMP_HUB"].startswith("std-lockurl:"):
                        lockfilename = os.environ["SAMP_HUB"][len("std-lockurl:"):]
                    else:
                        raise SAMPHubError("SAMP Hub profile not supported.")
                else:
                    if "HOME" in os.environ:
                        # UNIX
                        lockfilename = os.path.join(os.environ["HOME"], ".samp")
                    else:
                        # Windows
                        lockfilename = os.path.join(os.environ["USERPROFILE"], ".samp")
                hub_params = hubs[lockfilename]
            else:
                raise SAMPHubError("Unable to find a running SAMP Hub.")

        try:

            url = hub_params["samp.hub.xmlrpc.url"].replace("\\", "")

            # URL formatting for Basic Authentication parameters
            if user != None and password != None:
                trans, addr = url.split("://")
                url = "%s://%s:%s@%s" % (trans, user, password, addr)

            if SSL_SUPPORT and url[0:5] == "https":
                self.proxy = ServerProxyPool(pool_size, xmlrpc.ServerProxy,
                                             url, transport=SafeTransport(key_file, cert_file, cert_reqs,
                                                                          ca_certs, ssl_version),
                                             allow_none=1)
            else:
                self.proxy = ServerProxyPool(pool_size, xmlrpc.ServerProxy, url, allow_none=1)

            self.proxy.samp.hub.ping()

            self.lockfile = copy.deepcopy(hub_params)
            self._connected = True

        except xmlrpc.ProtocolError as p:
            # 401 Unauthorized
            if p.errcode == 401:
                raise SAMPHubError("Unauthorized access. Basic Authentication required or failed.")
            else:
                raise SAMPHubError("Protocol Error %d: %s" % (p.errcode, p.errmsg))
        except:
            err = StringIO()
            traceback.print_exc(file=err)
            txt = err.getvalue()
            if SSL_SUPPORT:
                if sys.exc_info()[0] == ssl.SSLError:
                    raise SAMPHubError("SSL Error: %s" % sys.exc_info()[1])
                else:
                    raise SAMPHubError("SAMP Hub connection refused.\n " + txt)
            else:
                raise SAMPHubError("SAMP Hub connection refused.\n" + txt)

    def disconnect(self):
        """Disconnect from the current SAMP Hub."""
        self.proxy = None
        self._connected = False
        self.lockfile = {}

    def ping(self):
        """Proxy to `ping` SAMP Hub method (Standard Profile only)."""
        return self.proxy.samp.hub.ping()

    def set_xmlrpc_callback(self, private_key, xmlrpc_addr):
        """Proxy to `setXmlrpcCallback` SAMP Hub method (Standard Profile only)."""
        return self.proxy.samp.hub.setXmlrpcCallback(private_key, xmlrpc_addr)

    def register(self, secret):
        """Proxy to `register` SAMP Hub method."""
        return self.proxy.samp.hub.register(secret)

    def unregister(self, private_key):
        """Proxy to `unregister` SAMP Hub method."""
        return self.proxy.samp.hub.unregister(private_key)

    def declare_metadata(self, private_key, metadata):
        """Proxy to `declareMetadata` SAMP Hub method."""
        return self.proxy.samp.hub.declareMetadata(private_key, metadata)

    def get_metadata(self, private_key, client_id):
        """Proxy to `getMetadata` SAMP Hub method."""
        return self.proxy.samp.hub.getMetadata(private_key, client_id)

    def declare_subscriptions(self, private_key, subscriptions):
        """Proxy to `declareSubscriptions` SAMP Hub method."""
        return self.proxy.samp.hub.declareSubscriptions(private_key, subscriptions)

    def get_subscriptions(self, private_key, client_id):
        """Proxy to `getSubscriptions` SAMP Hub method."""
        return self.proxy.samp.hub.getSubscriptions(private_key, client_id)

    def get_registered_clients(self, private_key):
        """Proxy to `getRegisteredClients` SAMP Hub method."""
        return self.proxy.samp.hub.getRegisteredClients(private_key)

    def get_subscribed_clients(self, private_key, mtype):
        """Proxy to `getSubscribedClients` SAMP Hub method."""
        return self.proxy.samp.hub.getSubscribedClients(private_key, mtype)

    def notify(self, private_key, recipient_id, message):
        """Proxy to `notify` SAMP Hub method."""
        # Add user in Basic Authentication case
        return self.proxy.samp.hub.notify(private_key, recipient_id, message)

    def notify_all(self, private_key, message):
        """Proxy to `notifyAll` SAMP Hub method."""
        return self.proxy.samp.hub.notifyAll(private_key, message)

    def call(self, private_key, recipient_id, msg_tag, message):
        """Proxy to `call` SAMP Hub method."""
        return self.proxy.samp.hub.call(private_key, recipient_id, msg_tag, message)

    def call_all(self, private_key, msg_tag, message):
        """Proxy to `callAll` SAMP Hub method."""
        return self.proxy.samp.hub.callAll(private_key, msg_tag, message)

    def call_and_wait(self, private_key, recipient_id, message, timeout):
        """Proxy to `callAndWait` SAMP Hub method.

        If timeout expires a `SAMPProxyError` instance is raised.
        """
        return self.proxy.samp.hub.callAndWait(private_key, recipient_id, message, timeout)

    def reply(self, private_key, msg_id, response):
        """Proxy to `reply` SAMP Hub method."""
        return self.proxy.samp.hub.reply(private_key, msg_id, response)

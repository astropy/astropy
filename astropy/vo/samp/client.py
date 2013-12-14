# Licensed under a 3-clause BSD style license - see LICENSE.rst
import copy
import inspect
import os
import re
import select
import socket
import sys
import threading
import traceback

try:
    import ssl
except ImportError:
    SSL_SUPPORT = False
else:
    SSL_SUPPORT = True

from ...extern.six import StringIO

from .constants import SAMP_STATUS_OK, SAMP_STATUS_WARNING
from .constants import _THREAD_STARTED_COUNT
from .errors import SAMPHubError, SAMPClientError, SAMPProxyError
from .utils import internet_on, ServerProxyPool, SafeTransport
from .utils import SecureXMLRPCServer, ThreadingXMLRPCServer
from .hub import SAMPHubServer

# Python 2 / 3 dependent imports ... for now get from utils to avoid code duplication
from .utils import xmlrpc


__all__ = ['SAMPHubProxy', 'SAMPIntegratedClient', 'SAMPClient']

__doctest_skip__ = ['SAMPIntegratedClient.*']


class SAMPHubProxy(object):
    """
    Proxy class to simplify the client interaction with a SAMP hub.
    """

    def __init__(self):
        self.proxy = None
        self._connected = False
        self.lockfile = {}

    def isConnected(self):
        """
        Testing method to verify the proxy connection with a running hub.

        Returns
        -------
        isConnected : bool
            True if the proxy is connected to a Hub, False otherwise
        """
        return self._connected

    @staticmethod
    def getRunningHubs():
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

        hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)

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
                    hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)
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
            hubs = SAMPHubProxy.getRunningHubs()
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

    def setXmlrpcCallback(self, private_key, xmlrpc_addr):
        """Proxy to `setXmlrpcCallback` SAMP Hub method (Standard Profile only)."""
        return self.proxy.samp.hub.setXmlrpcCallback(private_key, xmlrpc_addr)

    def register(self, secret):
        """Proxy to `register` SAMP Hub method."""
        return self.proxy.samp.hub.register(secret)

    def unregister(self, private_key):
        """Proxy to `unregister` SAMP Hub method."""
        return self.proxy.samp.hub.unregister(private_key)

    def declareMetadata(self, private_key, metadata):
        """Proxy to `declareMetadata` SAMP Hub method."""
        return self.proxy.samp.hub.declareMetadata(private_key, metadata)

    def getMetadata(self, private_key, client_id):
        """Proxy to `getMetadata` SAMP Hub method."""
        return self.proxy.samp.hub.getMetadata(private_key, client_id)

    def declareSubscriptions(self, private_key, subscriptions):
        """Proxy to `declareSubscriptions` SAMP Hub method."""
        return self.proxy.samp.hub.declareSubscriptions(private_key, subscriptions)

    def getSubscriptions(self, private_key, client_id):
        """Proxy to `getSubscriptions` SAMP Hub method."""
        return self.proxy.samp.hub.getSubscriptions(private_key, client_id)

    def getRegisteredClients(self, private_key):
        """Proxy to `getRegisteredClients` SAMP Hub method."""
        return self.proxy.samp.hub.getRegisteredClients(private_key)

    def getSubscribedClients(self, private_key, mtype):
        """Proxy to `getSubscribedClients` SAMP Hub method."""
        return self.proxy.samp.hub.getSubscribedClients(private_key, mtype)

    def notify(self, private_key, recipient_id, message):
        """Proxy to `notify` SAMP Hub method."""
        # Add user in Basic Authentication case
        return self.proxy.samp.hub.notify(private_key, recipient_id, message)

    def notifyAll(self, private_key, message):
        """Proxy to `notifyAll` SAMP Hub method."""
        return self.proxy.samp.hub.notifyAll(private_key, message)

    def call(self, private_key, recipient_id, msg_tag, message):
        """Proxy to `call` SAMP Hub method."""
        return self.proxy.samp.hub.call(private_key, recipient_id, msg_tag, message)

    def callAll(self, private_key, msg_tag, message):
        """Proxy to `callAll` SAMP Hub method."""
        return self.proxy.samp.hub.callAll(private_key, msg_tag, message)

    def callAndWait(self, private_key, recipient_id, message, timeout):
        """Proxy to `callAndWait` SAMP Hub method.

        If timeout expires a `SAMPProxyError` instance is raised.
        """
        return self.proxy.samp.hub.callAndWait(private_key, recipient_id, message, timeout)

    def reply(self, private_key, msg_id, response):
        """Proxy to `reply` SAMP Hub method."""
        return self.proxy.samp.hub.reply(private_key, msg_id, response)


class SAMPIntegratedClient(object):
    """
    A Simple SAMP client.

    This class is meant to simplify the client usage providing
    a proxy class that merges the `SAMPClient` and `SAMPHubProxy`
    functionalities in a simplified API.

    Parameters
    ----------
    name: str, optional
        Client name (corresponding to `samp.name` metadata keyword).

    description : str, optional
        Client description (corresponding to `samp.description.text` metadata keyword).

    metadata : dict, optional
        Client application metadata in the standard SAMP format. If present, `samp.name`
        keyword and `samp.description.text` keyword are overwritten by the parameters
        `name` and `description`.

    addr : string, optional
        Listening address (or IP).

    port : int, optional
        Listening XML-RPC server socket port.

    https : bool
        Set the callable client running on a Secure Sockets Layer connection (HTTPS).
        By default SSL is desabled.

    key_file : str
        Set the file containing the private key for SSL connections. If the
        certificate file (`cert_file`) contains the private key, then `key_file` can be omitted.

    cert_file : str
        Specify the file which contains a certificate to be used to identify the
        local side of the secure connection.

    cert_reqs : int
        The parameter `cert_reqs` specifies whether a certificate is required
        from the Hub side of the connection, and whether it will be validated if provided. It
        must be one of the three values `ssl.CERT_NONE` (certificates ignored), `ssl.CERT_OPTIONAL`
        (not required, but validated if provided), or `ssl.CERT_REQUIRED` (required and validated).
        If the value of this parameter is not `ssl.CERT_NONE`, then the `ca_certs` parameter must
        point to a file of CA certificates.

    ca_certs : str
        The `ca_certs` file contains a set of concatenated "Certification Authority"
        certificates, which are used to validate the certificate passed from the Hub end of the
        connection.

    ssl_version : int
        The `ssl_version` option specifies which version of the SSL protocol to use.
        Typically, the server chooses a particular protocol version, and the client must adapt to the
        server's choice. Most of the versions are not interoperable with the other versions. If not
        specified the default SSL version is `ssl.PROTOCOL_SSLv23`. This version provides the most
        compatibility with other versions Hub side. Other SSL protocol versions are:
        `ssl.PROTOCOL_SSLv2`, `ssl.PROTOCOL_SSLv3` and `ssl.PROTOCOL_TLSv1`.

    callable : bool
        Is the client callable?
    """

    def __init__(self, name=None, description=None, metadata=None,
                 addr=None, port=0, https=False, key_file=None,
                 cert_file=None, cert_reqs=0, ca_certs=None, ssl_version=2,
                 callable=True):
        self.hub = SAMPHubProxy()

        self.client = SAMPClient(self.hub, name, description, metadata, addr, port,
                                 https, key_file, cert_file, cert_reqs,
                                 ca_certs, ssl_version, callable)

    def __del__(self):
        try:
            self.disconnect()
        except:
            pass

    # GENERAL
    def isConnected(self):
        """
        Testing method to verify the client connection with a running Hub.

        Returns
        -------
        isConnected : bool
            True if the client is connected to a Hub, False otherwise.
        """
        return self.hub.isConnected() & self.client.isRunning()

    def connect(self, hub_params=None, user=None, password=None,
                key_file=None, cert_file=None, cert_reqs=0,
                ca_certs=None, ssl_version=1, pool_size=20):
        """
        Connect with the current or specified SAMP Hub, start and register the client.

        If a SAMP Hub is not running or refuses the connection, a `SAMPHubError` is raised.

        Parameters
        ----------
        hub_params : dict
            Optional dictionary containig the lock-file content of the Hub
            with which to connect. This dictionary has the form `{<token-name>: <token-string>, ...}`.

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
            Number of socket connections opened to communicate with the Hub.
        """
        self.hub.connect(hub_params, user, password, key_file, cert_file, cert_reqs,
                         ca_certs, ssl_version, pool_size)
        self.client.start()
        self.client.register()

    def disconnect(self):
        """
        Unregister the client from the current SAMP Hub, stop the client and
        disconnect from the Hub.
        """
        cliEx = None
        try:
            self.client.unregister()
        except SAMPClientError as cliEx:
            pass

        if self.client.isRunning():
            self.client.stop()
        self.hub.disconnect()

        if cliEx:
            raise cliEx

    # HUB
    def ping(self):
        """
        Proxy to `ping` SAMP Hub method (Standard Profile only).
        """
        return self.hub.ping()

    def declareMetadata(self, metadata):
        """
        Proxy to `declareMetadata` SAMP Hub method.
        """
        return self.client.declareMetadata(metadata)

    def getMetadata(self, client_id):
        """
        Proxy to `getMetadata` SAMP Hub method.
        """
        return self.hub.getMetadata(self.client.getPrivateKey(), client_id)

    def getSubscriptions(self, client_id):
        """
        Proxy to `getSubscriptions` SAMP Hub method.
        """
        return self.hub.getSubscriptions(self.client.getPrivateKey(), client_id)

    def getRegisteredClients(self):
        """
        Proxy to `getRegisteredClients` SAMP Hub method.
        """
        return self.hub.getRegisteredClients(self.client.getPrivateKey())

    def getSubscribedClients(self, mtype):
        """
        Proxy to `getSubscribedClients` SAMP Hub method.
        """
        return self.hub.getSubscribedClients(self.client.getPrivateKey(), mtype)

    def _format_easy_msg(self, mtype, params):

        msg = {}

        if "extra_kws" in params:
            extra = params["extra_kws"]
            del(params["extra_kws"])
            msg = {"samp.mtype": mtype, "samp.params": params}
            msg.update(extra)
        else:
            msg = {"samp.mtype": mtype, "samp.params": params}

        return msg

    def notify(self, recipient_id, message):
        """Proxy to `notify` SAMP Hub method."""
        return self.hub.notify(self.client.getPrivateKey(), recipient_id, message)

    def enotify(self, recipient_id, mtype, **params):
        """
        Easy to use `notify`.

        This is a proxy to `notify` method that allows to
        send the notification message in a simplified way.

        Note that reserved `extra_kws` keyword is a dictionary with the special meaning of
        being used to add extra keywords, in addition to the standard `samp.mtype`
        and `samp.params`, to the message sent.

        Parameters
        ----------
        recipient_id : str
            Recipient ID

        mtype : str
            the MType to be notified

        params : dict or set of keywords
            Variable keyword set which contains the list of parameters for
            the specified MType.

        Examples
        --------
        >>> import astropy.vo.samp as sampy
        >>> cli = sampy.SAMPIntegratedClient()
        >>> ...
        >>> cli.enotify("samp.msg.progress", msgid = "xyz", txt = "initialization",
        ...             percent = "10", extra_kws = {"my.extra.info": "just an example"})
        """
        return self.notify(recipient_id, self._format_easy_msg(mtype, params))

    def notifyAll(self, message):
        """
        Proxy to `notifyAll` SAMP Hub method.
        """
        return self.hub.notifyAll(self.client.getPrivateKey(), message)

    def enotifyAll(self, mtype, **params):
        """
        Easy to use `notify`.

        This is a proxy to `notifyAll` method that allows to
        send the notification message in a simplified way.


        Note that reserved `extra_kws` keyword is a dictionary with the special meaning of
        being used to add extra keywords, in addition to the standard `samp.mtype`
        and `samp.params`, to the message sent.

        Parameters
        ----------
        mtype : str
            MType to be notified.

        params : dict or set of keywords
            Variable keyword set which contains the list of parameters for
            the specified MType.

        Examples
        --------
        >>> import astropy.vo.samp as sampy
        >>> cli = sampy.SAMPIntegratedClient()
        >>> ...
        >>> cli.enotifyAll("samp.msg.progress", txt = "initialization",
        ...                percent = "10", extra_kws = {"my.extra.info": "just an example"})
        """
        return self.notifyAll(self._format_easy_msg(mtype, params))

    def call(self, recipient_id, msg_tag, message):
        """
        Proxy to `call` SAMP Hub method."""
        return self.hub.call(self.client.getPrivateKey(), recipient_id, msg_tag, message)

    def ecall(self, recipient_id, msg_tag, mtype, **params):
        """
        Easy to use `call`.

        This is a proxy to `call` method that allows to
        send a call message in a simplified way.

        Note that reserved `extra_kws` keyword is a dictionary with the special meaning of
        being used to add extra keywords, in addition to the standard `samp.mtype`
        and `samp.params`, to the message sent.

        Parameters
        ----------
        recipient_id : str
            Recipient ID

        msg_tag : str
            Message tag to use

        mtype : str
            MType to be sent

        params : dict of set of keywords
            Variable keyword set which contains the list of parameters for
            the specified MType.

        Examples
        --------
        >>> import astropy.vo.samp as sampy
        >>> cli = sampy.SAMPIntegratedClient()
        >>> ...
        >>> msgid = cli.ecall("abc", "xyz", "samp.msg.progress", txt = "initialization",
        ...                   percent = "10", extra_kws = {"my.extra.info": "just an example"})
        """

        return self.call(recipient_id, msg_tag, self._format_easy_msg(mtype, params))

    def callAll(self, msg_tag, message):
        """Proxy to `callAll` SAMP Hub method."""
        return self.hub.callAll(self.client.getPrivateKey(), msg_tag, message)

    def ecallAll(self, msg_tag, mtype, **params):
        """
        Easy to use `callAll`.

        This is a proxy to `callAll` method that allows to
        send the call message in a simplified way.

        Note that reserved `extra_kws` keyword is a dictionary with the special meaning of
        being used to add extra keywords, in addition to the standard `samp.mtype`
        and `samp.params`, to the message sent.

        Parameters
        ----------
        msg_tag : str
            Message tag to use

        mtype : str
            MType to be sent

        params : dict of set of keywords
            Variable keyword set which contains the list of parameters for
            the specified MType.

        Examples
        --------
        >>> import astropy.vo.samp as sampy
        >>> cli = sampy.SAMPIntegratedClient()
        >>> ...
        >>> msgid = cli.ecallAll("xyz", "samp.msg.progress", txt = "initialization",
        ...                      percent = "10", extra_kws = {"my.extra.info": "just an example"})
        """
        self.callAll(msg_tag, self._format_easy_msg(mtype, params))

    def callAndWait(self, recipient_id, message, timeout):
        """
        Proxy to `callAndWait` SAMP Hub method.

        If timeout expires a `SAMPProxyError` instance is raised.
        """
        return self.hub.callAndWait(self.client.getPrivateKey(), recipient_id, message, timeout)

    def ecallAndWait(self, recipient_id, mtype, timeout, **params):
        """
        Easy to use `callAndWait`.

        This is a proxy to `callAll` method that allows to
        send the call message in a simplified way.

        Note that reserved `extra_kws` keyword is a dictionary with the special meaning of
        being used to add extra keywords, in addition to the standard `samp.mtype`
        and `samp.params`, to the message sent.

        recipient_id : str
            Recipient ID

        mtype : str
            MType to be sent

        timeout : str
            Call timeout in seconds

        params : dict of set of keywords
            Variable keyword set which contains the list of parameters for
            the specified MType.

        Examples
        --------
        >>> import astropy.vo.samp as sampy
        >>> cli = sampy.SAMPIntegratedClient()
        >>> ...
        >>> cli.ecallAndWait("xyz", "samp.msg.progress", "5", txt = "initialization",
        ...                  percent = "10", extra_kws = {"my.extra.info": "just an example"})
        """
        return self.callAndWait(recipient_id, self._format_easy_msg(mtype, params), timeout)

    def reply(self, msg_id, response):
        """Proxy to `reply` SAMP Hub method."""
        return self.hub.reply(self.client.getPrivateKey(), msg_id, response)

    def _format_easy_response(self, status, result, error):

        msg = {"samp.status": status}
        if result != None:
            msg.update({"samp.result": result})
        if error != None:
            msg.update({"samp.error": error})

        return msg

    def ereply(self, msg_id, status, result=None, error=None):
        """
        Easy to use `reply`.

        This is a proxy to `callAll` method that allows to
        send a reply message in a simplified way.

        Parameters
        ----------
        @msg_id : str
            Message ID to which reply.

        status : str
            Content of the `samp.status` response keyword.

        result : dict
            Content of the `samp.result` response keyword.

        error : dict
            Content of the `samp.error` response keyword.

        Examples
        --------
        >>> import astropy.vo.samp as sampy
        >>> cli = sampy.SAMPIntegratedClient()
        >>> ...
        >>> cli.ereply("abd", sampy.SAMP_STATUS_ERROR, result={},
        ...            error={"samp.errortxt": "Test error message"})
        """
        return self.reply(msg_id, self._format_easy_response(status, result, error))

    # CLIENT
    def receiveNotification(self, private_key, sender_id, message):
        """
        Standard callable client `receiveNotification` method.

        This method is automatically handled when the `bindReceiveNotification`
        method is used to bind distinct operations to MTypes.
        In case of a customized callable client implementation that inherits from
        the `SAMPClient` class this method should be overwritten.

        ATTENTION: When overwritten, this method must always return
                   a string result (even empty).

        Parameters
        ----------
        private_key : str
            Client private key.

        sender_id : str
            Sender public ID.

        message : dict
            Received message.

        Returns
        -------
        confirmation : str
            Any confirmation string.
        """
        return self.client.receiveNotification(private_key, sender_id, message)

    def receiveCall(self, private_key, sender_id, msg_id, message):
        """
        Standard callable client C{receiveCall} method.

        This method is automatically handled when the `bindReceiveCall` method
        is used to bind distinct operations to MTypes.
        In case of a customized callable client implementation that inherits from
        the `SAMPClient` class this method should be overwritten.

        ATTENTION: When overwritten, this method must always return
                   a string result (even empty).

        Parameters
        ----------
        private_key : str
            Client private key.

        sender_id : str
            Sender public ID.

        msg_id : str
            Received message ID.

        message : dict
            Received message.

        Returns
        -------
        confirmation : str
            Any confirmation string.
        """
        return self.client.receiveCall(private_key, sender_id, msg_id, message)

    def receiveResponse(self, private_key, responder_id, msg_tag, response):
        """
        Standard callable client `receiveResponse` method.

        This method is automatically handled when the `bindReceiveResponse` method
        is used to bind distinct operations to MTypes.
        In case of a customized callable client implementation that inherits from
        the `SAMPClient` class this method should be overwritten.

        ATTENTION: When overwritten, this method must always return
                   a string result (even empty).

        Parameters
        ----------
        private_key : str
            Client private key.

        responder_id : str
            Responder public ID.

        msg_tag : str
            Response message tag.

        response : dict
            Received response.

        Returns
        -------
        confirmation : str
            Any confirmation string.
        """
        return self.client.receiveResponse(private_key, responder_id, msg_tag, response)

    def bindReceiveMessage(self, mtype, function, declare=True, metadata=None):
        """
        Bind a specific MType to a function or class method, being intended for
        a call or a notification.

        The function must be of the form:
        `def my_function_or_method(<self,> private_key, sender_id, msg_id, mtype, params, extra)`

        where `private_key` is the client private-key, the `sender_id` argument is the
        notification sender ID, `msg_id` is the Hub message-id (calls only, otherwise is None),
        `mtype` is the message MType, `params is the message parameter set (content of
        "samp.params") and `extra` is a dictionary containing any extra message map entry.
        The client is automatically declared subscribed to the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be catched.

        function : callable
            Application function to be used when C{mtype} is received.

        declare : bool
            Specify whether the client must be automatically declared as
            subscribed to the MType (see also `declareSubscriptions`).

        metadata : dict
            An optional map containing additional metadata to declare associated
            with the MType subscribed to (see also `declareSubscriptions`).
        """
        self.client.bindReceiveMessage(mtype, function, declare=True, metadata=None)

    def bindReceiveNotification(self, mtype, function, declare=True, metadata=None):
        """
        Bind a specific MType notification to a function or class method.

        The function must be of the form:

        `def my_function_or_method(<self,> private_key, sender_id, mtype, params, extra)`

        where `private_key` is the client private-key, `sender_id` argument is
        the notification sender ID, `mtype` is the message MType, `params` is
        the notified message parameter set (content of "samp.params") and `extra` is a
        dictionary containing any extra message map entry. The client is
        automatically declared subscribed to the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be catched.

        function : callable
            Application function to be used when C{mtype} is received.

        declare : bool
            Specify whether the client must be automatically declared as
            subscribed to the MType (see alse `declareSubscriptions`).

        metadata : dict, optional
            Dictionary containing additional metadata to declare associated
            with the MType subscribed to (see also `declareSubscriptions`).
        """
        self.client.bindReceiveNotification(mtype, function, declare, metadata)

    def bindReceiveCall(self, mtype, function, declare=True, metadata=None):
        """
        Bind a specific MType call to a function or class method.

        The function must be of the form:

        `def my_function_or_method(<self,> private_key, sender_id, msg_id, mtype, params, extra)`

        where `private_key` is the client private-key, `sender_id` argument is the
        notification sender ID, `msg_id` is the Hub message-id, `mtype` is the message MType,
        `params` is the message parameter set (content of "samp.params") and `extra` is a
        dictionary containing any extra message map entry. The client is
        automatically declared subscribed to the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be catched.

        function : callable
            Application function to be used when `mtype` is received.

        declare : bool
            Specify whether the client must be automatically declared as
            subscribed to the MType (see also L`declareSubscriptions`).

        metadata : dict, optional
            Dictionary containing additional metadata to declare associated
            with the MType subscribed to (see also `declareSubscriptions`).
        """
        self.client.bindReceiveCall(mtype, function, declare, metadata)

    def bindReceiveResponse(self, msg_tag, function):
        """
        Bind a specific msg-tag response to a function or class method.

        The function must be of the form:

        `def my_function_or_method(<self,> private_key, responder_id, msg_tag, response)`

        where `private_key` is the client private-key, `responder_id` argument is the message
        responder ID, `msg_tag` is the message-tag provided at call time and `response` is the
        response received.

        Parameters
        ----------
        msg_tag : str
            Message-tag to be catched.

        function : callable
            Application function to be used when `msg_tag` is received.
        """
        self.client.bindReceiveResponse(msg_tag, function)

    def unbindReceiveNotification(self, mtype, declare=True):
        """
        Remove from the notifications binding table the specified MType and unsubscribe
        the client from it (if required).

        Parameters
        ----------
        mtype : str
            MType to be removed.

        declare : bool
            Specify whether the client must be automatically declared as
            unsubscribed from the MType (see alse `declareSubscriptions`).
        """
        self.client.unbindReceiveNotification(mtype, declare)

    def unbindReceiveCall(self, mtype, declare=True):
        """
        Remove from the calls binding table the specified MType and unsubscribe
        the client from it (if required).

        Parameters
        ----------
        mtype : str
            MType to be removed.

        declare : bool
            Specify whether the client must be automatically declared as
            unsubscribed from the MType (see alse `declareSubscriptions`).
        """
        self.client.unbindReceiveCall(mtype, declare)

    def unbindReceiveResponse(self, msg_tag):
        """Remove from the responses binding table the specified message-tag.

        Parameters
        ----------
        msg_tag : str
            Message-tag to be removed.
        """
        self.client.unbindReceiveResponse(msg_tag)

    def declareSubscriptions(self, subscriptions=None):
        """
        Declares the MTypes the client wishes to subscribe to,
        implicitly defined with the MType binding methods `bindReceiveNotification` and `bindReceiveCall`.

        An optional `subscriptions` map can be added to the final map passed to
        the `SAMPHubProxy.declareSubscriptions` operation.

        Parameters
        ----------
        subscriptions: dict, optional
            Dictionary containing the list of MTypes to subscribe to,
            with the same format of the `subscriptions` map passed to the
            `SAMPHubProxy.declareSubscriptions` operation.
        """
        self.client.declareSubscriptions(subscriptions)

    def getPrivateKey(self):
        """
        Return the client private key used for the Standard Profile communications
        obtained at registration time (`samp.private-key`).

        Returns
        -------
        key : str
            The client private key
        """
        return self.client.getPrivateKey()

    def getPublicId(self):
        """
        Return public client ID obtained at registration time (`samp.self-id`).

        Returns
        -------
        id : str
            Client public ID
        """
        return self.client.getPublicId()


class SAMPClient(object):
    """
    Utility class which provides facilities to create and manage a SAMP compliant
    XML-RPC server that acts as SAMP callable client application.

    Parameters
    ----------
    hub : `SAMPHubProxy`
        An instance of `SAMPHubProxy` to be used for messaging with the SAMP Hub.

    name : str, optional
        Client name (corresponding to `samp.name` metadata keyword).

    description : str, optional
        Client description (corresponding to `samp.description.text` metadata keyword).

    metadata : dict, optional
        Client application metadata in the standard SAMP format.

    addr : str
        Listening address (or IP).

    port : int, optional
        Listening XML-RPC server socket port

    https : bool
        Set the callable client running on a Secure Sockets Layer connection (HTTPS)?
        By default SSL is disabled.

    key_file : str
        Set the file containing the private key for SSL connections. If the
        certificate file (`cert_file`) contains the private key, then `key_file` can be omitted.

    cert_file : str
        Specify the file which contains a certificate to be used to identify the
        local side of the secure connection.

    cert_reqs : int
        The parameter `cert_reqs` specifies whether a certificate is required
        from the Hub side of the connection, and whether it will be validated if provided. It
        must be one of the three values `ssl.CERT_NONE` (certificates ignored), `ssl.CERT_OPTIONAL`
        (not required, but validated if provided), or `ssl.CERT_REQUIRED` (required and validated).
        If the value of this parameter is not `ssl.CERT_NONE`, then the `ca_certs` parameter must
        point to a file of CA certificates.

    ca_certs : str
        The `ca_certs` file contains a set of concatenated "Certification Authority"
        certificates, which are used to validate the certificate passed from the Hub end of the
        connection.

    ssl_version : int
        The `ssl_version` option specifies which version of the SSL protocol to use.
        Typically, the server chooses a particular protocol version, and the client must adapt to the
        server's choice. Most of the versions are not interoperable with the other versions. If not
        specified the default SSL version is `ssl.PROTOCOL_SSLv23`. This version provides the most
        compatibility with other versions Hub side. Other SSL protocol versions are:
        `ssl.PROTOCOL_SSLv2`, `ssl.PROTOCOL_SSLv3` and `ssl.PROTOCOL_TLSv1`.

    callable : bool
        Is the client callable?
    """

    def __init__(self, hub, name=None, description=None, metadata=None,
                 addr=None, port=0, https=False, key_file=None, cert_file=None,
                 cert_reqs=0, ca_certs=None, ssl_version=2, callable=True):
        # GENERAL
        self._thread = None
        self._is_running = False

        if metadata == None:
            metadata = {}

        if name != None:
            metadata["samp.name"] = name

        if description != None:
            metadata["samp.description.text"] = description

        self._metadata = metadata

        self._addr = addr
        self._port = port
        self._xmlrpcAddr = None
        self._callable = callable

        # HUB INTERACTION
        self.client = None
        self._public_id = None
        self._private_key = None
        self._hub_id = None
        self._notification_bindings = {}
        self._call_bindings = {"samp.app.ping": [self._ping, {}],
                               "client.env.get": [self._client_env_get, {}]}
        self._response_bindings = {}

        self._host_name = "127.0.0.1"
        try:
            if internet_on():
                self._host_name = socket.getfqdn()
        except:
            pass

        self.hub = hub

        if self._callable:

            if SSL_SUPPORT and https:
                self.client = SecureXMLRPCServer((self._addr or self._host_name, self._port),
                                                 keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                                                 self._log, logRequests=False, allow_none=True)
            else:
                self.client = ThreadingXMLRPCServer((self._addr or self._host_name,
                                                     self._port), logRequests=False, allow_none=True)

            self.client.register_introspection_functions()
            self.client.register_function(self.receiveNotification, 'samp.client.receiveNotification')
            self.client.register_function(self.receiveCall, 'samp.client.receiveCall')
            self.client.register_function(self.receiveResponse, 'samp.client.receiveResponse')

            if self._port == 0:
                self._port = self.client.socket.getsockname()[1]

            self._xmlrpcAddr = "http://%s:%s" % (self._addr or
                                                 self._host_name,
                                                 self._port)

    def __del__(self):
        self.stop()

    def start(self):
        """Start the client in a non-blocking way."""
        global _THREAD_STARTED_COUNT
        _THREAD_STARTED_COUNT += 1
        self._is_running = True
        self._run_client()

    def stop(self, timeout=0.1):
        """
        Stop the client.

        Parameters
        ----------
        timeout : float
            Timeout after wich the client terminates even if the threading is still alive.
        """
        self._is_running = False
        if self._thread is not None:
            self._thread.join(timeout)
            self._thread = None

    def isRunning(self):
        """
        Return an information concerning the client running status.

        Returns
        -------
        running : bool
            `True` if the client is running, `False` otherwise.
        """
        return self._is_running != None

    def _run_client(self):
        if self._callable:
            self._thread = threading.Thread(target=self._serve_forever)
            self._thread.setDaemon(True)
            self._thread.start()

    def _serve_forever(self):
        while self._is_running:
            try:
                r = w = e = None
                r, w, e = select.select([self.client.socket], [], [], 0.1)
            except:
                pass
            if r:
                self.client.handle_request()

    def _ping(self, private_key, sender_id, msg_id, msg_mtype, msg_params, message):
        self.hub.reply(private_key, msg_id, {"samp.status": SAMP_STATUS_OK, "samp.result": {}})

    def _client_env_get(self, private_key, sender_id, msg_id, msg_mtype, msg_params, message):
        if msg_params["name"] in os.environ:
            self.hub.reply(private_key, msg_id, {"samp.status": SAMP_STATUS_OK,
                                                 "samp.result": {"value": os.environ[msg_params["name"]]}})
        else:
            self.hub.reply(private_key, msg_id, {"samp.status": SAMP_STATUS_WARNING,
                                                 "samp.result": {"value": ""},
                                                 "samp.error": {"samp.errortxt":
                                                                "Environment variable not defined."}})

    def _handle_notification(self, private_key, sender_id, message):

        if private_key == self.getPrivateKey() and "samp.mtype" in message:

            msg_mtype = message["samp.mtype"]
            del message["samp.mtype"]
            msg_params = message["samp.params"]
            del message["samp.params"]

            msubs = SAMPHubServer.getMTypeSubtypes(msg_mtype)
            for mtype in msubs:
                if mtype in self._notification_bindings:
                    bound_func = self._notification_bindings[mtype][0]
                    if ((inspect.ismethod(bound_func) and bound_func.__func__.__code__.co_argcount == 6) or
                        (inspect.isfunction(bound_func) and bound_func.__code__.co_argcount == 5)):
                        bound_func(private_key, sender_id, msg_mtype, msg_params, message)
                    else:
                        bound_func(private_key, sender_id, None, msg_mtype, msg_params, message)

        return ""

    def receiveNotification(self, private_key, sender_id, message):
        """
        Standard callable client `receiveNotification` method.

        This method is automatically handled when the `bindReceiveNotification`
        method is used to bind distinct operations to MTypes.
        In case of a customized callable client implementation that inherits from
        the `SAMPClient` class this method should be overwritten.

        ATTENTION: When overwritten, this method must always return
                   a string result (even empty).

        Parameters
        ----------
        private_key : str
            Client private key.

        sender_id : str
            Sender public ID.

        message : dict
            Received message.

        Returns
        -------
        confirmation : str
            Any confirmation string.
        """
        return self._handle_notification(private_key, sender_id, message)

    def _handle_call(self, private_key, sender_id, msg_id, message):
        if private_key == self.getPrivateKey() and "samp.mtype" in message:

            msg_mtype = message["samp.mtype"]
            del message["samp.mtype"]
            msg_params = message["samp.params"]
            del message["samp.params"]

            msubs = SAMPHubServer.getMTypeSubtypes(msg_mtype)

            for mtype in msubs:
                if mtype in self._call_bindings:
                    self._call_bindings[mtype][0](private_key, sender_id, msg_id,
                                                  msg_mtype, msg_params, message)

        return ""

    def receiveCall(self, private_key, sender_id, msg_id, message):
        """
        Standard callable client C{receiveCall} method.

        This method is automatically handled when the `bindReceiveCall` method
        is used to bind distinct operations to MTypes.
        In case of a customized callable client implementation that inherits from
        the `SAMPClient` class this method should be overwritten.

        ATTENTION: When overwritten, this method must always return
                   a string result (even empty).

        Parameters
        ----------
        private_key : str
            Client private key.

        sender_id : str
            Sender public ID.

        msg_id : str
            Message ID received.

        message : dict
            Received message.

        Returns
        -------
        confimation : str
            Any confirmation string.
        """
        return self._handle_call(private_key, sender_id, msg_id, message)

    def _handle_response(self, private_key, responder_id, msg_tag, response):
        if private_key == self.getPrivateKey() and msg_tag in self._response_bindings:
            self._response_bindings[msg_tag](private_key, responder_id, msg_tag, response)
        return ""

    def receiveResponse(self, private_key, responder_id, msg_tag, response):
        """
        Standard callable client `receiveResponse` method.

        This method is automatically handled when the `bindReceiveResponse` method
        is used to bind distinct operations to MTypes.
        In case of a customized callable client implementation that inherits from
        the `SAMPClient` class this method should be overwritten.

        ATTENTION: When overwritten, this method must always return
                   a string result (even empty).

        Parameters
        ----------
        private_key : str
            Client private key.

        responder_id : str
            Responder public ID.

        msg_tag : str
            Response message tag.

        response : dict
            Received response.

        Returns
        -------
        confirmation : str
            Any confirmation string.
        """
        return self._handle_response(private_key, responder_id, msg_tag, response)

    def bindReceiveMessage(self, mtype, function, declare=True, metadata=None):
        """
        Bind a specific MType to a function or class method, being intended for
        a call or a notification.

        The function must be of the form:
        `def my_function_or_method(<self,> private_key, sender_id, msg_id, mtype, params, extra)`

        where `private_key` is the client private-key, the `sender_id` argument is the
        notification sender ID, `msg_id` is the Hub message-id (calls only, otherwise is `None`),
        `mtype` is the message MType, `params` is the message parameter set (content of
        "samp.params") and `extra` is a dictionary containing any extra message map entry.
        The client is automatically declared subscribed to the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be catched.

        function : callable
            Application function to be used when `mtype` is received.

        declare : bool
            Specify whether the client must be automatically declared as
            subscribed to the MType (see also `declareSubscriptions`).

        metadata : dict, optional
            Dictionary containing additional metadata to declare associated
            with the MType subscribed to (see also `declareSubscriptions`).
        """
        self.bindReceiveCall(mtype, function, declare=True, metadata=None)
        self.bindReceiveNotification(mtype, function, declare=True, metadata=None)

    def bindReceiveNotification(self, mtype, function, declare=True, metadata=None):
        """
        Bind a specific MType notification to a function or class method.

        The function must be of the form:

        `def my_function_or_method(<self,> private_key, sender_id, mtype, params, extra)`

        where `private_key` is the client private-key, `sender_id` argument is
        the notification sender ID, `mtype` is the message MType, `params` is
        the notified message parameter set (content of "samp.params") and `extra` is a
        dictionary containing any extra message map entry. The client is
        automatically declared subscribed to the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be catched.

        function : callable
            Application function to be used when `mtype` is received.

        declare : bool
            Specify whether the client must be automatically declared as
            subscribed to the MType (see alse `declareSubscriptions`).

        metadata : dict, optional
            Dictionary containing additional metadata to declare associated
            with the MType subscribed to (see also `declareSubscriptions`).
        """
        if self._callable:
            if not metadata:
                metadata = {}
            self._notification_bindings[mtype] = [function, metadata]
            if declare:
                self._declareSubscriptions()
        else:
            raise SAMPClientError("Client not callable.")

    def bindReceiveCall(self, mtype, function, declare=True, metadata=None):
        """
        Bind a specific MType call to a function or class method.

        The function must be of the form:

        `def my_function_or_method(<self,> private_key, sender_id, msg_id, mtype, params, extra)`

        where `private_key` is the client private-key, `sender_id` argument is the
        notification sender ID, `msg_id` is the Hub message-id, `mtype` is the message MType,
        `params` is the message parameter set (content of "samp.params") and `extra` is a
        dictionary containing any extra message map entry. The client is
        automatically declared subscribed to the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be catched.

        function : callable
            Application function to be used when C{mtype} is received.

        declare : bool
            Specify whether the client must be automatically declared as
            subscribed to the MType (see also `declareSubscriptions`).

        metadata : dict, optional
            Dictionary containing additional metadata to declare associated
            with the MType subscribed to (see also `declareSubscriptions`).
        """
        if self._callable:
            if not metadata:
                metadata = {}
            self._call_bindings[mtype] = [function, metadata]
            if declare:
                self._declareSubscriptions()
        else:
            raise SAMPClientError("Client not callable.")

    def bindReceiveResponse(self, msg_tag, function):
        """
        Bind a specific msg-tag response to a function or class method.

        The function must be of the form:

        `def my_function_or_method(<self,> private_key, responder_id, msg_tag, response)`

        where `private_key` is the client private-key, `responder_id` argument is the message
        responder ID, `msg_tag` is the message-tag provided at call time and `response` is the
        response received.

        Parameters
        ----------
        msg_tag : str
            Message-tag to be catched.

        function : callable
            Application function to be used when `msg_tag` is received.
        """
        if self._callable:
            self._response_bindings[msg_tag] = function
        else:
            raise SAMPClientError("Client not callable.")

    def unbindReceiveNotification(self, mtype, declare=True):
        """
        Remove from the notifications binding table the specified MType and unsubscribe
        the client from it (if required).

        Parameters
        ----------
        mtype : str
            MType to be removed.

        declare : bool
            Specify whether the client must be automatically declared as
            unsubscribed from the MType (see alse `declareSubscriptions`).
        """
        if self._callable:
            del self._notification_bindings[mtype]
            if declare:
                self._declareSubscriptions()
        else:
            raise SAMPClientError("Client not callable.")

    def unbindReceiveCall(self, mtype, declare=True):
        """
        Remove from the calls binding table the specified MType and unsubscribe
        the client from it (if required).

        Parameters
        ----------
        mtype : str
            MType to be removed.

        declare : bool
            Specify whether the client must be automatically declared as
            unsubscribed from the MType (see alse `declareSubscriptions`).
        """
        if self._callable:
            del self._call_bindings[mtype]
            if declare:
                self._declareSubscriptions()
        else:
            raise SAMPClientError("Client not callable.")

    def unbindReceiveResponse(self, msg_tag):
        """
        Remove from the responses binding table the specified message-tag.

        Parameters
        ----------
        msg_tag : str
            Message-tag to be removed.
        """
        if self._callable:
            del self._response_bindings[msg_tag]
        else:
            raise SAMPClientError("Client not callable.")

    def declareSubscriptions(self, subscriptions=None):
        """
        Declares the MTypes the client wishes to subscribe to, implicitly defined
        with the MType binding methods `bindReceiveNotification` and `bindReceiveCall`.

        An optional `subscriptions` map can be added to the final map passed to
        the `SAMPHubProxy.declareSubscriptions` operation.

        Parameters
        ----------
        subscriptions : dict, optional
            Dictionary containing the list of MTypes to subscribe to,
            with the same format of the C{subscriptions} map passed to the
            `SAMPHubProxy.declareSubscriptions` operation.
        """
        if self._callable:
            self._declareSubscriptions(subscriptions)
        else:
            raise SAMPClientError("Client not callable.")

    def register(self):
        """
        Register the client to the SAMP Hub.

        If the registration fails a `SAMPClientError` is raised.
        """
        if self.hub.isConnected():

            if self._private_key != None:
                raise SAMPClientError("Client already registered")

            try:
                result = self.hub.register(self.hub.lockfile["samp.secret"])
                if result["samp.self-id"] == "" or result["samp.private-key"] == "":
                    raise SAMPClientError("Registation failed. Probably the secret code is wrong.")
                self._public_id = result["samp.self-id"]
                self._private_key = result["samp.private-key"]
                self._hub_id = result["samp.hub-id"]
                if self._callable:
                    self._setXmlrpcCallback()
                    self._declareSubscriptions()
                if self._metadata != {}:
                    self.declareMetadata()
            except SAMPProxyError as err:
                raise SAMPClientError(err.faultString)
            except:
                raise SAMPClientError("Unexpected error: registration failed")

        else:
            raise SAMPClientError("Unable to register to the SAMP Hub. Hub proxy not connected.")

    def unregister(self):
        """
        Unregister the client from the SAMP Hub.

        If the unregistration fails a `SAMPClientError` is raised.
        """
        if self.hub.isConnected():

            try:
                self.hub.unregister(self._private_key)
                self._hub_id = None
                self._public_id = None
                self._private_key = None
            except:
                raise SAMPClientError("Unable to unregister from the SAMP Hub.")
        else:
            raise SAMPClientError("Unable to unregister from the SAMP Hub. Hub proxy not connected.")

    def _setXmlrpcCallback(self):
        if self.hub.isConnected() and self._private_key != None:

            try:
                self.hub.setXmlrpcCallback(self._private_key,
                                           self._xmlrpcAddr)
            except:
                pass

    def _declareSubscriptions(self, subscriptions=None):
        if self.hub.isConnected() and self._private_key != None:

            try:
                mtypes_dict = {}
                # Collect notification mtypes and metadata
                for mtype in self._notification_bindings.keys():
                    mtypes_dict[mtype] = copy.deepcopy(self._notification_bindings[mtype][1])

                # Collect notification mtypes and metadata
                for mtype in self._call_bindings.keys():
                    mtypes_dict[mtype] = copy.deepcopy(self._call_bindings[mtype][1])

                # Add optional subscription map
                if subscriptions:
                    mtypes_dict.update(copy.deepcopy(subscriptions))

                self.hub.declareSubscriptions(self._private_key, mtypes_dict)

            except Exception as ex:
                raise SAMPClientError("Unable to declare subscriptions. Hub unreachable or not connected or client not registered (%s)." %str(ex))
        else:
            raise SAMPClientError("Unable to declare subscriptions. Hub unreachable or not connected or client not registered.")

    def declareMetadata(self, metadata=None):
        """
        Declare the client application metadata supported.

        Parameters
        ----------
        metadata : dict, optional
            Dictionary containig the client application metadata
            as defined in the SAMP definition document. If omitted, then none metadata are
            declared.
        """
        if self.hub.isConnected() and self._private_key != None:

            try:
                if metadata != None:
                    self._metadata.update(metadata)

                self.hub.declareMetadata(self._private_key, self._metadata)
            except:
                raise SAMPClientError("Unable to declare metadata. Hub unreachable or not connected or client not registered.")
        else:
            raise SAMPClientError("Unable to declare metadata. Hub unreachable or not connected or client not registered.")

    def getPrivateKey(self):
        """
        Return the client private key used for the Standard Profile communications
        obtained at registration time (`samp.private-key`).

        Returns
        -------
        key : str
            Client private key.
        """
        return self._private_key

    def getPublicId(self):
        """
        Return public client ID obtained at registration time (`samp.self-id`).

        Returns
        -------
        id : str
            Client public ID.
        """
        return self._public_id

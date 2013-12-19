# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .client import SAMPClient
from .errors import SAMPClientError
from .hub_proxy import SAMPHubProxy

__all__ = ['SAMPIntegratedClient']

__doctest_skip__ = ['SAMPIntegratedClient.*']


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
        Proxy to `call` SAMP Hub method.
        """
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
        return self.client.receiveNotification(private_key, sender_id, message)

    receiveNotification.__doc__ = SAMPClient.receiveNotification.__doc__

    def receiveCall(self, private_key, sender_id, msg_id, message):
        return self.client.receiveCall(private_key, sender_id, msg_id, message)

    receiveCall.__doc__ = SAMPClient.receiveCall.__doc__

    def receiveResponse(self, private_key, responder_id, msg_tag, response):
        return self.client.receiveResponse(private_key, responder_id, msg_tag, response)

    receiveResponse.__doc__ = SAMPClient.receiveResponse.__doc__

    def bindReceiveMessage(self, mtype, function, declare=True, metadata=None):
        self.client.bindReceiveMessage(mtype, function, declare=True, metadata=None)

    bindReceiveMessage.__doc__ = SAMPClient.bindReceiveMessage.__doc__

    def bindReceiveNotification(self, mtype, function, declare=True, metadata=None):
        self.client.bindReceiveNotification(mtype, function, declare, metadata)

    bindReceiveNotification.__doc__ = SAMPClient.bindReceiveNotification.__doc__

    def bindReceiveCall(self, mtype, function, declare=True, metadata=None):
        self.client.bindReceiveCall(mtype, function, declare, metadata)

    bindReceiveCall.__doc__ = SAMPClient.bindReceiveCall.__doc__

    def bindReceiveResponse(self, msg_tag, function):
        self.client.bindReceiveResponse(msg_tag, function)

    bindReceiveResponse.__doc__ = SAMPClient.bindReceiveResponse.__doc__

    def unbindReceiveNotification(self, mtype, declare=True):
        self.client.unbindReceiveNotification(mtype, declare)

    unbindReceiveNotification.__doc__ = SAMPClient.unbindReceiveNotification.__doc__

    def unbindReceiveCall(self, mtype, declare=True):
        self.client.unbindReceiveCall(mtype, declare)

    unbindReceiveCall.__doc__ = SAMPClient.unbindReceiveCall.__doc__

    def unbindReceiveResponse(self, msg_tag):
        self.client.unbindReceiveResponse(msg_tag)

    unbindReceiveResponse.__doc__ = SAMPClient.unbindReceiveResponse.__doc__

    def declareSubscriptions(self, subscriptions=None):
        self.client.declareSubscriptions(subscriptions)

    declareSubscriptions.__doc__ = SAMPClient.declareSubscriptions.__doc__

    def getPrivateKey(self):
        return self.client.getPrivateKey()

    getPrivateKey.__doc__ = SAMPClient.getPrivateKey.__doc__

    def getPublicId(self):
        return self.client.getPublicId()

    getPublicId.__doc__ = SAMPClient.getPublicId.__doc__

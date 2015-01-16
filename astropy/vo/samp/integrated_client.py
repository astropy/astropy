# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .client import SAMPClient
from .hub_proxy import SAMPHubProxy

__all__ = ['SAMPIntegratedClient']

__doctest_skip__ = ['SAMPIntegratedClient.*']


class SAMPIntegratedClient(object):
    """
    A Simple SAMP client.

    This class is meant to simplify the client usage providing a proxy class
    that merges the :class:`~astropy.vo.samp.SAMPClient` and
    :class:`~astropy.vo.samp.SAMPHubProxy` functionalities in a
    simplified API.

    Parameters
    ----------
    name : str, optional
        Client name (corresponding to ``samp.name`` metadata keyword).

    description : str, optional
        Client description (corresponding to ``samp.description.text`` metadata
        keyword).

    metadata : dict, optional
        Client application metadata in the standard SAMP format.

    addr : str, optional
        Listening address (or IP). This defaults to 127.0.0.1 if the internet
        is not reachable, otherwise it defaults to the host name.

    port : int, optional
        Listening XML-RPC server socket port. If left set to 0 (the default),
        the operating system will select a free port.

    https : bool, optional
        If `True`, set the callable client running on a Secure Sockets Layer
        (SSL) connection (HTTPS). By default SSL is disabled.

    key_file : str, optional
        The path to a file containing the private key for SSL connections. If
        the certificate file (``cert_file``) contains the private key, then
        ``key_file`` can be omitted.

    cert_file : str, optional
        The path to a file which contains a certificate to be used to identify
        the local side of the secure connection.

    cert_reqs : int, optional
        Whether a certificate is required from the server side of the
        connection, and whether it will be validated if provided. It must be
        one of the three values `ssl.CERT_NONE` (certificates ignored),
        `ssl.CERT_OPTIONAL` (not required, but validated if provided), or
        `ssl.CERT_REQUIRED` (required and validated). If the value of this
        parameter is not `ssl.CERT_NONE`, then the ``ca_certs`` parameter must
        point to a file of CA certificates.

    ca_certs : str, optional
        The path to a file containing a set of concatenated "Certification
        Authority" certificates, which are used to validate the certificate
        passed from the Hub end of the connection.

    ssl_version : int, optional
        Which version of the SSL protocol to use. Typically, the
        server chooses a particular protocol version, and the client
        must adapt to the server's choice. Most of the versions are
        not interoperable with the other versions. If not specified,
        the default SSL version is taken from the default in the
        installed version of the Python standard `ssl` library.  See
        the `ssl` documentation for more information.

    callable : bool, optional
        Whether the client can receive calls and notifications. If set to
        `False`, then the client can send notifications and calls, but can not
        receive any.
    """

    def __init__(self, name=None, description=None, metadata=None,
                 addr=None, port=0, https=False, key_file=None, cert_file=None,
                 cert_reqs=0, ca_certs=None, ssl_version=None, callable=True):

        self.hub = SAMPHubProxy()

        self.client_arguments = {
            'name': name,
            'description': description,
            'metadata': metadata,
            'addr': addr,
            'port': port,
            'https': https,
            'key_file': key_file,
            'cert_file': cert_file,
            'cert_reqs': cert_reqs,
            'ca_certs': ca_certs,
            'ssl_version': ssl_version,
            'callable': callable,
        }
        """
        Collected arguments that should be passed on to the SAMPClient below.
        The SAMPClient used to be instantiated in __init__; however, this
        caused problems with disconnecting and reconnecting to the HUB.
        The client_arguments is used to maintain backwards compatibility.
        """

        self.client = None
        "The client will be instantiated upon connect()."

    # GENERAL

    @property
    def is_connected(self):
        """
        Testing method to verify the client connection with a running Hub.

        Returns
        -------
        is_connected : bool
            True if the client is connected to a Hub, False otherwise.
        """
        return self.hub.is_connected and self.client.is_running

    def connect(self, hub=None, hub_params=None,
                key_file=None, cert_file=None, cert_reqs=0,
                ca_certs=None, ssl_version=None, pool_size=20):
        """
        Connect with the current or specified SAMP Hub, start and register the
        client.

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
        self.hub.connect(hub, hub_params, key_file, cert_file,
                         cert_reqs, ca_certs, ssl_version, pool_size)

        # The client has to be instantiated here and not in __init__() because
        # this allows disconnecting and reconnecting to the HUB. Nonetheless,
        # the client_arguments are set in __init__() because the
        # instantiation of the client used to happen there and this retains
        # backwards compatibility.
        self.client = SAMPClient(
            self.hub,
            **self.client_arguments
        )
        self.client.start()
        self.client.register()

    def disconnect(self):
        """
        Unregister the client from the current SAMP Hub, stop the client and
        disconnect from the Hub.
        """
        if self.is_connected:
            try:
                self.client.unregister()
            finally:
                if self.client.is_running:
                    self.client.stop()
                self.hub.disconnect()

    # HUB
    def ping(self):
        """
        Proxy to ``ping`` SAMP Hub method (Standard Profile only).
        """
        return self.hub.ping()

    def declare_metadata(self, metadata):
        """
        Proxy to ``declareMetadata`` SAMP Hub method.
        """
        return self.client.declare_metadata(metadata)

    def get_metadata(self, client_id):
        """
        Proxy to ``getMetadata`` SAMP Hub method.
        """
        return self.hub.get_metadata(self.get_private_key(), client_id)

    def get_subscriptions(self, client_id):
        """
        Proxy to ``getSubscriptions`` SAMP Hub method.
        """
        return self.hub.get_subscriptions(self.get_private_key(), client_id)

    def get_registered_clients(self):
        """
        Proxy to ``getRegisteredClients`` SAMP Hub method.

        This returns all the registered clients, excluding the current client.
        """
        return self.hub.get_registered_clients(self.get_private_key())

    def get_subscribed_clients(self, mtype):
        """
        Proxy to ``getSubscribedClients`` SAMP Hub method.
        """
        return self.hub.get_subscribed_clients(self.get_private_key(), mtype)

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
        """
        Proxy to ``notify`` SAMP Hub method.
        """
        return self.hub.notify(self.get_private_key(), recipient_id, message)

    def enotify(self, recipient_id, mtype, **params):
        """
        Easy to use version of :meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.notify`.

        This is a proxy to ``notify`` method that allows to send the
        notification message in a simplified way.

        Note that reserved ``extra_kws`` keyword is a dictionary with the
        special meaning of being used to add extra keywords, in addition to
        the standard ``samp.mtype`` and ``samp.params``, to the message sent.

        Parameters
        ----------
        recipient_id : str
            Recipient ID

        mtype : str
            the MType to be notified

        params : dict or set of keywords
            Variable keyword set which contains the list of parameters for the
            specified MType.

        Examples
        --------
        >>> from astropy.vo.samp import SAMPIntegratedClient
        >>> cli = SAMPIntegratedClient()
        >>> ...
        >>> cli.enotify("samp.msg.progress", msgid = "xyz", txt = "initialization",
        ...             percent = "10", extra_kws = {"my.extra.info": "just an example"})
        """
        return self.notify(recipient_id, self._format_easy_msg(mtype, params))

    def notify_all(self, message):
        """
        Proxy to ``notifyAll`` SAMP Hub method.
        """
        return self.hub.notify_all(self.get_private_key(), message)

    def enotify_all(self, mtype, **params):
        """
        Easy to use version of :meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.notify_all`.

        This is a proxy to ``notifyAll`` method that allows to send the
        notification message in a simplified way.

        Note that reserved ``extra_kws`` keyword is a dictionary with the
        special meaning of being used to add extra keywords, in addition to
        the standard ``samp.mtype`` and ``samp.params``, to the message sent.

        Parameters
        ----------
        mtype : str
            MType to be notified.

        params : dict or set of keywords
            Variable keyword set which contains the list of parameters for
            the specified MType.

        Examples
        --------
        >>> from astropy.vo.samp import SAMPIntegratedClient
        >>> cli = SAMPIntegratedClient()
        >>> ...
        >>> cli.enotify_all("samp.msg.progress", txt = "initialization",
        ...                 percent = "10",
        ...                 extra_kws = {"my.extra.info": "just an example"})
        """
        return self.notify_all(self._format_easy_msg(mtype, params))

    def call(self, recipient_id, msg_tag, message):
        """
        Proxy to ``call`` SAMP Hub method.
        """
        return self.hub.call(self.get_private_key(), recipient_id, msg_tag, message)

    def ecall(self, recipient_id, msg_tag, mtype, **params):
        """
        Easy to use version of :meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.call`.

        This is a proxy to ``call`` method that allows to send a call message
        in a simplified way.

        Note that reserved ``extra_kws`` keyword is a dictionary with the
        special meaning of being used to add extra keywords, in addition to
        the standard ``samp.mtype`` and ``samp.params``, to the message sent.

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
        >>> from astropy.vo.samp import SAMPIntegratedClient
        >>> cli = SAMPIntegratedClient()
        >>> ...
        >>> msgid = cli.ecall("abc", "xyz", "samp.msg.progress",
        ...                   txt = "initialization", percent = "10",
        ...                   extra_kws = {"my.extra.info": "just an example"})
        """

        return self.call(recipient_id, msg_tag, self._format_easy_msg(mtype, params))

    def call_all(self, msg_tag, message):
        """
        Proxy to ``callAll`` SAMP Hub method.
        """
        return self.hub.call_all(self.get_private_key(), msg_tag, message)

    def ecall_all(self, msg_tag, mtype, **params):
        """
        Easy to use version of :meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.call_all`.

        This is a proxy to ``callAll`` method that allows to send the call
        message in a simplified way.

        Note that reserved ``extra_kws`` keyword is a dictionary with the
        special meaning of being used to add extra keywords, in addition to
        the standard ``samp.mtype`` and ``samp.params``, to the message sent.

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
        >>> from astropy.vo.samp import SAMPIntegratedClient
        >>> cli = SAMPIntegratedClient()
        >>> ...
        >>> msgid = cli.ecall_all("xyz", "samp.msg.progress",
        ...                       txt = "initialization", percent = "10",
        ...                       extra_kws = {"my.extra.info": "just an example"})
        """
        self.call_all(msg_tag, self._format_easy_msg(mtype, params))

    def call_and_wait(self, recipient_id, message, timeout):
        """
        Proxy to ``callAndWait`` SAMP Hub method.
        """
        return self.hub.call_and_wait(self.get_private_key(), recipient_id, message, timeout)

    def ecall_and_wait(self, recipient_id, mtype, timeout, **params):
        """
        Easy to use version of :meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.call_and_wait`.

        This is a proxy to ``callAndWait`` method that allows to send the call
        message in a simplified way.

        Note that reserved ``extra_kws`` keyword is a dictionary with the
        special meaning of being used to add extra keywords, in addition to
        the standard ``samp.mtype`` and ``samp.params``, to the message sent.

        Parameters
        ----------
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
        >>> from astropy.vo.samp import SAMPIntegratedClient
        >>> cli = SAMPIntegratedClient()
        >>> ...
        >>> cli.ecall_and_wait("xyz", "samp.msg.progress", "5",
        ...                    txt = "initialization", percent = "10",
        ...                    extra_kws = {"my.extra.info": "just an example"})
        """
        return self.call_and_wait(recipient_id, self._format_easy_msg(mtype, params), timeout)

    def reply(self, msg_id, response):
        """
        Proxy to ``reply`` SAMP Hub method.
        """
        return self.hub.reply(self.get_private_key(), msg_id, response)

    def _format_easy_response(self, status, result, error):

        msg = {"samp.status": status}
        if result is not None:
            msg.update({"samp.result": result})
        if error is not None:
            msg.update({"samp.error": error})

        return msg

    def ereply(self, msg_id, status, result=None, error=None):
        """
        Easy to use version of :meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.reply`.

        This is a proxy to ``reply`` method that allows to send a reply
        message in a simplified way.

        Parameters
        ----------
        msg_id : str
            Message ID to which reply.

        status : str
            Content of the ``samp.status`` response keyword.

        result : dict
            Content of the ``samp.result`` response keyword.

        error : dict
            Content of the ``samp.error`` response keyword.

        Examples
        --------
        >>> from astropy.vo.samp import SAMPIntegratedClient, SAMP_STATUS_ERROR
        >>> cli = SAMPIntegratedClient()
        >>> ...
        >>> cli.ereply("abd", SAMP_STATUS_ERROR, result={},
        ...            error={"samp.errortxt": "Test error message"})
        """
        return self.reply(msg_id, self._format_easy_response(status, result, error))

    # CLIENT

    def receive_notification(self, private_key, sender_id, message):
        return self.client.receive_notification(private_key, sender_id, message)

    receive_notification.__doc__ = SAMPClient.receive_notification.__doc__

    def receive_call(self, private_key, sender_id, msg_id, message):
        return self.client.receive_call(private_key, sender_id, msg_id, message)

    receive_call.__doc__ = SAMPClient.receive_call.__doc__

    def receive_response(self, private_key, responder_id, msg_tag, response):
        return self.client.receive_response(private_key, responder_id, msg_tag, response)

    receive_response.__doc__ = SAMPClient.receive_response.__doc__

    def bind_receive_message(self, mtype, function, declare=True, metadata=None):
        self.client.bind_receive_message(mtype, function, declare=True, metadata=None)

    bind_receive_message.__doc__ = SAMPClient.bind_receive_message.__doc__

    def bind_receive_notification(self, mtype, function, declare=True, metadata=None):
        self.client.bind_receive_notification(mtype, function, declare, metadata)

    bind_receive_notification.__doc__ = SAMPClient.bind_receive_notification.__doc__

    def bind_receive_call(self, mtype, function, declare=True, metadata=None):
        self.client.bind_receive_call(mtype, function, declare, metadata)

    bind_receive_call.__doc__ = SAMPClient.bind_receive_call.__doc__

    def bind_receive_response(self, msg_tag, function):
        self.client.bind_receive_response(msg_tag, function)

    bind_receive_response.__doc__ = SAMPClient.bind_receive_response.__doc__

    def unbind_receive_notification(self, mtype, declare=True):
        self.client.unbind_receive_notification(mtype, declare)

    unbind_receive_notification.__doc__ = SAMPClient.unbind_receive_notification.__doc__

    def unbind_receive_call(self, mtype, declare=True):
        self.client.unbind_receive_call(mtype, declare)

    unbind_receive_call.__doc__ = SAMPClient.unbind_receive_call.__doc__

    def unbind_receive_response(self, msg_tag):
        self.client.unbind_receive_response(msg_tag)

    unbind_receive_response.__doc__ = SAMPClient.unbind_receive_response.__doc__

    def declare_subscriptions(self, subscriptions=None):
        self.client.declare_subscriptions(subscriptions)

    declare_subscriptions.__doc__ = SAMPClient.declare_subscriptions.__doc__

    def get_private_key(self):
        return self.client.get_private_key()

    get_private_key.__doc__ = SAMPClient.get_private_key.__doc__

    def get_public_id(self):
        return self.client.get_public_id()

    get_public_id.__doc__ = SAMPClient.get_public_id.__doc__

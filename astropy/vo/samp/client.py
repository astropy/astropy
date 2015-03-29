# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy
import os
import select
import socket
import threading
import warnings

from ... import log
from ...extern.six.moves.urllib.parse import urlunparse

from .constants import SAMP_STATUS_OK, SAMP_STATUS_WARNING
from .hub import SAMPHubServer
from .errors import SAMPClientError, SAMPWarning
from .utils import internet_on, get_num_args

from .constants import SSL_SUPPORT
from .standard_profile import ThreadingXMLRPCServer

if SSL_SUPPORT:
    import ssl
    from .ssl_utils import SecureXMLRPCServer

__all__ = ['SAMPClient']


class SAMPClient(object):
    """
    Utility class which provides facilities to create and manage a SAMP
    compliant XML-RPC server that acts as SAMP callable client application.

    Parameters
    ----------
    hub : :class:`~astropy.vo.samp.SAMPHubProxy`
        An instance of :class:`~astropy.vo.samp.SAMPHubProxy` to be
        used for messaging with the SAMP Hub.

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

    # TODO: define what is meant by callable

    def __init__(self, hub, name=None, description=None, metadata=None,
                 addr=None, port=0, https=False, key_file=None, cert_file=None,
                 cert_reqs=0, ca_certs=None, ssl_version=None, callable=True):

        # GENERAL
        self._is_running = False
        self._is_registered = False

        if metadata is None:
            metadata = {}

        if name is not None:
            metadata["samp.name"] = name

        if description is not None:
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
        if internet_on():
            try:
                self._host_name = socket.getfqdn()
                socket.getaddrinfo(self._addr or self._host_name, self._port or 0)
            except socket.error:
                self._host_name = "127.0.0.1"

        self.hub = hub

        if self._callable:

            self._thread = threading.Thread(target=self._serve_forever)
            self._thread.daemon = True

            if SSL_SUPPORT and https:
                self.client = SecureXMLRPCServer((self._addr or self._host_name, self._port),
                                                 key_file, cert_file, cert_reqs, ca_certs, ssl_version,
                                                 log, logRequests=False, allow_none=True)
            else:
                self.client = ThreadingXMLRPCServer((self._addr or self._host_name,
                                                     self._port), logRequests=False, allow_none=True)

            self.client.register_introspection_functions()
            self.client.register_function(self.receive_notification, 'samp.client.receiveNotification')
            self.client.register_function(self.receive_call, 'samp.client.receiveCall')
            self.client.register_function(self.receive_response, 'samp.client.receiveResponse')

            # If the port was set to zero, then the operating system has
            # selected a free port. We now check what this port number is.
            if self._port == 0:
                self._port = self.client.socket.getsockname()[1]

            if SSL_SUPPORT and https:
                protocol = 'https'
            else:
                protocol = 'http'

            self._xmlrpcAddr = urlunparse((protocol,
                                           '{0}:{1}'.format(self._addr or self._host_name,
                                                            self._port),
                                           '', '', '', ''))

    def start(self):
        """
        Start the client in a separate thread (non-blocking).

        This only has an effect if ``callable`` was set to `True` when
        initializing the client.
        """
        if self._callable:
            self._is_running = True
            self._run_client()

    def stop(self, timeout=10.):
        """
        Stop the client.

        Parameters
        ----------
        timeout : float
            Timeout after which to give up if the client cannot be cleanly
            shut down.
        """
        # Setting _is_running to False causes the loop in _serve_forever to
        # exit. The thread should then stop running. We wait for the thread to
        # terminate until the timeout, then we continue anyway.
        self._is_running = False
        if self._callable and self._thread.is_alive():
            self._thread.join(timeout)
        if self._thread.is_alive():
            raise SAMPClientError("Client was not shut down successfully "
                                  "(timeout={0}s)".format(timeout))

    @property
    def is_running(self):
        """
        Whether the client is currently running.
        """
        return self._is_running

    @property
    def is_registered(self):
        """
        Whether the client is currently registered.
        """
        return self._is_registered

    def _run_client(self):
        if self._callable:
            self._thread.start()

    def _serve_forever(self):
        while self._is_running:
            try:
                read_ready = select.select([self.client.socket], [], [], 0.1)[0]
            except select.error as exc:
                warnings.warn("Call to select in SAMPClient failed: {0}".format(exc),
                              SAMPWarning)
            else:
                if read_ready:
                    self.client.handle_request()

    def _ping(self, private_key, sender_id, msg_id, msg_mtype, msg_params,
              message):

        reply = {"samp.status": SAMP_STATUS_OK, "samp.result": {}}

        self.hub.reply(private_key, msg_id, reply)

    def _client_env_get(self, private_key, sender_id, msg_id, msg_mtype,
                        msg_params, message):

        if msg_params["name"] in os.environ:
            reply = {"samp.status": SAMP_STATUS_OK,
                     "samp.result": {"value": os.environ[msg_params["name"]]}}
        else:
            reply = {"samp.status": SAMP_STATUS_WARNING,
                     "samp.result": {"value": ""},
                     "samp.error": {"samp.errortxt":
                                    "Environment variable not defined."}}

        self.hub.reply(private_key, msg_id, reply)

    def _handle_notification(self, private_key, sender_id, message):

        if private_key == self.get_private_key() and "samp.mtype" in message:

            msg_mtype = message["samp.mtype"]
            del message["samp.mtype"]
            msg_params = message["samp.params"]
            del message["samp.params"]

            msubs = SAMPHubServer.get_mtype_subtypes(msg_mtype)
            for mtype in msubs:
                if mtype in self._notification_bindings:
                    bound_func = self._notification_bindings[mtype][0]
                    if get_num_args(bound_func) == 5:
                        bound_func(private_key, sender_id, msg_mtype,
                                   msg_params, message)
                    else:
                        bound_func(private_key, sender_id, None, msg_mtype,
                                   msg_params, message)

        return ""

    def receive_notification(self, private_key, sender_id, message):
        """
        Standard callable client ``receive_notification`` method.

        This method is automatically handled when the
        :meth:`~astropy.vo.samp.client.SAMPClient.bind_receive_notification`
        method is used to bind distinct operations to MTypes. In case of a
        customized callable client implementation that inherits from the
        :class:`~astropy.vo.samp.SAMPClient` class this method should be
        overwritten.

        .. note:: When overwritten, this method must always return
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

        if private_key == self.get_private_key() and "samp.mtype" in message:

            msg_mtype = message["samp.mtype"]
            del message["samp.mtype"]
            msg_params = message["samp.params"]
            del message["samp.params"]

            msubs = SAMPHubServer.get_mtype_subtypes(msg_mtype)

            for mtype in msubs:
                if mtype in self._call_bindings:
                    self._call_bindings[mtype][0](private_key, sender_id,
                                                  msg_id, msg_mtype,
                                                  msg_params, message)

        return ""

    def receive_call(self, private_key, sender_id, msg_id, message):
        """
        Standard callable client ``receive_call`` method.

        This method is automatically handled when the
        :meth:`~astropy.vo.samp.client.SAMPClient.bind_receive_call` method is
        used to bind distinct operations to MTypes. In case of a customized
        callable client implementation that inherits from the
        :class:`~astropy.vo.samp.SAMPClient` class this method should be
        overwritten.

        .. note:: When overwritten, this method must always return
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
        confirmation : str
            Any confirmation string.
        """
        return self._handle_call(private_key, sender_id, msg_id, message)

    def _handle_response(self, private_key, responder_id, msg_tag, response):
        if (private_key == self.get_private_key() and
            msg_tag in self._response_bindings):
            self._response_bindings[msg_tag](private_key, responder_id,
                                    msg_tag, response)
        return ""

    def receive_response(self, private_key, responder_id, msg_tag, response):
        """
        Standard callable client ``receive_response`` method.

        This method is automatically handled when the
        :meth:`~astropy.vo.samp.client.SAMPClient.bind_receive_response` method
        is used to bind distinct operations to MTypes. In case of a customized
        callable client implementation that inherits from the
        :class:`~astropy.vo.samp.SAMPClient` class this method should be
        overwritten.

        .. note:: When overwritten, this method must always return
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
        return self._handle_response(private_key, responder_id, msg_tag,
                                     response)

    def bind_receive_message(self, mtype, function, declare=True,
                             metadata=None):
        """
        Bind a specific MType to a function or class method, being intended for
        a call or a notification.

        The function must be of the form::

            def my_function_or_method(<self,> private_key, sender_id, msg_id,
                                      mtype, params, extra)

        where ``private_key`` is the client private-key, ``sender_id`` is the
        notification sender ID, ``msg_id`` is the Hub message-id (calls only,
        otherwise is `None`), ``mtype`` is the message MType, ``params`` is the
        message parameter set (content of ``"samp.params"``) and ``extra`` is a
        dictionary containing any extra message map entry. The client is
        automatically declared subscribed to the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be catched.

        function : callable
            Application function to be used when ``mtype`` is received.

        declare : bool, optional
            Specify whether the client must be automatically declared as
            subscribed to the MType (see also
            :meth:`~astropy.vo.samp.client.SAMPClient.declare_subscriptions`).

        metadata : dict, optional
            Dictionary containing additional metadata to declare associated
            with the MType subscribed to (see also
            :meth:`~astropy.vo.samp.client.SAMPClient.declare_subscriptions`).
        """

        self.bind_receive_call(mtype, function, declare=declare,
                               metadata=metadata)

        self.bind_receive_notification(mtype, function, declare=declare,
                                       metadata=metadata)

    def bind_receive_notification(self, mtype, function, declare=True, metadata=None):
        """
        Bind a specific MType notification to a function or class method.

        The function must be of the form::

            def my_function_or_method(<self,> private_key, sender_id, mtype,
                                      params, extra)

        where ``private_key`` is the client private-key, ``sender_id`` is the
        notification sender ID, ``mtype`` is the message MType, ``params`` is
        the notified message parameter set (content of ``"samp.params"``) and
        ``extra`` is a dictionary containing any extra message map entry. The
        client is automatically declared subscribed to the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be caught.

        function : callable
            Application function to be used when ``mtype`` is received.

        declare : bool, optional
            Specify whether the client must be automatically declared as
            subscribed to the MType (see also
            :meth:`~astropy.vo.samp.client.SAMPClient.declare_subscriptions`).

        metadata : dict, optional
            Dictionary containing additional metadata to declare associated
            with the MType subscribed to (see also
            :meth:`~astropy.vo.samp.client.SAMPClient.declare_subscriptions`).
        """
        if self._callable:
            if not metadata:
                metadata = {}
            self._notification_bindings[mtype] = [function, metadata]
            if declare:
                self._declare_subscriptions()
        else:
            raise SAMPClientError("Client not callable.")

    def bind_receive_call(self, mtype, function, declare=True, metadata=None):
        """
        Bind a specific MType call to a function or class method.

        The function must be of the form::

            def my_function_or_method(<self,> private_key, sender_id, msg_id,
                                      mtype, params, extra)

        where ``private_key`` is the client private-key, ``sender_id`` is the
        notification sender ID, ``msg_id`` is the Hub message-id, ``mtype`` is
        the message MType, ``params`` is the message parameter set (content of
        ``"samp.params"``) and ``extra`` is a dictionary containing any extra
        message map entry. The client is automatically declared subscribed to
        the MType by default.

        Parameters
        ----------
        mtype : str
            MType to be caught.

        function : callable
            Application function to be used when ``mtype`` is received.

        declare : bool, optional
            Specify whether the client must be automatically declared as
            subscribed to the MType (see also
            :meth:`~astropy.vo.samp.client.SAMPClient.declare_subscriptions`).

        metadata : dict, optional
            Dictionary containing additional metadata to declare associated
            with the MType subscribed to (see also
            :meth:`~astropy.vo.samp.client.SAMPClient.declare_subscriptions`).
        """
        if self._callable:
            if not metadata:
                metadata = {}
            self._call_bindings[mtype] = [function, metadata]
            if declare:
                self._declare_subscriptions()
        else:
            raise SAMPClientError("Client not callable.")

    def bind_receive_response(self, msg_tag, function):
        """
        Bind a specific msg-tag response to a function or class method.

        The function must be of the form::

            def my_function_or_method(<self,> private_key, responder_id,
                                      msg_tag, response)

        where ``private_key`` is the client private-key, ``responder_id`` is
        the message responder ID, ``msg_tag`` is the message-tag provided at
        call time and ``response`` is the response received.

        Parameters
        ----------
        msg_tag : str
            Message-tag to be caught.

        function : callable
            Application function to be used when ``msg_tag`` is received.
        """
        if self._callable:
            self._response_bindings[msg_tag] = function
        else:
            raise SAMPClientError("Client not callable.")

    def unbind_receive_notification(self, mtype, declare=True):
        """
        Remove from the notifications binding table the specified MType and
        unsubscribe the client from it (if required).

        Parameters
        ----------
        mtype : str
            MType to be removed.

        declare : bool
            Specify whether the client must be automatically declared as
            unsubscribed from the MType (see also
            :meth:`~astropy.vo.samp.client.SAMPClient.declare_subscriptions`).
        """
        if self._callable:
            del self._notification_bindings[mtype]
            if declare:
                self._declare_subscriptions()
        else:
            raise SAMPClientError("Client not callable.")

    def unbind_receive_call(self, mtype, declare=True):
        """
        Remove from the calls binding table the specified MType and unsubscribe
        the client from it (if required).

        Parameters
        ----------
        mtype : str
            MType to be removed.

        declare : bool
            Specify whether the client must be automatically declared as
            unsubscribed from the MType (see also
            :meth:`~astropy.vo.samp.client.SAMPClient.declare_subscriptions`).
        """
        if self._callable:
            del self._call_bindings[mtype]
            if declare:
                self._declare_subscriptions()
        else:
            raise SAMPClientError("Client not callable.")

    def unbind_receive_response(self, msg_tag):
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

    def declare_subscriptions(self, subscriptions=None):
        """
        Declares the MTypes the client wishes to subscribe to, implicitly
        defined with the MType binding methods
        :meth:`~astropy.vo.samp.client.SAMPClient.bind_receive_notification`
        and :meth:`~astropy.vo.samp.client.SAMPClient.bind_receive_call`.

        An optional ``subscriptions`` map can be added to the final map passed
        to the :meth:`~astropy.vo.samp.hub_proxy.SAMPHubProxy.declare_subscriptions`
        method.

        Parameters
        ----------
        subscriptions : dict, optional
            Dictionary containing the list of MTypes to subscribe to, with the
            same format of the ``subscriptions`` map passed to the
            :meth:`~astropy.vo.samp.hub_proxy.SAMPHubProxy.declare_subscriptions`
            method.
        """
        if self._callable:
            self._declare_subscriptions(subscriptions)
        else:
            raise SAMPClientError("Client not callable.")

    def register(self):
        """
        Register the client to the SAMP Hub.
        """
        if self.hub.is_connected:

            if self._private_key is not None:
                raise SAMPClientError("Client already registered")

            result = self.hub.register(self.hub.lockfile["samp.secret"])

            if result["samp.self-id"] == "":
                raise SAMPClientError("Registration failed - "
                                      "samp.self-id was not set by the hub.")

            if result["samp.private-key"] == "":
                raise SAMPClientError("Registration failed - "
                                      "samp.private-key was not set by the hub.")

            self._public_id = result["samp.self-id"]
            self._private_key = result["samp.private-key"]
            self._hub_id = result["samp.hub-id"]

            if self._callable:
                self._set_xmlrpc_callback()
                self._declare_subscriptions()

            if self._metadata != {}:
                self.declare_metadata()

            self._is_registered = True

        else:
            raise SAMPClientError("Unable to register to the SAMP Hub. "
                                  "Hub proxy not connected.")

    def unregister(self):
        """
        Unregister the client from the SAMP Hub.
        """
        if self.hub.is_connected:
            self._is_registered = False
            self.hub.unregister(self._private_key)
            self._hub_id = None
            self._public_id = None
            self._private_key = None
        else:
            raise SAMPClientError("Unable to unregister from the SAMP Hub. "
                                  "Hub proxy not connected.")

    def _set_xmlrpc_callback(self):
        if self.hub.is_connected and self._private_key is not None:
            self.hub.set_xmlrpc_callback(self._private_key,
                                         self._xmlrpcAddr)

    def _declare_subscriptions(self, subscriptions=None):
        if self.hub.is_connected and self._private_key is not None:

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

            self.hub.declare_subscriptions(self._private_key, mtypes_dict)

        else:
            raise SAMPClientError("Unable to declare subscriptions. Hub "
                                  "unreachable or not connected or client "
                                  "not registered.")

    def declare_metadata(self, metadata=None):
        """
        Declare the client application metadata supported.

        Parameters
        ----------
        metadata : dict, optional
            Dictionary containing the client application metadata as defined in
            the SAMP definition document. If omitted, then no metadata are
            declared.
        """
        if self.hub.is_connected and self._private_key is not None:
            if metadata is not None:
                self._metadata.update(metadata)
            self.hub.declare_metadata(self._private_key, self._metadata)
        else:
            raise SAMPClientError("Unable to declare metadata. Hub "
                                  "unreachable or not connected or client "
                                  "not registered.")

    def get_private_key(self):
        """
        Return the client private key used for the Standard Profile
        communications obtained at registration time (``samp.private-key``).

        Returns
        -------
        key : str
            Client private key.
        """
        return self._private_key

    def get_public_id(self):
        """
        Return public client ID obtained at registration time
        (``samp.self-id``).

        Returns
        -------
        id : str
            Client public ID.
        """
        return self._public_id

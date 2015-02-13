# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy
import os
import select
import socket
import threading
import time
import uuid
import warnings

from ...extern.six.moves import queue
from ...extern.six.moves import xmlrpc_client as xmlrpc
from ...extern.six.moves.urllib.parse import urlunparse
from ... import log

from .constants import SAMP_STATUS_OK
from .constants import __profile_version__
from .errors import SAMPWarning, SAMPHubError, SAMPProxyError
from .utils import internet_on, ServerProxyPool, _HubAsClient
from .lockfile_helpers import read_lockfile, create_lock_file

from .standard_profile import ThreadingXMLRPCServer
from .web_profile import WebProfileXMLRPCServer, web_profile_text_dialog

from .constants import SSL_SUPPORT
from ...utils.exceptions import AstropyWarning

if SSL_SUPPORT:
    import ssl
    from .ssl_utils import SafeTransport, SecureXMLRPCServer

__all__ = ['SAMPHubServer', 'WebProfileDialog']

__doctest_skip__ = ['.', 'SAMPHubServer.*']


class SAMPHubServer(object):
    """
    SAMP Hub Server.

    The SSL parameters are usable only if Python has been compiled with SSL
    support and/or `ssl <http://docs.python.org/dev/library/ssl.html>`_ module
    is installed (available by default since Python 2.6).

    Parameters
    ----------
    secret : str, optional
        The secret code to use for the SAMP lockfile. If none is is specified,
        the :func:`uuid.uuid1` function is used to generate one.

    addr : str, optional
        Listening address (or IP). This defaults to 127.0.0.1 if the internet
        is not reachable, otherwise it defaults to the host name.

    port : int, optional
        Listening XML-RPC server socket port. If left set to 0 (the default),
        the operating system will select a free port.

    lockfile : str, optional
        Custom lockfile name.

    timeout : int, optional
        Hub inactivity timeout. If ``timeout > 0`` then the Hub automatically
        stops after an inactivity period longer than ``timeout`` seconds. By
        default ``timeout`` is set to 0 (Hub never expires).

    client_timeout : int, optional
        Client inactivity timeout. If ``client_timeout > 0`` then the Hub
        automatically unregisters the clients which result inactive for a
        period longer than ``client_timeout`` seconds. By default
        ``client_timeout`` is set to 0 (clients never expire).

    mode : str, optional
        Defines the Hub running mode. If ``mode`` is ``'single'`` then the Hub
        runs using the standard ``.samp`` lock-file, having a single instance
        for user desktop session. Otherwise, if ``mode`` is ``'multiple'``,
        then the Hub runs using a non-standard lock-file, placed in
        ``.samp-1`` directory, of the form ``samp-hub-<UUID>``, where
        ``<UUID>`` is a unique UUID assigned to the hub.

    label : str, optional
        A string used to label the Hub with a human readable name. This string
        is written in the lock-file assigned to the ``hub.label`` token.

    https : bool, optional
        Set the Hub running on a Secure Sockets Layer connection (HTTPS). By
        default SSL is disabled.

    key_file : str, optional
        The path to a file containing the private key for SSL connections. If
        the certificate file (``cert_file``) contains the private key, then
        ``key_file`` can be omitted.

    cert_file : str, optional
        The path to a file containing a certificate to be used to identify the
        local side of the secure connection.

    cert_reqs : int, optional
        The parameter ``cert_reqs`` specifies whether a certificate is
        required from the client side of the connection, and whether it will
        be validated if provided. It must be one of the three values
        `ssl.CERT_NONE` (certificates ignored), `ssl.CERT_OPTIONAL` (not
        required, but validated if provided), or `ssl.CERT_REQUIRED` (required
        and validated). If the value of this parameter is not `ssl.CERT_NONE`,
        then the ``ca_certs`` parameter must point to a file of CA
        certificates.

    ca_certs : str, optional
        The path to a file containing a set of concatenated "Certification
        Authority" certificates, which are used to validate the certificate
        passed from the Hub end of the connection.

    ssl_version : int, optional
        The ``ssl_version`` option specifies which version of the SSL
        protocol to use. Typically, the server chooses a particular
        protocol version, and the client must adapt to the server's
        choice. Most of the versions are not interoperable with the
        other versions. If not specified, the default SSL version is
        taken from the default in the installed version of the Python
        standard `ssl` library.  See the `ssl` documentation for more
        information.

    web_profile : bool, optional
        Enables or disables the Web Profile support.

    web_profile_dialog : class, optional
        Allows a class instance to be specified using ``web_profile_dialog``
        to replace the terminal-based message with e.g. a GUI pop-up. Two
        `queue.Queue` instances will be added to the instance as attributes
        ``queue_request`` and ``queue_result``. When a request is received via
        the ``queue_request`` queue, the pop-up should be displayed, and a
        value of `True` or `False` should be added to ``queue_result``
        depending on whether the user accepted or refused the connection.

    web_port : int, optional
        The port to use for web SAMP. This should not be changed except for
        testing purposes, since web SAMP should always use port 21012.

    pool_size : int, optional
        The number of socket connections opened to communicate with the
        clients.
    """

    def __init__(self, secret=None, addr=None, port=0, lockfile=None,
                 timeout=0, client_timeout=0, mode='single', label="",
                 https=False, key_file=None, cert_file=None, cert_reqs=0,
                 ca_certs=None, ssl_version=None, web_profile=True,
                 web_profile_dialog=None, web_port=21012, pool_size=20):

        # Generate random ID for the hub
        self._id = str(uuid.uuid1())

        # General settings
        self._is_running = False
        self._customlockfilename = lockfile
        self._lockfile = None
        self._addr = addr
        self._port = port
        self._mode = mode
        self._label = label
        self._timeout = timeout
        self._client_timeout = client_timeout
        self._pool_size = pool_size

        self._web_profile = web_profile
        self._web_profile_server = None
        self._web_profile_callbacks = {}
        self._web_profile_requests_queue = queue.Queue(1)
        self._web_profile_requests_result = queue.Queue(1)
        self._web_profile_requests_semaphore = queue.Queue(1)

        self._web_profile_dialog = web_profile_dialog
        if self._web_profile_dialog is not None:
            # TODO: Some sort of duck-typing on the web_profile_dialog object
            self._web_profile_dialog.queue_request = self._web_profile_requests_queue
            self._web_profile_dialog.queue_result = self._web_profile_requests_result
        self._web_port = web_port

        if web_profile:
            try:
                self._web_profile_server = WebProfileXMLRPCServer(('localhost', self._web_port), log,
                                                                  logRequests=False,
                                                                  allow_none=True)
                self._web_port = self._web_profile_server.socket.getsockname()[1]
                self._web_profile_server.register_introspection_functions()
                log.info("Hub set to run with Web Profile support enabled.")
            except socket.error:
                warnings.warn("Port {0} already in use. Impossible to run the "
                              "Hub with Web Profile support.".format(self._web_port),
                              SAMPWarning)
                self._web_profile = web_profile = False

        # SSL general settings
        self._https = https
        self._key_file = key_file
        self._cert_file = cert_file
        self._cert_reqs = cert_reqs
        self._ca_certs = ca_certs
        self._ssl_version = ssl_version

        self._host_name = "127.0.0.1"
        if internet_on():
            try:
                self._host_name = socket.getfqdn()
                socket.getaddrinfo(self._addr or self._host_name,
                                   self._port or 0)
            except socket.error:
                self._host_name = "127.0.0.1"

        # XML-RPC server settings
        if https:

            if key_file is not None and not os.path.isfile(key_file):
                raise SAMPHubError("Unable to load SSL private key file!")

            if cert_file is None or not os.path.isfile(cert_file):
                raise SAMPHubError("Unable to load SSL cert file!")

            log.info("Hub set for using SSL.")
            self._server = SecureXMLRPCServer((self._addr or self._host_name,
                                               self._port or 0),
                                              key_file, cert_file, cert_reqs,
                                              ca_certs, ssl_version, log,
                                              logRequests=False,
                                              allow_none=True)

            self._port = self._server.socket.getsockname()[1]
            self._url = urlunparse(('https',
                                    "{0}:{1}".format(self._addr or self._host_name,
                                                     self._port),
                                    '', '', '', ''))
        else:

            self._server = ThreadingXMLRPCServer((self._addr or self._host_name,
                                                  self._port or 0),
                                                 log, logRequests=False,
                                                 allow_none=True)

            self._port = self._server.socket.getsockname()[1]
            self._url = urlunparse(('http',
                                    "{0}:{1}".format(self._addr or self._host_name,
                                                     self._port),
                                    '', '', '', ''))

        self._server.register_introspection_functions()

        # Threading stuff
        self._thread_lock = threading.Lock()
        self._thread_run = threading.Thread(target=self._serve_forever)
        self._thread_run.daemon = True

        if self._timeout > 0:
            self._thread_hub_timeout = threading.Thread(target=self._timeout_test_hub,
                                                        name="Hub timeout test")
            self._thread_hub_timeout.daemon = True
        else:
            self._thread_hub_timeout = None

        if self._client_timeout > 0:
            self._thread_client_timeout = threading.Thread(target=self._timeout_test_client,
                                                               name="Client timeout test")
            self._thread_client_timeout.daemon = True
        else:
            self._thread_client_timeout = None

        self._launched_threads = []

        # Variables for timeout testing:
        self._last_activity_time = None
        self._client_activity_time = {}

        # Hub message id counter, used to create hub msg ids
        self._hub_msg_id_counter = 0

        # Hub secret code
        self._hub_secret_code_customized = secret
        self._hub_secret = self._create_secret_code()

        # Hub public id (as SAMP client)
        self._hub_public_id = ""

        # Client ids
        # {private_key: (public_id, timestamp)}
        self._private_keys = {}

        # Metadata per client
        # {private_key: metadata}
        self._metadata = {}

        # List of subscribed clients per MType
        # {mtype: private_key list}
        self._mtype2ids = {}

        # List of subscribed MTypes per client
        # {private_key: mtype list}
        self._id2mtypes = {}

        # List of XML-RPC addresses per client
        # {public_id: (XML-RPC address, ServerProxyPool instance)}
        self._xmlrpc_endpoints = {}

        # Synchronous message id heap
        self._sync_msg_ids_heap = {}

        # Public ids counter
        self._client_id_counter = -1

        # Standard Profile only operations
        self._server.register_function(self._ping, 'samp.hub.ping')
        self._server.register_function(self._set_xmlrpc_callback, 'samp.hub.setXmlrpcCallback')

        # Standard API operations
        self._server.register_function(self._register, 'samp.hub.register')
        self._server.register_function(self._unregister, 'samp.hub.unregister')
        self._server.register_function(self._declare_metadata, 'samp.hub.declareMetadata')
        self._server.register_function(self._get_metadata, 'samp.hub.getMetadata')
        self._server.register_function(self._declare_subscriptions, 'samp.hub.declareSubscriptions')
        self._server.register_function(self._get_subscriptions, 'samp.hub.getSubscriptions')
        self._server.register_function(self._get_registered_clients, 'samp.hub.getRegisteredClients')
        self._server.register_function(self._get_subscribed_clients, 'samp.hub.getSubscribedClients')
        self._server.register_function(self._notify, 'samp.hub.notify')
        self._server.register_function(self._notify_all, 'samp.hub.notifyAll')
        self._server.register_function(self._call, 'samp.hub.call')
        self._server.register_function(self._call_all, 'samp.hub.callAll')
        self._server.register_function(self._call_and_wait, 'samp.hub.callAndWait')
        self._server.register_function(self._reply, 'samp.hub.reply')

        if web_profile:

            # Web Profile methods like Standard Profile
            self._web_profile_server.register_function(self._ping, 'samp.webhub.ping')
            self._web_profile_server.register_function(self._unregister, 'samp.webhub.unregister')
            self._web_profile_server.register_function(self._declare_metadata, 'samp.webhub.declareMetadata')
            self._web_profile_server.register_function(self._get_metadata, 'samp.webhub.getMetadata')
            self._web_profile_server.register_function(self._declare_subscriptions, 'samp.webhub.declareSubscriptions')
            self._web_profile_server.register_function(self._get_subscriptions, 'samp.webhub.getSubscriptions')
            self._web_profile_server.register_function(self._get_registered_clients, 'samp.webhub.getRegisteredClients')
            self._web_profile_server.register_function(self._get_subscribed_clients, 'samp.webhub.getSubscribedClients')
            self._web_profile_server.register_function(self._notify, 'samp.webhub.notify')
            self._web_profile_server.register_function(self._notify_all, 'samp.webhub.notifyAll')
            self._web_profile_server.register_function(self._call, 'samp.webhub.call')
            self._web_profile_server.register_function(self._call_all, 'samp.webhub.callAll')
            self._web_profile_server.register_function(self._call_and_wait, 'samp.webhub.callAndWait')
            self._web_profile_server.register_function(self._reply, 'samp.webhub.reply')

            # Methods peculiar for Web Profile
            self._web_profile_server.register_function(self._web_profile_register, 'samp.webhub.register')
            self._web_profile_server.register_function(self._web_profile_allowReverseCallbacks, 'samp.webhub.allowReverseCallbacks')
            self._web_profile_server.register_function(self._web_profile_pullCallbacks, 'samp.webhub.pullCallbacks')

    @property
    def id(self):
        """
        The unique hub ID.
        """
        return self._id

    def _launch_thread(self, group=None, target=None, name=None, args=None):

        # Remove inactive threads
        remove = []
        for t in self._launched_threads:
            if not t.is_alive():
                remove.append(t)
        for t in remove:
            self._launched_threads.remove(t)

        # Start new thread
        t = threading.Thread(group=group, target=target, name=name, args=args)
        t.start()

        # Add to list of launched threads
        self._launched_threads.append(t)

    def _join_launched_threads(self, timeout=None):
        for t in self._launched_threads:
            t.join(timeout=timeout)

    def _timeout_test_hub(self):

        if self._timeout == 0:
            return

        last = time.time()
        while self._is_running:
            time.sleep(0.05)  # keep this small to check _is_running often
            now = time.time()
            if now - last > 1.:
                with self._thread_lock:
                    if self._last_activity_time is not None:
                        if now - self._last_activity_time >= self._timeout:
                            warnings.warn("Timeout expired, Hub is shutting down!",
                                          SAMPWarning)
                            self.stop()
                            return
                last = now

    def _timeout_test_client(self):

        if self._client_timeout == 0:
            return

        last = time.time()
        while self._is_running:
            time.sleep(0.05)  # keep this small to check _is_running often
            now = time.time()
            if now - last > 1.:
                for private_key in self._client_activity_time.keys():
                    if (now - self._client_activity_time[private_key] > self._client_timeout
                        and private_key != self._hub_private_key):
                        warnings.warn("Client %s timeout expired!"
                                      % private_key,
                                      SAMPWarning)
                        self._notify_disconnection(private_key)
                        self._unregister(private_key)
                last = now

    def _hub_as_client_request_handler(self, method, args):
        if method == 'samp.client.receiveCall':
            return self._receive_call(*args)
        elif method == 'samp.client.receiveNotification':
            return self._receive_notification(*args)
        elif method == 'samp.client.receiveResponse':
            return self._receive_response(*args)
        elif method == 'samp.app.ping':
            return self._ping(*args)

    def _setup_hub_as_client(self):

        hub_metadata = {"samp.name": "Astropy SAMP Hub",
                        "samp.description.text": self._label,
                        "author.name": "The Astropy Collaboration",
                        "samp.documentation.url": "http://docs.astropy.org/en/stable/vo/samp",
                        "samp.icon.url": self._url + "/samp/icon"}

        result = self._register(self._hub_secret)
        self._hub_public_id = result["samp.self-id"]
        self._hub_private_key = result["samp.private-key"]
        self._set_xmlrpc_callback(self._hub_private_key, self._url)
        self._declare_metadata(self._hub_private_key, hub_metadata)
        self._declare_subscriptions(self._hub_private_key,
                                    {"samp.app.ping": {},
                                     "x-samp.query.by-meta": {}})

    def start(self, wait=False):
        """
        Start the current SAMP Hub instance and create the lock file. Hub
        start-up can be blocking or non blocking depending on the ``wait``
        parameter.

        Parameters
        ----------
        wait : bool
            If `True` then the Hub process is joined with the caller, blocking
            the code flow. Usually `True` option is used to run a stand-alone
            Hub in an executable script. If `False` (default), then the Hub
            process runs in a separated thread. `False` is usually used in a
            Python shell.
        """
        if self._is_running:
            raise SAMPHubError("Hub is already running")

        if self._lockfile is not None:
            raise SAMPHubError("Hub is not running but lockfile is set")

        self._lockfile = create_lock_file(lockfilename=self._customlockfilename,
                                          mode=self._mode, hub_id=self.id,
                                          hub_params=self.params)

        self._is_running = True
        self._update_last_activity_time()

        self._setup_hub_as_client()
        self._start_threads()

        log.info("Hub started")

        if wait and self._is_running:
            self._thread_run.join()

    @property
    def params(self):
        """
        The hub parameters (which are written to the logfile)
        """

        params = {}

        # Keys required by standard profile

        params['samp.secret'] = self._hub_secret
        params['samp.hub.xmlrpc.url'] = self._url
        params['samp.profile.version'] = __profile_version__

        # Custom keys

        params['hub.id'] = self.id
        params['hub.label'] = self._label or "Hub {0}".format(self.id)

        if SSL_SUPPORT and self._https:

            # Certificate request
            cert_reqs_types = ["NONE", "OPTIONAL", "REQUIRED"]
            params['hub.ssl.certificate'] = cert_reqs_types[self._cert_reqs]

            # SSL protocol version
            ssl_protocol_types = ["SSLv2", "SSLv3", "SSLv23", "TLSv1"]
            params['hub.ssl.protocol'] = ssl_protocol_types[self._ssl_version]

        return params

    def _start_threads(self):
        self._thread_run.start()
        if self._thread_hub_timeout is not None:
            self._thread_hub_timeout.start()
        if self._thread_client_timeout is not None:
            self._thread_client_timeout.start()

    def _create_secret_code(self):
        if self._hub_secret_code_customized is not None:
            return self._hub_secret_code_customized
        else:
            return str(uuid.uuid1())

    def stop(self):
        """
        Stop the current SAMP Hub instance and delete the lock file.
        """

        if not self._is_running:
            return

        log.info("Hub is stopping...")

        self._notify_shutdown()

        self._is_running = False

        if (os.path.isfile(self._lockfile)):
            lockfiledict = read_lockfile(self._lockfile)
            if lockfiledict['samp.secret'] == self._hub_secret:
                os.remove(self._lockfile)
        self._lockfile = None

        # Reset vaiables
        self._join_all_threads(timeout=10.)

        self._hub_msg_id_counter = 0
        self._hub_secret = self._create_secret_code()
        self._hub_public_id = ""
        self._metadata = {}
        self._private_keys = {}
        self._mtype2ids = {}
        self._id2mtypes = {}
        self._xmlrpc_endpoints = {}
        self._last_activity_time = None

        log.info("Hub stopped.")

    def _join_all_threads(self, timeout=None):
        # In some cases, ``stop`` may be called from some of the sub-threads,
        # so we just need to make sure that we don't try and shut down the
        # calling thread.
        current_thread = threading.current_thread()
        if self._thread_run is not current_thread:
            self._thread_run.join(timeout=timeout)
        if self._thread_hub_timeout is not None and self._thread_hub_timeout is not current_thread:
            self._thread_hub_timeout.join(timeout=timeout)
        if self._thread_client_timeout is not None and self._thread_client_timeout is not current_thread:
            self._thread_client_timeout.join(timeout=timeout)
        self._join_launched_threads(timeout=timeout)

    @property
    def is_running(self):
        """Return an information concerning the Hub running status.

        Returns
        -------
        running : bool
            Is the hub running?
        """
        return self._is_running

    def _serve_forever(self):

        while self._is_running:

            try:
                read_ready = select.select([self._server.socket], [], [], 0.01)[0]
            except select.error as exc:
                warnings.warn("Call to select() in SAMPHubServer failed: {0}".format(exc),
                              SAMPWarning)
            else:
                if read_ready:
                    self._server.handle_request()

            if self._web_profile:

                # We now check if there are any connection requests from the
                # web profile, and if so, we initialize the pop-up.
                if self._web_profile_dialog is None:
                    try:
                        request = self._web_profile_requests_queue.get_nowait()
                    except queue.Empty:
                        pass
                    else:
                        web_profile_text_dialog(request, self._web_profile_requests_result)

                # We now check for requests over the web profile socket, and we
                # also update the pop-up in case there are any changes.
                try:
                    read_ready = select.select([self._web_profile_server.socket], [], [], 0.01)[0]
                except select.error as exc:
                    warnings.warn("Call to select() in SAMPHubServer failed: {0}".format(exc),
                                  SAMPWarning)
                else:
                    if read_ready:
                        self._web_profile_server.handle_request()

    def _notify_shutdown(self):
        msubs = SAMPHubServer.get_mtype_subtypes("samp.hub.event.shutdown")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                for key in self._mtype2ids[mtype]:
                    self._notify_(self._hub_private_key,
                                  self._private_keys[key][0],
                                  {"samp.mtype": "samp.hub.event.shutdown",
                                   "samp.params": {}})

    def _notify_register(self, private_key):
        msubs = SAMPHubServer.get_mtype_subtypes("samp.hub.event.register")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                public_id = self._private_keys[private_key][0]
                for key in self._mtype2ids[mtype]:
                    # if key != private_key:
                    self._notify(self._hub_private_key,
                                 self._private_keys[key][0],
                                 {"samp.mtype": "samp.hub.event.register",
                                  "samp.params": {"id": public_id}})

    def _notify_unregister(self, private_key):
        msubs = SAMPHubServer.get_mtype_subtypes("samp.hub.event.unregister")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                public_id = self._private_keys[private_key][0]
                for key in self._mtype2ids[mtype]:
                    if key != private_key:
                        self._notify(self._hub_private_key,
                                     self._private_keys[key][0],
                                     {"samp.mtype": "samp.hub.event.unregister",
                                      "samp.params": {"id": public_id}})

    def _notify_metadata(self, private_key):
        msubs = SAMPHubServer.get_mtype_subtypes("samp.hub.event.metadata")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                public_id = self._private_keys[private_key][0]
                for key in self._mtype2ids[mtype]:
                    # if key != private_key:
                    self._notify(self._hub_private_key,
                                 self._private_keys[key][0],
                                 {"samp.mtype": "samp.hub.event.metadata",
                                  "samp.params": {"id": public_id,
                                                  "metadata": self._metadata[private_key]}
                                  })

    def _notify_subscriptions(self, private_key):
        msubs = SAMPHubServer.get_mtype_subtypes("samp.hub.event.subscriptions")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                public_id = self._private_keys[private_key][0]
                for key in self._mtype2ids[mtype]:
                    self._notify(self._hub_private_key,
                                 self._private_keys[key][0],
                                 {"samp.mtype": "samp.hub.event.subscriptions",
                                  "samp.params": {"id": public_id,
                                                  "subscriptions": self._id2mtypes[private_key]}
                                  })

    def _notify_disconnection(self, private_key):

        def _xmlrpc_call_disconnect(endpoint, private_key, hub_public_id, message):
            endpoint.samp.client.receiveNotification(private_key, hub_public_id, message)

        msubs = SAMPHubServer.get_mtype_subtypes("samp.hub.disconnect")
        public_id = self._private_keys[private_key][0]
        endpoint = self._xmlrpc_endpoints[public_id][1]

        for mtype in msubs:
            if mtype in self._mtype2ids and private_key in self._mtype2ids[mtype]:
                log.debug("notify disconnection to %s" % (public_id))
                self._launch_thread(target=_xmlrpc_call_disconnect,
                                   args=(endpoint, private_key,
                                         self._hub_public_id,
                                         {"samp.mtype": "samp.hub.disconnect",
                                          "samp.params": {"reason": "Timeout expired!"}}))

    def _ping(self):
        self._update_last_activity_time()
        log.debug("ping")
        return "1"

    def _query_by_metadata(self, key, value):
        public_id_list = []
        for private_id in self._metadata:
            if key in self._metadata[private_id]:
                if self._metadata[private_id][key] == value:
                    public_id_list.append(self._private_keys[private_id][0])

        return public_id_list

    def _set_xmlrpc_callback(self, private_key, xmlrpc_addr):
        self._update_last_activity_time(private_key)
        if private_key in self._private_keys:
            if private_key == self._hub_private_key:
                public_id = self._private_keys[private_key][0]
                self._xmlrpc_endpoints[public_id] = \
                    (xmlrpc_addr, _HubAsClient(self._hub_as_client_request_handler))
                return ""

            # Dictionary stored with the public id

            log.debug("set_xmlrpc_callback: %s %s" % (private_key,
                                                      xmlrpc_addr))

            server_proxy_pool = None
            if SSL_SUPPORT and xmlrpc_addr[0:5] == "https":

                transport = SafeTransport(key_file=self._key_file,
                                          cert_file=self._cert_file,
                                          cert_reqs=self._cert_reqs,
                                          ca_certs=self._ca_certs,
                                          ssl_version=self._ssl_version)

                server_proxy_pool = ServerProxyPool(self._pool_size,
                                                    xmlrpc.ServerProxy,
                                                    xmlrpc_addr,
                                                    transport=transport,
                                                    allow_none=1)
            else:

                server_proxy_pool = ServerProxyPool(self._pool_size,
                                                    xmlrpc.ServerProxy,
                                                    xmlrpc_addr, allow_none=1)

            public_id = self._private_keys[private_key][0]
            self._xmlrpc_endpoints[public_id] = (xmlrpc_addr,
                                                server_proxy_pool)
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

        return ""

    def _perform_standard_register(self):

        with self._thread_lock:
            private_key, public_id = self._get_new_ids()
        self._private_keys[private_key] = (public_id, time.time())
        self._update_last_activity_time(private_key)
        self._notify_register(private_key)
        log.debug("register: private-key = %s and self-id = %s"
                  % (private_key, public_id))
        return {"samp.self-id": public_id,
                "samp.private-key": private_key,
                "samp.hub-id": self._hub_public_id}

    def _register(self, secret):
        self._update_last_activity_time()
        if secret == self._hub_secret:
            return self._perform_standard_register()
        else:
            # return {"samp.self-id": "", "samp.private-key": "", "samp.hub-id": ""}
            raise SAMPProxyError(7, "Bad secret code")

    def _get_new_ids(self):
        private_key = str(uuid.uuid1())
        self._client_id_counter += 1
        public_id = 'cli#hub'
        if self._client_id_counter > 0:
            public_id = "cli#%d" % (self._client_id_counter)

        return private_key, public_id

    def _unregister(self, private_key):

        self._update_last_activity_time()

        public_key = ""

        self._notify_unregister(private_key)

        with self._thread_lock:

            if private_key in self._private_keys:
                public_key = self._private_keys[private_key][0]
                del self._private_keys[private_key]
            else:
                return ""

            if private_key in self._metadata:
                del self._metadata[private_key]

            if private_key in self._id2mtypes:
                del self._id2mtypes[private_key]

            for mtype in self._mtype2ids.keys():
                if private_key in self._mtype2ids[mtype]:
                    self._mtype2ids[mtype].remove(private_key)

            if public_key in self._xmlrpc_endpoints:
                del self._xmlrpc_endpoints[public_key]

            if private_key in self._client_activity_time:
                del self._client_activity_time[private_key]

            if self._web_profile:
                if private_key in self._web_profile_callbacks:
                    del self._web_profile_callbacks[private_key]
                self._web_profile_server.remove_client(private_key)

        log.debug("unregister %s (%s)" % (public_key, private_key))

        return ""

    def _declare_metadata(self, private_key, metadata):
        self._update_last_activity_time(private_key)
        if private_key in self._private_keys:
            log.debug("declare_metadata: private-key = %s metadata = %s"
                      % (private_key, str(metadata)))
            self._metadata[private_key] = metadata
            self._notify_metadata(private_key)
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)
        return ""

    def _get_metadata(self, private_key, client_id):
        self._update_last_activity_time(private_key)
        if private_key in self._private_keys:
            client_private_key = self._public_id_to_private_key(client_id)
            log.debug("get_metadata: private-key = %s client-id = %s"
                      % (private_key, client_id))
            if client_private_key is not None:
                if client_private_key in self._metadata:
                    log.debug("--> metadata = %s"
                              % self._metadata[client_private_key])
                    return self._metadata[client_private_key]
                else:
                    return {}
            else:
                raise SAMPProxyError(6, "Invalid client ID")
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    def _declare_subscriptions(self, private_key, mtypes):

        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:

            log.debug("declare_subscriptions: private-key = %s mtypes = %s"
                      % (private_key, str(mtypes)))

            # remove subscription to previous mtypes
            if private_key in self._id2mtypes:

                prev_mtypes = self._id2mtypes[private_key]

                for mtype in prev_mtypes:
                    try:
                        self._mtype2ids[mtype].remove(private_key)
                    except ValueError:  # private_key is not in list
                        pass

            self._id2mtypes[private_key] = copy.deepcopy(mtypes)

            # remove duplicated MType for wildcard overwriting
            original_mtypes = copy.deepcopy(mtypes)

            for mtype in original_mtypes:
                if mtype.endswith("*"):
                    for mtype2 in original_mtypes:
                        if mtype2.startswith(mtype[:-1]) and \
                           mtype2 != mtype:
                            if mtype2 in mtypes:
                                del(mtypes[mtype2])

            log.debug("declare_subscriptions: subscriptions accepted from %s => %s"
                      % (private_key, str(mtypes)))

            for mtype in mtypes:

                if mtype in self._mtype2ids:
                    if not private_key in self._mtype2ids[mtype]:
                        self._mtype2ids[mtype].append(private_key)
                else:
                    self._mtype2ids[mtype] = [private_key]

            self._notify_subscriptions(private_key)

        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

        return ""

    def _get_subscriptions(self, private_key, client_id):

        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:
            client_private_key = self._public_id_to_private_key(client_id)
            if client_private_key is not None:
                if client_private_key in self._id2mtypes:
                    log.debug("get_subscriptions: client-id = %s mtypes = %s"
                              % (client_id,
                                 str(self._id2mtypes[client_private_key])))
                    return self._id2mtypes[client_private_key]
                else:
                    log.debug("get_subscriptions: client-id = %s mtypes = "
                              "missing" % client_id)
                    return {}
            else:
                raise SAMPProxyError(6, "Invalid client ID")
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    def _get_registered_clients(self, private_key):

        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:
            reg_clients = []
            for pkey in self._private_keys.keys():
                if pkey != private_key:
                    reg_clients.append(self._private_keys[pkey][0])
            log.debug("get_registered_clients: private_key = %s clients = %s"
                      % (private_key, reg_clients))
            return reg_clients
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    def _get_subscribed_clients(self, private_key, mtype):

        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:
            sub_clients = {}

            for pkey in self._private_keys.keys():
                if pkey != private_key and self._is_subscribed(pkey, mtype):
                    sub_clients[self._private_keys[pkey][0]] = {}

            log.debug("get_subscribed_clients: private_key = %s mtype = %s clients = %s"
                      % (private_key, mtype, sub_clients))
            return sub_clients
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    @staticmethod
    def get_mtype_subtypes(mtype):
        """
        Return a list containing all the possible wildcarded subtypes of MType.

        Parameters
        ----------
        mtype : str
            MType to be parsed.

        Returns
        -------
        types : list
            List of subtypes

        Examples
        --------
        >>> from astropy.vo.samp import SAMPHubServer
        >>> SAMPHubServer.get_mtype_subtypes("samp.app.ping")
        ['samp.app.ping', 'samp.app.*', 'samp.*', '*']
        """

        subtypes = []

        msubs = mtype.split(".")
        indexes = list(range(len(msubs)))
        indexes.reverse()
        indexes.append(-1)

        for i in indexes:
            tmp_mtype = ".".join(msubs[:i + 1])
            if tmp_mtype != mtype:
                if tmp_mtype != "":
                    tmp_mtype = tmp_mtype + ".*"
                else:
                    tmp_mtype = "*"
            subtypes.append(tmp_mtype)

        return subtypes

    def _is_subscribed(self, private_key, mtype):

        subscribed = False

        msubs = SAMPHubServer.get_mtype_subtypes(mtype)

        for msub in msubs:
            if msub in self._mtype2ids:
                if private_key in self._mtype2ids[msub]:
                    subscribed = True

        return subscribed

    def _notify(self, private_key, recipient_id, message):
        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:
            if self._is_subscribed(self._public_id_to_private_key(recipient_id),
                                   message["samp.mtype"]) is False:
                raise SAMPProxyError(2, "Client %s not subscribed to MType %s"
                                        % (recipient_id, message["samp.mtype"]))

            self._launch_thread(target=self._notify_, args=(private_key,
                                                            recipient_id,
                                                            message))
            return {}
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    def _notify_(self, sender_private_key, recipient_public_id, message):

        if sender_private_key not in self._private_keys:
            return

        sender_public_id = self._private_keys[sender_private_key][0]

        try:

            log.debug("notify %s from %s to %s"
                      % (message["samp.mtype"], sender_public_id,
                         recipient_public_id))

            recipient_private_key = self._public_id_to_private_key(recipient_public_id)

            if recipient_private_key is None:
                raise SAMPHubError("Invalid client ID")

            for attempt in range(10):

                if not self._is_running:
                    time.sleep(0.01)
                    continue

                try:

                    if (self._web_profile and
                        recipient_private_key in self._web_profile_callbacks):

                        # Web Profile
                        callback = {"samp.methodName": "receiveNotification",
                                    "samp.params": [sender_public_id, message]}
                        self._web_profile_callbacks[recipient_private_key].put(callback)

                    else:

                        # Standard Profile
                        hub = self._xmlrpc_endpoints[recipient_public_id][1]
                        hub.samp.client.receiveNotification(recipient_private_key,
                                                            sender_public_id,
                                                            message)

                except xmlrpc.Fault as exc:
                    log.debug("%s XML-RPC endpoint error (attempt %d): %s"
                              % (recipient_public_id, attempt + 1,
                                 exc.faultString))
                    time.sleep(0.01)
                else:
                    return

            # If we are here, then the above attempts failed
            raise SAMPHubError("notification failed after 10 attempts")

        except Exception as exc:
            warnings.warn("%s notification from client %s to client %s failed [%s]"
                          % (message["samp.mtype"], sender_public_id,
                             recipient_public_id, exc),
                          SAMPWarning)

    def _notify_all(self, private_key, message):
        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:
            if not "samp.mtype" in message:
                raise SAMPProxyError(3, "samp.mtype keyword is missing")
            recipient_ids = self._notify_all_(private_key, message)
            return recipient_ids
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    def _notify_all_(self, sender_private_key, message):

        recipient_ids = []
        msubs = SAMPHubServer.get_mtype_subtypes(message["samp.mtype"])

        for mtype in msubs:
            if mtype in self._mtype2ids:
                for key in self._mtype2ids[mtype]:
                    if key != sender_private_key:
                        _recipient_id = self._private_keys[key][0]
                        recipient_ids.append(_recipient_id)
                        self._launch_thread(target=self._notify,
                                         args=(sender_private_key,
                                               _recipient_id, message)
                                         )

        return recipient_ids

    def _call(self, private_key, recipient_id, msg_tag, message):
        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:
            if self._is_subscribed(self._public_id_to_private_key(recipient_id),
                                   message["samp.mtype"]) is False:
                raise SAMPProxyError(2, "Client %s not subscribed to MType %s"
                                        % (recipient_id, message["samp.mtype"]))
            public_id = self._private_keys[private_key][0]
            msg_id = self._get_new_hub_msg_id(public_id, msg_tag)
            self._launch_thread(target=self._call_, args=(private_key, public_id,
                                                          recipient_id, msg_id,
                                                          message))
            return msg_id
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    def _call_(self, sender_private_key, sender_public_id,
               recipient_public_id, msg_id, message):

        if sender_private_key not in self._private_keys:
            return

        try:

            log.debug("call %s from %s to %s (%s)"
                      % (msg_id.split(";;")[0], sender_public_id,
                         recipient_public_id, message["samp.mtype"]))

            recipient_private_key = self._public_id_to_private_key(recipient_public_id)

            if recipient_private_key is None:
                raise SAMPHubError("Invalid client ID")

            for attempt in range(10):

                if not self._is_running:
                    time.sleep(0.01)
                    continue

                try:

                    if (self._web_profile and
                        recipient_private_key in self._web_profile_callbacks):

                        # Web Profile
                        callback = {"samp.methodName": "receiveCall",
                                    "samp.params": [sender_public_id, msg_id, message]}
                        self._web_profile_callbacks[recipient_private_key].put(callback)

                    else:

                        # Standard Profile
                        hub = self._xmlrpc_endpoints[recipient_public_id][1]
                        hub.samp.client.receiveCall(recipient_private_key,
                                                    sender_public_id, msg_id,
                                                    message)

                except xmlrpc.Fault as exc:
                    log.debug("%s XML-RPC endpoint error (attempt %d): %s"
                              % (recipient_public_id, attempt + 1,
                                 exc.faultString))
                    time.sleep(0.01)
                else:
                    return

            # If we are here, then the above attempts failed
            raise SAMPHubError("call failed after 10 attempts")

        except Exception as exc:
            warnings.warn("%s call %s from client %s to client %s failed [%s,%s]"
                          % (message["samp.mtype"], msg_id.split(";;")[0],
                             sender_public_id, recipient_public_id, type(exc), exc),
                          SAMPWarning)

    def _call_all(self, private_key, msg_tag, message):
        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:
            if not "samp.mtype" in message:
                raise SAMPProxyError(3, "samp.mtype keyword is missing in "
                                        "message tagged as %s" % msg_tag)

            public_id = self._private_keys[private_key][0]
            msg_id = self._call_all_(private_key, public_id, msg_tag, message)
            return msg_id
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    def _call_all_(self, sender_private_key, sender_public_id, msg_tag,
                   message):

        msg_id = {}
        msubs = SAMPHubServer.get_mtype_subtypes(message["samp.mtype"])

        for mtype in msubs:
            if mtype in self._mtype2ids:
                for key in self._mtype2ids[mtype]:
                    if key != sender_private_key:
                        _msg_id = self._get_new_hub_msg_id(sender_public_id,
                                                           msg_tag)
                        receiver_public_id = self._private_keys[key][0]
                        msg_id[receiver_public_id] = _msg_id
                        self._launch_thread(target=self._call_,
                                            args=(sender_private_key,
                                                  sender_public_id,
                                                  receiver_public_id, _msg_id,
                                                  message))
        return msg_id

    def _call_and_wait(self, private_key, recipient_id, message, timeout):
        self._update_last_activity_time(private_key)

        if private_key in self._private_keys:
            timeout = int(timeout)

            now = time.time()
            response = {}

            msg_id = self._call(private_key, recipient_id, "samp::sync::call",
                                message)
            self._sync_msg_ids_heap[msg_id] = None

            while self._is_running:
                if timeout > 0 and time.time() - now >= timeout:
                    del(self._sync_msg_ids_heap[msg_id])
                    raise SAMPProxyError(1, "Timeout expired!")

                if self._sync_msg_ids_heap[msg_id] is not None:
                    response = copy.deepcopy(self._sync_msg_ids_heap[msg_id])
                    del(self._sync_msg_ids_heap[msg_id])
                    break
                time.sleep(0.01)

            return response
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

    def _reply(self, private_key, msg_id, response):
        """
        The main method that gets called for replying. This starts up an
        asynchronous reply thread and returns.
        """
        self._update_last_activity_time(private_key)
        if private_key in self._private_keys:
            self._launch_thread(target=self._reply_, args=(private_key, msg_id,
                                                           response))
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)

        return {}

    def _reply_(self, responder_private_key, msg_id, response):

        if responder_private_key not in self._private_keys or not msg_id:
            return

        responder_public_id = self._private_keys[responder_private_key][0]
        counter, hub_public_id, recipient_public_id, recipient_msg_tag = msg_id.split(";;", 3)

        try:

            log.debug("reply %s from %s to %s"
                      % (counter, responder_public_id, recipient_public_id))

            if recipient_msg_tag == "samp::sync::call":

                if msg_id in self._sync_msg_ids_heap.keys():
                    self._sync_msg_ids_heap[msg_id] = response

            else:

                recipient_private_key = self._public_id_to_private_key(recipient_public_id)

                if recipient_private_key is None:
                    raise SAMPHubError("Invalid client ID")

                for attempt in range(10):

                    if not self._is_running:
                        time.sleep(0.01)
                        continue

                    try:

                        if (self._web_profile and
                            recipient_private_key in self._web_profile_callbacks):

                            # Web Profile
                            callback = {"samp.methodName": "receiveResponse",
                                        "samp.params": [responder_public_id,
                                                        recipient_msg_tag,
                                                        response]}
                            self._web_profile_callbacks[recipient_private_key].put(callback)

                        else:

                            # Standard Profile
                            hub = self._xmlrpc_endpoints[recipient_public_id][1]
                            hub.samp.client.receiveResponse(recipient_private_key,
                                                            responder_public_id,
                                                            recipient_msg_tag,
                                                            response)

                    except xmlrpc.Fault as exc:
                        log.debug("%s XML-RPC endpoint error (attempt %d): %s"
                                  % (recipient_public_id, attempt + 1,
                                     exc.faultString))
                        time.sleep(0.01)
                    else:
                        return

                # If we are here, then the above attempts failed
                raise SAMPHubError("reply failed after 10 attempts")

        except Exception as exc:
            warnings.warn("%s reply from client %s to client %s failed [%s]"
                          % (recipient_msg_tag, responder_public_id,
                             recipient_public_id, exc),
                          SAMPWarning)

    def _public_id_to_private_key(self, public_id):

        for private_key in self._private_keys.keys():
            if self._private_keys[private_key][0] == public_id:
                return private_key
        return None

    def _get_new_hub_msg_id(self, sender_public_id, sender_msg_id):
        with self._thread_lock:
            self._hub_msg_id_counter += 1
        return "msg#%d;;%s;;%s;;%s" % \
               (self._hub_msg_id_counter, self._hub_public_id,
                sender_public_id, sender_msg_id)

    def _update_last_activity_time(self, private_key=None):
        with self._thread_lock:
            self._last_activity_time = time.time()
            if private_key is not None:
                self._client_activity_time[private_key] = time.time()

    def _receive_notification(self, private_key, sender_id, message):
        return ""

    def _receive_call(self, private_key, sender_id, msg_id, message):
        if private_key == self._hub_private_key:

            if "samp.mtype" in message and message["samp.mtype"] == "samp.app.ping":
                self._reply(self._hub_private_key, msg_id,
                            {"samp.status": SAMP_STATUS_OK, "samp.result": {}})

            elif ("samp.mtype" in message and
                 (message["samp.mtype"] == "x-samp.query.by-meta" or
                  message["samp.mtype"] == "samp.query.by-meta")):

                ids_list = self._query_by_metadata(message["samp.params"]["key"],
                                                   message["samp.params"]["value"])
                self._reply(self._hub_private_key, msg_id,
                            {"samp.status": SAMP_STATUS_OK,
                             "samp.result": {"ids": ids_list}})

            return ""
        else:
            return ""

    def _receive_response(self, private_key, responder_id, msg_tag, response):
        return ""

    def _web_profile_register(self, identity_info,
                              client_address=("unknown", 0),
                              origin="unknown"):

        self._update_last_activity_time()

        if not client_address[0] in ["localhost", "127.0.0.1"]:
            raise SAMPProxyError(403, "Request of registration rejected "
                                      "by the Hub.")

        if not origin:
            origin = "unknown"

        if isinstance(identity_info, dict):
            # an old version of the protocol provided just a string with the app name
            if "samp.name" not in identity_info:
                raise SAMPProxyError(403, "Request of registration rejected "
                                          "by the Hub (application name not "
                                          "provided).")

        # Red semaphore for the other threads
        self._web_profile_requests_semaphore.put("wait")
        # Set the request to be displayed for the current thread
        self._web_profile_requests_queue.put((identity_info, client_address,
                                              origin))
        # Get the popup dialogue response
        response = self._web_profile_requests_result.get()
        # OK, semaphore green
        self._web_profile_requests_semaphore.get()

        if response:
            register_map = self._perform_standard_register()
            translator_url = ("http://localhost:%d/translator/%s?ref="
                              % (self._web_port, register_map["samp.private-key"]))
            register_map["samp.url-translator"] = translator_url
            self._web_profile_server.add_client(register_map["samp.private-key"])
            return register_map
        else:
            raise SAMPProxyError(403, "Request of registration rejected by "
                                      "the user.")

    def _web_profile_allowReverseCallbacks(self, private_key, allow):
        self._update_last_activity_time()
        if private_key in self._private_keys:
            if allow == "0":
                if private_key in self._web_profile_callbacks:
                    del self._web_profile_callbacks[private_key]
            else:
                self._web_profile_callbacks[private_key] = queue.Queue()
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)
        return ""

    def _web_profile_pullCallbacks(self, private_key, timeout_secs):
        self._update_last_activity_time()
        if private_key in self._private_keys:
            callback = []
            callback_queue = self._web_profile_callbacks[private_key]
            try:
                while self._is_running:
                    item_queued = callback_queue.get_nowait()
                    callback.append(item_queued)
            except queue.Empty:
                pass
            return callback
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid."
                                    % private_key)


class WebProfileDialog(object):
    """
    A base class to make writing Web Profile GUI consent dialogs
    easier.

    The concrete class must:

        1) Poll ``handle_queue`` periodically, using the timer services
           of the GUI's event loop.  This function will call
           ``self.show_dialog`` when a request requires authorization.
           ``self.show_dialog`` will be given the arguments:

              - ``samp_name``: The name of the application making the request.

              - ``details``: A dictionary of details about the client
                making the request.

              - ``client``: A hostname, port pair containing the client
                address.

              - ``origin``: A string containing the origin of the
                request.

        2) Call ``consent`` or ``reject`` based on the user's response to
           the dialog.
    """

    def handle_queue(self):
        try:
            request = self.queue_request.get_nowait()
        except queue.Empty:  # queue is set but empty
            pass
        except AttributeError:  # queue has not been set yet
            pass
        else:
            if isinstance(request[0], str):  # To support the old protocol version
                samp_name = request[0]
            else:
                samp_name = request[0]["samp.name"]

            self.show_dialog(samp_name, request[0], request[1], request[2])

    def consent(self):
        self.queue_result.put(True)

    def reject(self):
        self.queue_result.put(False)

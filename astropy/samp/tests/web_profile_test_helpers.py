import time
import threading
import xmlrpc.client as xmlrpc

from ..hub import WebProfileDialog
from ..hub_proxy import SAMPHubProxy
from ..client import SAMPClient
from ..integrated_client import SAMPIntegratedClient
from ..utils import ServerProxyPool
from ..errors import SAMPClientError, SAMPHubError


class AlwaysApproveWebProfileDialog(WebProfileDialog):

    def __init__(self):
        self.polling = True
        WebProfileDialog.__init__(self)

    def show_dialog(self, *args):
        self.consent()

    def poll(self):
        while self.polling:
            self.handle_queue()
            time.sleep(0.1)

    def stop(self):
        self.polling = False


class SAMPWebHubProxy(SAMPHubProxy):
    """
    Proxy class to simplify the client interaction with a SAMP hub (via the web
    profile).

    In practice web clients should run from the browser, so this is provided as
    a means of testing a hub's support for the web profile from Python.
    """

    def connect(self, pool_size=20, web_port=21012):
        """
        Connect to the current SAMP Hub on localhost:web_port

        Parameters
        ----------
        pool_size : int, optional
            The number of socket connections opened to communicate with the
            Hub.
        """

        self._connected = False

        try:
            self.proxy = ServerProxyPool(pool_size, xmlrpc.ServerProxy,
                                         'http://127.0.0.1:{0}'.format(web_port),
                                         allow_none=1)
            self.ping()
            self._connected = True
        except xmlrpc.ProtocolError as p:
            raise SAMPHubError("Protocol Error {}: {}".format(p.errcode, p.errmsg))

    @property
    def _samp_hub(self):
        """
        Property to abstract away the path to the hub, which allows this class
        to be used for both the standard and the web profile.
        """
        return self.proxy.samp.webhub

    def set_xmlrpc_callback(self, private_key, xmlrpc_addr):
        raise NotImplementedError("set_xmlrpc_callback is not defined for the "
                                  "web profile")

    def register(self, identity_info):
        """
        Proxy to ``register`` SAMP Hub method.
        """
        return self._samp_hub.register(identity_info)

    def allow_reverse_callbacks(self, private_key, allow):
        """
        Proxy to ``allowReverseCallbacks`` SAMP Hub method.
        """
        return self._samp_hub.allowReverseCallbacks(private_key, allow)

    def pull_callbacks(self, private_key, timeout):
        """
        Proxy to ``pullCallbacks`` SAMP Hub method.
        """
        return self._samp_hub.pullCallbacks(private_key, timeout)


class SAMPWebClient(SAMPClient):
    """
    Utility class which provides facilities to create and manage a SAMP
    compliant XML-RPC server that acts as SAMP callable web client application.

    In practice web clients should run from the browser, so this is provided as
    a means of testing a hub's support for the web profile from Python.

    Parameters
    ----------
    hub : :class:`~astropy.samp.hub_proxy.SAMPWebHubProxy`
        An instance of :class:`~astropy.samp.hub_proxy.SAMPWebHubProxy` to
        be used for messaging with the SAMP Hub.

    name : str, optional
        Client name (corresponding to ``samp.name`` metadata keyword).

    description : str, optional
        Client description (corresponding to ``samp.description.text`` metadata
        keyword).

    metadata : dict, optional
        Client application metadata in the standard SAMP format.

    callable : bool, optional
        Whether the client can receive calls and notifications. If set to
        `False`, then the client can send notifications and calls, but can not
        receive any.
    """

    def __init__(self, hub, name=None, description=None, metadata=None,
                 callable=True):

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

        self.hub = hub

        if self._callable:
            self._thread = threading.Thread(target=self._serve_forever)
            self._thread.daemon = True

    def _serve_forever(self):
        while self.is_running:
            # Watch for callbacks here
            if self._is_registered:
                results = self.hub.pull_callbacks(self.get_private_key(), 0)
                for result in results:
                    if result['samp.methodName'] == 'receiveNotification':
                        self.receive_notification(self._private_key,
                                                  *result['samp.params'])
                    elif result['samp.methodName'] == 'receiveCall':
                        self.receive_call(self._private_key,
                                          *result['samp.params'])
                    elif result['samp.methodName'] == 'receiveResponse':
                        self.receive_response(self._private_key,
                                              *result['samp.params'])

        self.hub.server_close()

    def register(self):
        """
        Register the client to the SAMP Hub.
        """
        if self.hub.is_connected:

            if self._private_key is not None:
                raise SAMPClientError("Client already registered")

            result = self.hub.register("Astropy SAMP Web Client")

            if result["samp.self-id"] == "":
                raise SAMPClientError("Registation failed - samp.self-id "
                                      "was not set by the hub.")

            if result["samp.private-key"] == "":
                raise SAMPClientError("Registation failed - samp.private-key "
                                      "was not set by the hub.")

            self._public_id = result["samp.self-id"]
            self._private_key = result["samp.private-key"]
            self._hub_id = result["samp.hub-id"]

            if self._callable:
                self._declare_subscriptions()
                self.hub.allow_reverse_callbacks(self._private_key, True)

            if self._metadata != {}:
                self.declare_metadata()

            self._is_registered = True

        else:
            raise SAMPClientError("Unable to register to the SAMP Hub. Hub "
                                  "proxy not connected.")


class SAMPIntegratedWebClient(SAMPIntegratedClient):
    """
    A Simple SAMP web client.

    In practice web clients should run from the browser, so this is provided as
    a means of testing a hub's support for the web profile from Python.

    This class is meant to simplify the client usage providing a proxy class
    that merges the :class:`~astropy.samp.client.SAMPWebClient` and
    :class:`~astropy.samp.hub_proxy.SAMPWebHubProxy` functionalities in a
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

    callable : bool, optional
        Whether the client can receive calls and notifications. If set to
        `False`, then the client can send notifications and calls, but can not
        receive any.
    """

    def __init__(self, name=None, description=None, metadata=None,
                 callable=True):

        self.hub = SAMPWebHubProxy()

        self.client = SAMPWebClient(self.hub, name, description, metadata,
                                    callable)

    def connect(self, pool_size=20, web_port=21012):
        """
        Connect with the current or specified SAMP Hub, start and register the
        client.

        Parameters
        ----------
        pool_size : int, optional
            The number of socket connections opened to communicate with the
            Hub.
        """
        self.hub.connect(pool_size, web_port=web_port)
        self.client.start()
        self.client.register()

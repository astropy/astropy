# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import socket
import sys
import time
import stat
import copy
import uuid
import re
import select
import threading
import traceback
import datetime
from astropy.extern.six import StringIO
from astropy.extern.six.moves import queue

try:
    import ssl
except ImportError:
    SSL_SUPPORT = False
else:
    SSL_SUPPORT = True

try:
    import bsddb
except ImportError:
    BDB_SUPPORT = False
else:
    BDB_SUPPORT = True

from .constants import SAMP_HUB_SINGLE_INSTANCE, SAMP_RESTRICT_GROUP, SAMP_RESTRICT_OWNER
from .constants import SAMP_HUB_MULTIPLE_INSTANCE, SAMP_STATUS_OK
from .constants import _THREAD_STARTED_COUNT, __profile_version__, __release__
from .errors import SAMPHubError, SAMPProxyError
from .utils import internet_on, SAMPLog, SafeTransport, ServerProxyPool, _HubAsClient
from .utils import WebProfileXMLRPCServer, SecureXMLRPCServer
from .utils import ThreadingXMLRPCServer, WebProfilePopupDialogue

# Python 2 / 3 dependent imports ... for now get from utils to avoid code duplication
from .utils import xmlrpc, urlopen
from .utils import urlparse

if BDB_SUPPORT:
    from .utils import BasicAuthXMLRPCServer

if SSL_SUPPORT and BDB_SUPPORT:
    from .utils import BasicAuthSecureXMLRPCServer

__all__ = ['SAMPHubServer', 'main']

__doctest_skip__ = ['.', 'SAMPHubServer.*']


class SAMPHubServer(object):

    """SAMP Hub Server.

      The SSL parameters are usable only if Python has
      been compiled with SSL support and/or U{ssl <http://docs.python.org/dev/library/ssl.html>}
      module is installed (available by default since Python 2.6).

      Parameters
      ----------
      secret : str
          Secret code.

      addr : str
          Listening address (or IP)

      port : int
          Listening port number.

      lockfile : str
          Custom lockfile name.

      timeout : int
          Hub inactivity timeout. If `timeout` > 0 then the Hub automatically
          stops after an inactivity period longer than `timeout` seconds. By default `timeout`
          is set to 0 (Hub never expires).

      client_timeout : int
          Client inactivity timeout. If `client_timeout` > 0 then the
          Hub automatically unregisters the clients which result inactive for a period longer
          than `client_timeout` seconds. By default `client_timeout` is set to 0 (clients never
          expire).

      log : `SAMPLog`
          A `SAMPLog` instance for the Hub logging. If missing then
          a standard `SAMPLog` instance is used.

      mode : str
          Defines the Hub running mode. If `mode` is 'single' then the Hub runs
          using the standard `.samp` lock-file, having a single instance for user desktop
          session. Otherwise, if `mode` is 'multiple', then the Hub runs using a non-standard
          lock-file, placed in `.samp-1` directory, of the form `samp-hub-<PID>-<ID>`, where
          `<PID>` is the process id and `<ID>` is a general sub-id (integer number).

      label : str
          A string used to label the Hub with a human readable name. This string
          is written in the lock-file assigned to the `hub.label` token.

      owner : str
          General purpose Hub owner name. This value is written in the lock-file
          and assigned to the `hub.owner.name` token.

      owner_group : str
          General purpose Hub owner group name. This value is written in the
          lock-file and assigned to the `hub.owner.group` token.

      auth_file : str
          Authentication file path used for Basic Authentication. The authentication file
          must be a Berkeley DB file in Hash format containing a set of
          `<user name>=md5(<password>)<group 1>,<group 2>,<group 3>,...)` key/value pairs.

      access_restrict : str
          Define whether the Hub access must be restricted to the Hub owner, to a certain owner
          group or not restricted at all.
          Values accepted: `SAMP_RESTRICT_OWNER`, `SAMP_RESTRICT_GROUP`, `None`.

      admin : str
          Define the name of the administrator user in case of restricted access. The administrator user
          can always access the hub instance even if it is running with `SAMP_RESTRICT_OWNER` policy.
          The administrator must be declared in the authentication file.

      https : bool
          Set the Hub running on a Secure Sockets Layer connection (HTTPS)?
          By default SSL is disabled.

      keyfile : str
          Set the file containing the private key for SSL connections. If the
          certificate file (`certfile`) contains the private key, then `keyfile` can be omitted.

      certfile : str
          Specify the file which contains a certificate to be used to identify the
          local side of the secure connection.

      cert_reqs : int
          The parameter `cert_reqs` specifies whether a certificate is required
          from the client side of the connection, and whether it will be validated if provided. It
          must be one of the three values `ssl.CERT_NONE` (certificates ignored), `ssl.CERT_OPTIONAL`
          (not required, but validated if provided), or `ssl.CERT_REQUIRED` (required and validated).
          If the value of this parameter is not `ssl.CERT_NONE`, then the `ca_certs` parameter must
          point to a file of CA certificates.

      ca_certs : str
          The `ca_certs` file contains a set of concatenated "Certification Authority"
          certificates, which are used to validate certificates passed from the client end of the
          connection.

      ssl_version : int
          The `ssl_version` option specifies which version of the SSL protocol to use.
          Typically, the server chooses    a particular protocol version, and the client must adapt to the
          server's choice. Most of the versions are    not interoperable with the other versions. If not
          specified the default SSL version is  `ssl.PROTOCOL_SSLv23`. This version provides the most
          compatibility with other versions client side. Other SSL protocol versions are:
          `ssl.PROTOCOL_SSLv2`, `ssl.PROTOCOL_SSLv3` and `ssl.PROTOCOL_TLSv1`.

      web_profile : bool
          The `web_profile` option enables/disables the Web Profile support.

      pool_size : int
          The number of socket connections opened to communicate with the clients.
    """

    def __init__(self, secret=None, addr=None, port=0, lockfile=None, timeout=0,
                 client_timeout=0, log=None,
                 mode=SAMP_HUB_SINGLE_INSTANCE, label="",
                 owner="", owner_group="", auth_file=None,
                 access_restrict=None, admin="admin", https=False,
                 keyfile=None, certfile=None,
                 cert_reqs=0, ca_certs=None, ssl_version=2, web_profile=True,
                 pool_size=20):
    # General settings
        self._is_running = False
        self._lockfilename = lockfile
        self._admin = admin
        self._addr = addr
        self._port = port
        self._mode = mode
        self._label = label
        self._owner = owner
        self._owner_group = owner_group
        self._timeout = timeout
        self._client_timeout = client_timeout
        if log == None:
            self._log = SAMPLog()
        else:
            self._log = log
        self._pool_size = pool_size

        self._web_profile = web_profile
        self._web_profile_server = None
        self._web_profile_callbacks = {}
        self._web_profile_popup_dialogue = None
        self._web_profile_requests_queue = queue.Queue(1)
        self._web_profile_requests_result = queue.Queue(1)
        self._web_profile_requests_semaphore = queue.Queue(1)
        if web_profile:
            try:
                self._web_profile_server = WebProfileXMLRPCServer(('localhost', 21012), self._log,
                                                                  logRequests=False, allow_none=True)
                self._web_profile_server.register_introspection_functions()
                self._log.info("Hub set to run with Web Profile support enabled.")
            except:
                self._log.error("Port 21012 already in use. Impossible to run the Hub with Web Profile support.")
                self._web_profile = web_profile = False

        # SSL general settings
        self._https = https
        self._keyfile = keyfile
        self._certfile = certfile
        self._cert_reqs = cert_reqs
        self._ca_certs = cert_reqs
        self._ssl_version = ssl_version
        # Basic Authentication settings
        self._auth_file = auth_file
        self._access_restrict = access_restrict

        # Reformat access_restrict string to suitable dictionary
        if access_restrict != None:
            if access_restrict == SAMP_RESTRICT_GROUP:
                access_restrict = {"group": owner_group, "admin": admin}
            elif access_restrict == SAMP_RESTRICT_OWNER:
                access_restrict = {"user": owner, "admin": admin}
            else:
                access_restrict = None

        # Athentication file test
        if auth_file != None:
            if not os.path.isfile(auth_file):
                self._log.error("Unable to load authentication file!")
                raise SAMPHubError("Unable to load authentication file!")

        self._host_name = "127.0.0.1"
        if internet_on():
            try:
                self._host_name = socket.getfqdn()
                socket.getaddrinfo(self._addr or self._host_name, self._port or 0)
            except:
                pass

        # XML-RPC server settings
        if https:

            if keyfile != None and not os.path.isfile(keyfile):
                self._log.error("Unable to load SSL private key file!")
                raise SAMPHubError("Unable to load SSL private key file!")

            if certfile == None or not os.path.isfile(certfile):
                self._log.error("Unable to load SSL cert file!")
                raise SAMPHubError("Unable to load SSL cert file!")

            if auth_file != None:
                self._log.info("Hub set for Basic Authentication using SSL.")
                self._server = BasicAuthSecureXMLRPCServer((self._addr or self._host_name, self._port or 0),
                                                           keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                                                           auth_file, access_restrict, self._log,
                                                           logRequests=False, allow_none=True)
            else:
                self._log.info("Hub set for using SSL.")
                self._server = SecureXMLRPCServer((self._addr or self._host_name, self._port or 0),
                                                  keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                                                  self._log, logRequests=False, allow_none=True)

            self._port = self._server.socket.getsockname()[1]
            self._url = "https://%s:%s" %(self._addr or self._host_name,
                                          self._port)
        else:

            if auth_file != None:
                self._log.info("Hub set for Basic Authentication.")
                self._server = BasicAuthXMLRPCServer((self._addr or self._host_name, self._port or 0),
                                                     auth_file, access_restrict, self._log,
                                                     logRequests=False, allow_none=True)
            else:
                self._server = ThreadingXMLRPCServer((self._addr or self._host_name, self._port or 0),
                                                     self._log, logRequests=False, allow_none=True)

            self._port = self._server.socket.getsockname()[1]
            self._url = "http://%s:%s" %(self._addr or self._host_name,
                                         self._port)

        self._server.register_introspection_functions()

        # Threading stuff
        self._thread_lock = threading.Lock()
        self._thread_run = None
        self._thread_hub_timeout = None
        self._thread_client_timeout = None

        # Variables for timeout testing:
        self._last_activity_time = None
        self._client_activity_time = {}

        # Hub message id counter, used to create hub msg ids
        self._hub_msg_id_counter = 0
        # Hub secred code
        self._hub_secret_code_customized = secret
        self._hub_secret = self._createSecretCode()
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
        self._xmlrpcEndpoints = {}
        # Synchronous message id heap
        self._sync_msg_ids_heap = {}
        # Public ids counter
        self._client_id_counter = -1

        # Standard Profile only operations
        self._server.register_function(self._ping, 'samp.hub.ping')
        self._server.register_function(self._setXmlrpcCallback, 'samp.hub.setXmlrpcCallback')
        # Stadard API operations
        self._server.register_function(self._register, 'samp.hub.register')
        self._server.register_function(self._unregister, 'samp.hub.unregister')
        self._server.register_function(self._declareMetadata, 'samp.hub.declareMetadata')
        self._server.register_function(self._getMetadata, 'samp.hub.getMetadata')
        self._server.register_function(self._declareSubscriptions, 'samp.hub.declareSubscriptions')
        self._server.register_function(self._getSubscriptions, 'samp.hub.getSubscriptions')
        self._server.register_function(self._getRegisteredClients, 'samp.hub.getRegisteredClients')
        self._server.register_function(self._getSubscribedClients, 'samp.hub.getSubscribedClients')
        self._server.register_function(self._notify, 'samp.hub.notify')
        self._server.register_function(self._notifyAll, 'samp.hub.notifyAll')
        self._server.register_function(self._call, 'samp.hub.call')
        self._server.register_function(self._callAll, 'samp.hub.callAll')
        self._server.register_function(self._callAndWait, 'samp.hub.callAndWait')
        self._server.register_function(self._reply, 'samp.hub.reply')
        # Hub as client operations (see _hubAsClientRequestHandler)
        # self._server.register_function(self._receiveNotification, 'samp.client.receiveNotification')
        # self._server.register_function(self._receiveCall, 'samp.client.receiveCall')
        # self._server.register_function(self._receiveResponse, 'samp.client.receiveResponse')

        if web_profile:
            # Web Profile methods like Standard Profile
            self._web_profile_server.register_function(self._ping, 'samp.webhub.ping')
            self._web_profile_server.register_function(self._unregister, 'samp.webhub.unregister')
            self._web_profile_server.register_function(self._declareMetadata, 'samp.webhub.declareMetadata')
            self._web_profile_server.register_function(self._getMetadata, 'samp.webhub.getMetadata')
            self._web_profile_server.register_function(self._declareSubscriptions, 'samp.webhub.declareSubscriptions')
            self._web_profile_server.register_function(self._getSubscriptions, 'samp.webhub.getSubscriptions')
            self._web_profile_server.register_function(self._getRegisteredClients, 'samp.webhub.getRegisteredClients')
            self._web_profile_server.register_function(self._getSubscribedClients, 'samp.webhub.getSubscribedClients')
            self._web_profile_server.register_function(self._notify, 'samp.webhub.notify')
            self._web_profile_server.register_function(self._notifyAll, 'samp.webhub.notifyAll')
            self._web_profile_server.register_function(self._call, 'samp.webhub.call')
            self._web_profile_server.register_function(self._callAll, 'samp.webhub.callAll')
            self._web_profile_server.register_function(self._callAndWait, 'samp.webhub.callAndWait')
            self._web_profile_server.register_function(self._reply, 'samp.webhub.reply')
            # Methods peculiar for Web Profile
            self._web_profile_server.register_function(self._web_profile_register, 'samp.webhub.register')
            self._web_profile_server.register_function(self._web_profile_allowReverseCallbacks, 'samp.webhub.allowReverseCallbacks')
            self._web_profile_server.register_function(self._web_profile_pullCallbacks, 'samp.webhub.pullCallbacks')

    def __del__(self):
        self.stop()

    def _timeoutTestHub(self):

        while self._is_running:
            time.sleep(1)
            self._thread_lock.acquire()
            if self._timeout > 0 and self._last_activity_time != None:
                if time.time() - self._last_activity_time >= self._timeout:
                    self._thread_lock.release()
                    self._log.warning("Timeout expired, Hub is shutting down!")
                    self.stop()
                    break
            if self._thread_lock.locked() == True:
                self._thread_lock.release()

    def _timeoutTestClient(self):

        while self._is_running:
            time.sleep(1)
            if self._client_timeout > 0:
                now = time.time()
                for private_key in self._client_activity_time.keys():
                    if now - self._client_activity_time[private_key] > self._client_timeout \
                       and private_key != self._hub_private_key:
                        self._log.warning("Client %s timeout expired!" % private_key)
                        self._notifyDisconnection(private_key)
                        self._unregister(private_key)

    def _hubAsClientRequestHandler(self, method, args):
        if method == 'samp.client.receiveCall':
            return self._receiveCall(*args)
        elif method == 'samp.client.receiveNotification':
            return self._receiveNotification(*args)
        elif method == 'samp.client.receiveResponse':
            return self._receiveResponse(*args)
        elif method == 'samp.app.ping':
            return self._ping(*args)
        else:
            return _hubAsClientRequestHandler

    def _setupHubAsClient(self):
        result = self._register(self._hub_secret)
        self._hub_public_id = result["samp.self-id"]
        self._hub_private_key = result["samp.private-key"]
        self._setXmlrpcCallback(self._hub_private_key, self._url)
        self._declareMetadata(self._hub_private_key, {"samp.name": "Hub",
                                                      "samp.description.text": self._label,
                                                      "author.name": "Luigi Paioro",
                                                      "author.email": "luigi@iasf-milano.inaf.it",
                                                      "author.affiliation": "INAF-IASF Milano",
                                                      "samp.documentation.url": "http://packages.python.org/sampy/",
                                                      "samp.icon.url": self._url + "/sampy/icon"})
        self._declareSubscriptions(self._hub_private_key, {"samp.app.ping":{}, "x-samp.query.by-meta":{}})

    def start(self, wait=False):
        """Start the current SAMP Hub instance and create the lock file. Hub start-up can
        be blocking or non blocking depending on the `wait` parameter.

        Parameters
        ----------
        wait : bool
            If `True` then the Hub process is joined with the caller, blocking the
            code flow. Usually `True` option is used to run a stand-alone Hub in
            an executable script. If `False` (default), then the Hub process runs in a
            separated thread. `False` is usually used in a Python shell.
        """
        if self._is_running == False:

            self._is_running = True
            self._updateLastActivityTime()

            if self._createLockFile() == False:
                self._is_running = False
                return

            self._setupHubAsClient()
            self._startThreads()

            self._log.info("Hub started")

        if wait and self._is_running:
            self._thread_run.join()

    def _createLockFile(self):

        # Remove lock-files of dead hubs
        self._removeGarbageLockFiles()

        lockfilename = ""
        lockfiledir = ""

        # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
        if "SAMP_HUB" in os.environ:
            # For the time being I assume just the std profile supported.
            if os.environ["SAMP_HUB"].startswith("std-lockurl:"):

                lockfilename = os.environ["SAMP_HUB"][len("std-lockurl:"):]
                lockfile_parsed = urlparse.urlparse(lockfilename)

                if lockfile_parsed[0] != 'file':
                    self._log.warning("Unable to start a Hub with lockfile %s. Start-up process aborted." % lockfilename)
                    return False
                else:
                    lockfilename = lockfile_parsed[2]
        else:
            # If it is a fresh Hub instance
            if self._lockfilename is None:

                self._log.debug("Running mode: " + self._mode)

                if self._mode == SAMP_HUB_SINGLE_INSTANCE:
                    lockfilename = ".samp"
                else:
                    lockfilename = "samp-hub-%d-%s" % (os.getpid(), threading._counter + 1)

                if "HOME" in os.environ:
                    # UNIX
                    lockfiledir = os.environ["HOME"]
                else:
                    # Windows
                    lockfiledir = os.environ["USERPROFILE"]

                if self._mode == SAMP_HUB_MULTIPLE_INSTANCE:
                    lockfiledir = os.path.join(lockfiledir, ".samp-1")

                # If missing create .samp-1 directory
                if not os.path.isdir(lockfiledir):
                    os.mkdir(lockfiledir)
                    os.chmod(lockfiledir, stat.S_IREAD + stat.S_IWRITE + stat.S_IEXEC)

                lockfilename = os.path.join(lockfiledir, lockfilename)

            else:
                self._log.debug("Running mode: multiple")
                lockfilename = self._lockfilename

        hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)

        if hub_is_running:
            self._log.warning("Another SAMP Hub is already running. Start-up process aborted.")
            return False

        self._log.debug("Lock-file: " + lockfilename)

        result = self._new_lockfile(lockfilename)
        if result:
            self._lockfilename = lockfilename

        return result

    def _new_lockfile(self, lockfilename):

        lockfile = open(lockfilename, "w")
        lockfile.close()
        os.chmod(lockfilename, stat.S_IREAD + stat.S_IWRITE)
        lockfile = open(lockfilename, "w")
        lockfile.write("# SAMP lockfile written on %s\n" % datetime.datetime.now().isoformat())
        lockfile.write("# Standard Profile required keys\n")
        lockfile.write("samp.secret=%s\n" % self._hub_secret)
        lockfile.write("samp.hub.xmlrpc.url=%s\n" % self._url)
        lockfile.write("samp.profile.version=%s\n" % __profile_version__)

        # Custom tokens

        lockfile.write("hub.id=%d-%s\n" % (os.getpid(), _THREAD_STARTED_COUNT))

        if self._label == "":
            self._label = "Hub %d-%s" % (os.getpid(), _THREAD_STARTED_COUNT)
        if self._label != "":
            lockfile.write("hub.label=%s\n" % self._label)
        if self._owner != "":
            lockfile.write("hub.owner.name=%s\n" % self._owner)
        if self._owner_group != "":
            lockfile.write("hub.owner.group=%s\n" % self._owner_group)

        if self._auth_file != None:
            lockfile.write("hub.access.auth.file=%s\n" % self._auth_file)

        if self._access_restrict != None:
            lockfile.write("hub.access.auth.restrict=%s\n" % self._access_restrict)

        if SSL_SUPPORT and self._https:
            # Certificate request
            cert_reqs_types = ["NONE", "OPTIONAL", "REQUIRED"]
            lockfile.write("hub.ssl.certificate=%s\n" % cert_reqs_types[self._cert_reqs])
            # SSL protocol version
            ssl_protocol_types = ["SSLv2", "SSLv3", "SSLv23", "TLSv1"]
            lockfile.write("hub.ssl.protocol=%s\n" % ssl_protocol_types[self._ssl_version])

        lockfile.close()

        return True

    def _removeGarbageLockFiles(self):

        lockfilename = ""

        # HUB SINGLE INSTANCE MODE

        if "HOME" in os.environ:
            # UNIX
            lockfilename = os.path.join(os.environ["HOME"], ".samp")
        else:
            # Windows
            lockfilename = os.path.join(os.environ["USERPROFILE"], ".samp")

        hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)

        if not hub_is_running:
            # If lockfilename belongs to a dead hub, then it is deleted
            if os.path.isfile(lockfilename):
                try:
                    os.remove(lockfilename)
                except:
                    pass

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
                    if not hub_is_running:
                        # If lockfilename belongs to a dead hub, then it is deleted
                        if os.path.isfile(lockfilename):
                            try:
                                os.remove(lockfilename)
                            except:
                                pass

    def _startThreads(self):
        self._thread_run = threading.Thread(target=self._serve_forever)
        self._thread_run.setDaemon(True)
        self._thread_hub_timeout = threading.Thread(target=self._timeoutTestHub,
                                                    name="Hub timeout test")
        self._thread_client_timeout = threading.Thread(target=self._timeoutTestClient,
                                                       name="Client timeout test")
        self._thread_run.start()
        self._thread_hub_timeout.start()
        self._thread_client_timeout.start()

    @staticmethod
    def checkRunningHub(lockfilename):
        """Test whether a Hub identified by `lockfilename` is running or not.

        Parameters
        ----------
        lockfilename : str
            Lock-file name (path + file name) of the Hub to be tested.
        """

        is_running = False
        lockfiledict = {}

        # Check whether a lockfile alredy exists
        try:
            if not (lockfilename.startswith("file:") or
                    lockfilename.startswith("http:") or
                    lockfilename.startswith("https:")):
                lockfilename = "file://" + lockfilename

            lockfile = urlopen(lockfilename)
            lockfile_content = lockfile.readlines()
            lockfile.close()
        except:
            return is_running, lockfiledict

        for line in lockfile_content:
            if not line.startswith(b"#"):
                kw, val = line.split(b"=")
                lockfiledict[kw.decode().strip()] = val.decode().strip()

        if "samp.hub.xmlrpc.url" in lockfiledict:
            try:
                proxy = xmlrpc.ServerProxy(lockfiledict["samp.hub.xmlrpc.url"].replace("\\", ""),
                                           allow_none=1)
                proxy.samp.hub.ping()
                is_running = True
            except xmlrpc.ProtocolError:
                # There is a protocol error (e.g. for authentication required),
                # but the server is alive
                is_running = True
            except:
                if SSL_SUPPORT:
                    if sys.exc_info()[0] == ssl.SSLError:
                        # SSL connection refused for certifcate reasons...
                        # anyway the server is alive
                        is_running = True

        return is_running, lockfiledict

    def _createSecretCode(self):
        if self._hub_secret_code_customized != None:
            return self._hub_secret_code_customized
        else:
            return str(uuid.uuid1())

    def stop(self):
        """Stop the current SAMP Hub instance and delete the lock file."""
        if self._is_running:

            self._log.info("Hub is stopping...")

            self._notifyShutdown()

            self._is_running = False

            if (os.path.isfile(self._lockfilename)):
                lockfile = open(self._lockfilename, "r")
                lockfile_content = lockfile.readlines()
                lockfile.close()
                for line in lockfile_content:
                    if line.strip()[0] != "#":
                        kw, val = line.split("=")
                        if kw.strip() == "samp.secret" and val.strip() == self._hub_secret:
                            os.remove(self._lockfilename)
                            break

        # Reset vaiables
        self._joinAllThreads()

        self._hub_msg_id_counter = 0
        self._hub_secret = self._createSecretCode()
        self._hub_public_id = ""
        self._metadata = {}
        self._private_keys = {}
        self._mtype2ids = {}
        self._id2mtypes = {}
        self._xmlrpcEndpoints = {}
        self._last_activity_time = None
        self._lockfilename = None

        self._log.info("Hub stopped.")

    def _joinOneThread(self, thread_name, timeout):
        t = getattr(self, thread_name)
        if t is None:
            return
        t.join(timeout)
        setattr(self, thread_name, None)

    def _joinAllThreads(self, timeout=1):
        for thread_name in [
            "_thread_run",
            "_thread_hub_timeout",
                "_thread_client_timeout"]:
            self._joinOneThread(thread_name, timeout)

    def isRunning(self):
        """Return an information concerning the Hub running status.

        Returns
        -------
        running : bool
            Is the hub running?
        """
        return self._is_running

    def _serve_forever(self):

        if self._web_profile:
            self._web_profile_popup_dialogue = \
                WebProfilePopupDialogue(self._web_profile_requests_result)

        while self._is_running:

            try:
                r = w = e = None
                r, w, e = select.select([self._server.socket], [], [], 0.1)
            except:
                pass
            if r:
                self._server.handle_request()

            if self._web_profile:

                try:
                    request = self._web_profile_requests_queue.get_nowait()
                    self._web_profile_popup_dialogue.showPopup(request)
                except queue.Empty:
                    pass

                try:
                    r = w = e = None
                    r, w, e = select.select([self._web_profile_server.socket], [], [], 0.01)
                    self._web_profile_popup_dialogue.update()
                except:
                    pass
                if r:
                    self._web_profile_server.handle_request()

    def _notifyShutdown(self):
        msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.shutdown")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                for key in self._mtype2ids[mtype]:
                    self._notify_(self._hub_private_key, self._private_keys[key][0],
                                  {"samp.mtype":"samp.hub.event.shutdown",
                                   "samp.params": {}})

    def _notifyRegister(self, private_key):
        msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.register")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                public_id = self._private_keys[private_key][0]
                for key in self._mtype2ids[mtype]:
                    # if key != private_key:
                    self._notify(self._hub_private_key, self._private_keys[key][0],
                                 {"samp.mtype":"samp.hub.event.register",
                                  "samp.params": {"id": public_id}})

    def _notifyUnregister(self, private_key):
        msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.unregister")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                public_id = self._private_keys[private_key][0]
                for key in self._mtype2ids[mtype]:
                    if key != private_key:
                        self._notify(self._hub_private_key, self._private_keys[key][0],
                                     {"samp.mtype":"samp.hub.event.unregister",
                                      "samp.params": {"id": public_id}})

    def _notifyMetadata(self, private_key):
        msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.metadata")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                public_id = self._private_keys[private_key][0]
                for key in self._mtype2ids[mtype]:
                    # if key != private_key:
                    self._notify(self._hub_private_key, self._private_keys[key][0],
                                 {"samp.mtype":"samp.hub.event.metadata",
                                  "samp.params": {"id": public_id,
                                                  "metadata": self._metadata[private_key]}
                                  })

    def _notifySubscriptions(self, private_key):
        msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.subscriptions")
        for mtype in msubs:
            if mtype in self._mtype2ids:
                public_id = self._private_keys[private_key][0]
                for key in self._mtype2ids[mtype]:
                    # if key != private_key:
                    self._notify(self._hub_private_key, self._private_keys[key][0],
                                 {"samp.mtype":"samp.hub.event.subscriptions",
                                  "samp.params": {"id": public_id,
                                                  "subscriptions": self._id2mtypes[private_key]}
                                  })

    def _notifyDisconnection(self, private_key):

        def _xmlrpc_call_disconnect(endpoint, private_key, hub_public_id, message):
            try:
                endpoint.samp.client.receiveNotification(private_key, hub_public_id, message)
            except:
                pass

        msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.disconnect")
        public_id = self._private_keys[private_key][0]
        endpoint = self._xmlrpcEndpoints[public_id][1]

        for mtype in msubs:
            if mtype in self._mtype2ids and private_key in self._mtype2ids[mtype]:
                try:
                    self._log.debug("notify disconnection to %s" % (public_id))
                    threading.Thread(target=_xmlrpc_call_disconnect,
                                     args=(endpoint, private_key, self._hub_public_id,
                                    {"samp.mtype":"samp.hub.disconnect",
                                    "samp.params": {"reason": "Timeout expired!"}})).start()
                except:
                    self._log.warning("disconnection notification to client %s failed\n" % (public_id))

    def _ping(self):
        self._updateLastActivityTime()
        self._log.debug("ping")
        return "1"

    def _query_by_metadata(self, key, value):
        public_id_list = []
        for private_id in self._metadata:
            if key in self._metadata[private_id]:
                if self._metadata[private_id][key] == value:
                    public_id_list.append(self._private_keys[private_id][0])

        return public_id_list

    def _setXmlrpcCallback(self, private_key, xmlrpc_addr):
        self._updateLastActivityTime(private_key)
        if private_key in self._private_keys:
            if private_key == self._hub_private_key:
                self._xmlrpcEndpoints[self._private_keys[private_key][0]] = (xmlrpc_addr, _HubAsClient(self._hubAsClientRequestHandler))
                return ""
            # Dictionary stored with the public id
            self._log.debug("setXmlrpcCallback: %s %s" % (private_key, xmlrpc_addr))
            server_proxy_pool = None
            if SSL_SUPPORT and xmlrpc_addr[0:5] == "https":
                server_proxy_pool = ServerProxyPool(self._pool_size, xmlrpc.ServerProxy,
                                                    xmlrpc_addr, transport=SafeTransport(key_file=self._keyfile,
                                                   cert_file=self._certfile,
                                                   cert_reqs=self._cert_reqs,
                                                   ca_certs=self._ca_certs,
                                                   ssl_version=ssl.PROTOCOL_SSLv3),
                                                    allow_none=1)
            else:
                server_proxy_pool = ServerProxyPool(self._pool_size, xmlrpc.ServerProxy,
                                                    xmlrpc_addr, allow_none=1)

            self._xmlrpcEndpoints[self._private_keys[private_key][0]] = (xmlrpc_addr, server_proxy_pool)
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

        return ""

    def _perform_standard_register(self):

        self._thread_lock.acquire()
        private_key, public_id = self._getNewIds()
        self._thread_lock.release()
        self._private_keys[private_key] = (public_id, time.time())
        self._updateLastActivityTime(private_key)
        self._notifyRegister(private_key)
        self._log.debug("register: private-key = %s and self-id = %s" % (private_key, public_id))
        return {"samp.self-id": public_id,
                "samp.private-key": private_key,
                "samp.hub-id": self._hub_public_id}

    def _register(self, secret):
        self._updateLastActivityTime()
        if secret == self._hub_secret:
            return self._perform_standard_register()
        else:
            # return {"samp.self-id": "", "samp.private-key": "", "samp.hub-id": ""}
            raise SAMPProxyError(7, "Bad secret code")

    def _getNewIds(self):
        private_key = str(uuid.uuid1())
        self._client_id_counter += 1
        public_id = 'cli#hub'
        if self._client_id_counter > 0:
            public_id = "cli#%d" % (self._client_id_counter)

        return private_key, public_id

    def _unregister(self, private_key):

        self._updateLastActivityTime()

        public_key = ""

        self._notifyUnregister(private_key)

        self._thread_lock.acquire()

        if private_key in self._private_keys:
            public_key = self._private_keys[private_key][0]
            del self._private_keys[private_key]
        else:
            self._thread_lock.release()
            return ""

        if private_key in self._metadata:
            del self._metadata[private_key]

        if private_key in self._id2mtypes:
            del self._id2mtypes[private_key]

        for mtype in self._mtype2ids.keys():
            if private_key in self._mtype2ids[mtype]:
                self._mtype2ids[mtype].remove(private_key)

        if public_key in self._xmlrpcEndpoints:
            del self._xmlrpcEndpoints[public_key]

        if private_key in self._client_activity_time:
            del self._client_activity_time[private_key]

        if self._web_profile:
            if private_key in self._web_profile_callbacks:
                del self._web_profile_callbacks[private_key]
            self._web_profile_server.remove_client(private_key)

        self._thread_lock.release()

        self._log.debug("unregister %s (%s)" % (public_key, private_key))

        return ""

    def _declareMetadata(self, private_key, metadata):
        self._updateLastActivityTime(private_key)
        if private_key in self._private_keys:
            self._log.debug("declareMetadata: private-key = %s metadata = %s" % (private_key, str(metadata)))
            self._metadata[private_key] = metadata
            self._notifyMetadata(private_key)
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)
        return ""

    def _getMetadata(self, private_key, client_id):
        self._updateLastActivityTime(private_key)
        if private_key in self._private_keys:
            client_private_key = self._getPrivateKeyFromPublicId(client_id)
            self._log.debug("getMetadata: private-key = %s client-id = %s" %
                            (private_key, client_id))
            if client_private_key != None:
                if client_private_key in self._metadata:
                    self._log.debug("--> metadata = %s" % self._metadata[client_private_key])
                    return self._metadata[client_private_key]
                else:
                    return {}
            else:
                self._log.warning("Invalid client ID %s" % client_id)
                raise SAMPProxyError(6, "Invalid client ID")
        else:
            self._log.warning("Private-key %s expired or invalid." % private_key)
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    def _declareSubscriptions(self, private_key, mtypes):

        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:

            self._log.debug("declareSubscriptions: private-key = %s mtypes = %s" % (private_key, str(mtypes)))

            # remove subscription to previous mtypes
            if private_key in self._id2mtypes:

                prev_mtypes = self._id2mtypes[private_key]

                for mtype in prev_mtypes:
                    try:
                        self._mtype2ids[mtype].remove(private_key)
                    except:
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

            self._log.debug("declareSubscriptions: subscriptions accepted from %s => %s" % (private_key, str(mtypes)))

            for mtype in mtypes:

                if mtype in self._mtype2ids:
                    if not private_key in self._mtype2ids[mtype]:
                        self._mtype2ids[mtype].append(private_key)
                else:
                    self._mtype2ids[mtype] = [private_key]

            self._notifySubscriptions(private_key)

        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

        return ""

    def _getSubscriptions(self, private_key, client_id):

        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:
            client_private_key = self._getPrivateKeyFromPublicId(client_id)
            if client_private_key != None:
                if client_private_key in self._id2mtypes:
                    self._log.debug("getSubscriptions: client-id = %s mtypes = %s" %
                                    (client_id, str(self._id2mtypes[client_private_key])))
                    return self._id2mtypes[client_private_key]
                else:
                    self._log.debug("getSubscriptions: client-id = %s mtypes = missing" % client_id)
                    return {}
            else:
                raise SAMPProxyError(6, "Invalid client ID")
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    def _getRegisteredClients(self, private_key):

        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:
            reg_clients = []
            for pkey in self._private_keys.keys():
                if pkey != private_key:
                    reg_clients.append(self._private_keys[pkey][0])
            self._log.debug("getRegisteredClients: private_key = %s clients = %s" % (private_key, reg_clients))
            return reg_clients
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    def _getSubscribedClients(self, private_key, mtype):

        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:
            sub_clients = {}

            for pkey in self._private_keys.keys():
                if pkey != private_key and self._isSubscribed(pkey, mtype):
                    sub_clients[self._private_keys[pkey][0]] = {}

            self._log.debug("getSubscribedClients: private_key = %s mtype = %s clients = %s" %
                            (private_key, mtype, sub_clients))
            return sub_clients
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    @staticmethod
    def getMTypeSubtypes(mtype):
        """Return a list containing all the possible wildcarded subtypes of MType.

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
        >>> import astropy.vo.samp as sampy
        >>> sampy.SAMPHubServer.getMTypeSubtypes("samp.app.ping")
        ['samp.app.ping', 'samp.app.*', 'samp.*', '*']
        """

        subtypes = []

        msubs = mtype.split(".")
        indexes = list(range(len(msubs)))
        indexes.reverse()
        indexes.append(-1)

        for i in indexes:
            tmp_mtype = ".".join(msubs[:i+1])
            if tmp_mtype != mtype:
                if tmp_mtype != "":
                    tmp_mtype = tmp_mtype + ".*"
                else:
                    tmp_mtype = "*"
            subtypes.append(tmp_mtype)

        return subtypes

    def _isSubscribed(self, private_key, mtype):

        subscribed = False

        msubs = SAMPHubServer.getMTypeSubtypes(mtype)

        for msub in msubs:
            if msub in self._mtype2ids:
                if private_key in self._mtype2ids[msub]:
                    subscribed = True

        return subscribed

    def _notify(self, private_key, recipient_id, message):
        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:
            if self._isSubscribed(self._getPrivateKeyFromPublicId(recipient_id),
                                  message["samp.mtype"]) == False:
                raise SAMPProxyError(2, "Client %s not subscribed to MType %s" % (recipient_id, message["samp.mtype"]))

            threading.Thread(target=self._notify_, args=(private_key, recipient_id, message)).start()
            return {}
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    def _notify_(self, sender_private_key, recipient_public_id, message):
        if sender_private_key in self._private_keys:
            sender_public_id = self._private_keys[sender_private_key][0]
            try:
                self._log.debug("notify %s from %s to %s" % (message["samp.mtype"], sender_public_id, recipient_public_id))
                recipient_private_key = self._getPrivateKeyFromPublicId(recipient_public_id)
                if recipient_private_key != None:
                    for attempt in range(10):
                        if self._is_running:
                            try:
                                if self._web_profile and recipient_private_key in self._web_profile_callbacks:
                                    # Web Profile
                                    self._web_profile_callbacks[recipient_private_key].put({"samp.methodName": "receiveNotification",
                                                                                            "samp.params": [sender_public_id, message]})
                                    return
                                # Standard Profile
                                self._xmlrpcEndpoints[recipient_public_id][1].samp.client.receiveNotification(
                                    recipient_private_key, sender_public_id, message)
                                break
                            except:
                                err = StringIO()
                                traceback.print_exc(file=err)
                                txt = err.getvalue()
                                self._log.debug("%s XML-RPC endpoint error (attempt %d): \n%s" % (recipient_public_id, attempt + 1, txt))
                                time.sleep(0.01)
                else:
                    raise SAMPProxyError(6, "Invalid client ID")
            except:
                self._log.warning("%s notification from client %s to client %s failed\n" % (message["samp.mtype"], sender_public_id, recipient_public_id))

    def _notifyAll(self, private_key, message):
        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:
            if not "samp.mtype" in message:
                raise SAMPProxyError(3, "samp.mtype keyword is missing")
            recipient_ids = self._notifyAll(private_key, message)
            return recipient_ids
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    def _notifyAll(self, sender_private_key, message):

        recipient_ids = []
        msubs = SAMPHubServer.getMTypeSubtypes(message["samp.mtype"])

        for mtype in msubs:
            if mtype in self._mtype2ids:
                for key in self._mtype2ids[mtype]:
                    if key != sender_private_key:
                        _recipient_id = self._private_keys[key][0]
                        recipient_ids.append(_recipient_id)
                        threading.Thread(target=self._notify,
                                         args=(sender_private_key,
                                               _recipient_id, message)
                                         ).start()

        return recipient_ids

    def _call(self, private_key, recipient_id, msg_tag, message):
        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:
            if self._isSubscribed(self._getPrivateKeyFromPublicId(recipient_id),
                                  message["samp.mtype"]) == False:
                raise SAMPProxyError(2, "Client %s not subscribed to MType %s" % (recipient_id, message["samp.mtype"]))
            public_id = self._private_keys[private_key][0]
            msg_id = self._getNewHubMsgId(public_id, msg_tag)
            threading.Thread(target=self._call_, args=(private_key, public_id, recipient_id, msg_id, message)).start()
            return msg_id
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    def _call_(self, sender_private_key, sender_public_id, recipient_public_id, msg_id, message):
        if sender_private_key in self._private_keys:
            try:
                self._log.debug("call %s from %s to %s (%s)" % (msg_id.split(";;")[0], sender_public_id, recipient_public_id, message["samp.mtype"]))
                recipient_private_key = self._getPrivateKeyFromPublicId(recipient_public_id)
                if recipient_private_key != None:
                    for attempt in range(10):
                        if self._is_running:
                            try:
                                if self._web_profile and recipient_private_key in self._web_profile_callbacks:
                                    # Web Profile
                                    self._web_profile_callbacks[recipient_private_key].put({"samp.methodName": "receiveCall",
                                                                                            "samp.params": [sender_public_id, msg_id, message]})
                                    return
                                # Standard Profile
                                self._xmlrpcEndpoints[recipient_public_id][1].samp.client.receiveCall(
                                    recipient_private_key, sender_public_id, msg_id, message)
                                break
                            except:
                                err = StringIO()
                                traceback.print_exc(file=err)
                                txt = err.getvalue()
                                self._log.debug("%s XML-RPC endpoint error (attempt %d): \n%s" % (recipient_public_id, attempt + 1, txt))
                                time.sleep(0.01)
                else:
                    raise SAMPProxyError(6, "Invalid client ID")
            except:
                self._log.warning("%s call %s from client %s to client %s failed\n" % (message["samp.mtype"], msg_id.split(";;")[0], sender_public_id, recipient_public_id))

    def _callAll(self, private_key, msg_tag, message):
        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:
            if not "samp.mtype" in message:
                raise SAMPProxyError(3, "samp.mtype keyword is missing in message tagged as %s" % msg_tag)

            public_id = self._private_keys[private_key][0]
            msg_id = self._callAll_(private_key, public_id, msg_tag, message)
            return msg_id
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    def _callAll_(self, sender_private_key, sender_public_id, msg_tag, message):

        msg_id = {}
        msubs = SAMPHubServer.getMTypeSubtypes(message["samp.mtype"])

        for mtype in msubs:
            if mtype in self._mtype2ids:
                for key in self._mtype2ids[mtype]:
                    if key != sender_private_key:
                        _msg_id = self._getNewHubMsgId(sender_public_id, msg_tag)
                        receiver_public_id = self._private_keys[key][0]
                        msg_id[receiver_public_id] = _msg_id
                        threading.Thread(target=self._call_,
                                         args=(sender_private_key, sender_public_id,
                                               receiver_public_id, _msg_id, message)
                                         ).start()
        return msg_id

    def _callAndWait(self, private_key, recipient_id, message, timeout):
        self._updateLastActivityTime(private_key)

        if private_key in self._private_keys:
            timeout = int(timeout)

            now = time.time()
            response = {}

            msg_id = self._call(private_key, recipient_id, "sampy::sync::call", message)
            self._sync_msg_ids_heap[msg_id] = None

            while self._is_running:
                if timeout > 0 and time.time() - now >= timeout:
                    del(self._sync_msg_ids_heap[msg_id])
                    raise SAMPProxyError(1, "Timeout expired!")

                if self._sync_msg_ids_heap[msg_id] != None:
                    response = copy.deepcopy(self._sync_msg_ids_heap[msg_id])
                    del(self._sync_msg_ids_heap[msg_id])
                    break
                time.sleep(0.01)

            return response
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    def _reply(self, private_key, msg_id, response):
        self._updateLastActivityTime(private_key)
        if private_key in self._private_keys:
            threading.Thread(target=self._reply_, args=(private_key, msg_id, response)).start()
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

        return {}

    def _reply_(self, responder_private_key, msg_id, response):

        if (responder_private_key in self._private_keys) and msg_id:
            responder_public_id = self._private_keys[responder_private_key][0]
            counter, hub_public_id, recipient_public_id, \
                recipient_msg_tag = msg_id.split(";;", 3)

            try:
                self._log.debug("reply %s from %s to %s" % (counter, responder_public_id, recipient_public_id))
                if recipient_msg_tag == "sampy::sync::call":
                    if msg_id in self._sync_msg_ids_heap.keys():
                        self._sync_msg_ids_heap[msg_id] = response
                else:
                    recipient_private_key = self._getPrivateKeyFromPublicId(recipient_public_id)
                    if recipient_private_key != None:
                        for attempt in range(10):
                            if self._is_running:
                                try:
                                    if self._web_profile and recipient_private_key in self._web_profile_callbacks:
                                        # Web Profile
                                        self._web_profile_callbacks[recipient_private_key].put({"samp.methodName": "receiveResponse",
                                                                                                "samp.params": [responder_public_id,
                                                                                                                recipient_msg_tag, response]})
                                        return
                                    # Standard Profile
                                    self._xmlrpcEndpoints[recipient_public_id][1].samp.client.receiveResponse(
                                        recipient_private_key, responder_public_id, recipient_msg_tag, response)
                                    break
                                except:
                                    err = StringIO()
                                    traceback.print_exc(file=err)
                                    txt = err.getvalue()
                                    self._log.debug("%s XML-RPC endpoint error (attempt %d): \n%s" % (recipient_public_id, attempt + 1, txt))
                                    time.sleep(0.01)
                    else:
                        raise SAMPProxyError(6, "Invalid client ID")
            except:
                self._log.warning("%s reply from client %s to client %s failed\n" % (recipient_msg_tag, responder_public_id, recipient_public_id))

    def _getPrivateKeyFromPublicId(self, public_id):

        for private_key in self._private_keys.keys():
            if self._private_keys[private_key][0] == public_id:
                return private_key
        return None

    def _getNewHubMsgId(self, sender_public_id, sender_msg_id):
        self._thread_lock.acquire()
        self._hub_msg_id_counter += 1
        self._thread_lock.release()
        return "msg#%d;;%s;;%s;;%s" % \
               (self._hub_msg_id_counter, self._hub_public_id,
                sender_public_id, sender_msg_id)

    def _updateLastActivityTime(self, private_key=None):
        self._thread_lock.acquire()
        self._last_activity_time = time.time()
        if private_key != None:
            self._client_activity_time[private_key] = time.time()
        self._thread_lock.release()

    def _receiveNotification(self, private_key, sender_id, message):
        return ""

    def _receiveCall(self, private_key, sender_id, msg_id, message):
        if private_key == self._hub_private_key:

            if "samp.mtype" in message and message["samp.mtype"] == "samp.app.ping":
                self._reply(self._hub_private_key, msg_id, {"samp.status": SAMP_STATUS_OK, "samp.result": {}})

            elif "samp.mtype" in message and \
                 (message["samp.mtype"] == "x-samp.query.by-meta" or
                  message["samp.mtype"] == "samp.query.by-meta"):

                ids_list = self._query_by_metadata(message["samp.params"]["key"], message["samp.params"]["value"])
                self._reply(self._hub_private_key, msg_id,
                            {"samp.status": SAMP_STATUS_OK,
                             "samp.result": {"ids": ids_list}})

            return ""
        else:
            return ""

    def _receiveResponse(self, private_key, responder_id, msg_tag, response):
        return ""

    def _web_profile_register(self, identity_info, client_address=("uknown", 0), origin = "unknown"):

        self._updateLastActivityTime()

        if not client_address[0] in ["localhost", "127.0.0.1"]:
            raise SAMPProxyError(403, "Request of registration rejected by the Hub.")

        if not origin:
            origin = "unknown"

        if isinstance(identity_info, dict):
            # an old version of the protocol provided just a string with the app name
            if "samp.name" not in identity_info:
                raise SAMPProxyError(403, "Request of registration rejected by the Hub (application name not provided).")

        # Red semaphore for the other threads
        self._web_profile_requests_semaphore.put("wait")
        # Set the request to be displayed for the current thread
        self._web_profile_requests_queue.put((identity_info, client_address, origin))
        # Get the popup dialogue response
        response = self._web_profile_requests_result.get()
        # OK, semaphore green
        self._web_profile_requests_semaphore.get()

        if response:
            register_map = self._perform_standard_register()
            register_map["samp.url-translator"] = "http://localhost:21012/translator/%s?ref=" % register_map["samp.private-key"]
            self._web_profile_server.add_client(register_map["samp.private-key"])
            return register_map
        else:
            raise SAMPProxyError(403, "Request of registration rejected by the user.")

    def _web_profile_allowReverseCallbacks(self, private_key, allow):
        self._updateLastActivityTime()
        if private_key in self._private_keys:
            if allow == "0":
                if private_key in self._web_profile_callbacks:
                    del self._web_profile_callbacks[private_key]
            else:
                self._web_profile_callbacks[private_key] = queue.Queue()
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)
        return ""

    def _web_profile_pullCallbacks(self, private_key, timeout_secs):
        self._updateLastActivityTime()
        if private_key in self._private_keys:
            callback = []
            try:
                while self._is_running:
                    item_queued = self._web_profile_callbacks[private_key].get(block=True, timeout=int(timeout_secs))
                    callback.append(item_queued)
                    if self._web_profile_callbacks[private_key].empty():
                        break
            except queue.Empty:
                pass
            return callback
        else:
            raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)


def main(timeout=0):

    from optparse import OptionParser, OptionGroup
    parser = OptionParser(version="%prog " + __release__)

    parser.disable_interspersed_args()

    parser.add_option("-k", "--secret", dest="secret", metavar="CODE",
                      help="custom secret code.")

    parser.add_option("-d", "--addr", dest="addr", metavar="ADDR",
                      help="listening address (or IP).")

    parser.add_option("-p", "--port", dest="port", metavar="PORT", type="int",
                      help="listening port number.")

    parser.add_option("-f", "--lockfile", dest="lockfile", metavar="FILE",
                      help="custom lockfile.")

    parser.add_option("-w", "--no-web-profile", dest="web_profile", action="store_false",
                      help="run the Hub disabling the Web Profile.")

    parser.add_option("-P", "--pool-size", dest="pool_size", metavar="SIZE", type="int",
                      help="the socket connections pool size.")

    timeout_group = OptionGroup(parser, "Timeout group",
                                "Special options to setup hub and client timeouts."
                                "It contains a set of special options that allows to set up the Hub and "
                                "clients inactivity timeouts, that is the Hub or client inactivity time "
                                "interval after which the Hub shuts down or unregisters the client. "
                                "Notification of samp.hub.disconnect MType is sent to the clients "
                                "forcibly unregistered for timeout expiration.")

    timeout_group.add_option("-t", "--timeout", dest="timeout", metavar="SECONDS",
                             help="set the Hub inactivity timeout in SECONDS. By default it "
                             "is set to 0, that is the Hub never expires.", type="int")

    timeout_group.add_option("-c", "--client-timeout", dest="client_timeout", metavar="SECONDS",
                             help="set the client inactivity timeout in SECONDS. By default it "
                             "is set to 0, that is the client never expires.", type="int")

    parser.add_option_group(timeout_group)

    log_group = OptionGroup(parser, "Logging options",
                            "Additional options which allow to customize the logging output. By "
                            "default SAMPy Hub uses the standard output and standard error "
                            "devices to print out INFO level logging messages. Using the options "
                            "here below it is possible to modify the logging level and also "
                            "specify the output files where redirect the logging messages.")

    log_group.add_option("-L", "--log-level", dest="loglevel", metavar="LEVEL",
                         help="set the Hub instance log level (OFF, ERROR, WARNING, INFO, DEBUG).",
                         type="choice", choices=["OFF", "ERROR", "WARNING", "INFO", "DEBUG"])

    log_group.add_option("-O", "--log-output", dest="logout", metavar="FILE",
                         help="set the output file for INFO and DEBUG messages.")

    log_group.add_option("-E", "--log-error", dest="logerr", metavar="FILE",
                         help="set the output file for WARNING and ERROR messages.")

    parser.add_option_group(log_group)

    adv_group = OptionGroup(parser, "Advanced group",
                            "Advanced options addressed to facilitate administrative tasks and "
                            "allow new non-standard Hub behaviors. In particular the --label "
                            "options is used to assign a value to hub.label token and is used to "
                            "assign a name to the Hub instance. Option --auth-file sets the "
                            "authentication file enabling the access control through Basic "
                            "Authentication. The authentication file is a Berkeley DB file "
                            "in Hash format containing a set of <user name>=md5(<password>)<group 1>,<group 2>,<group 3>,... "
                            "key=value pairs. Options --owner and --group are additional tokens "
                            "(hub.owner.name and hub.owner.group respectively) available for "
                            "possible administrative tasks and usable in conjunction with --restrict "
                            "option for authentication tasks. --restrict option allows to restrict "
                            "the Hub access only to OWNER or to a certain owner GROUP. Access restriction "
                            "is actually enabled only if an authentication file is provided. "
                            "--admin option defines the name of the administrator user in case "
                            "of restricted access. The administrator user can always access the "
                            "hub instance even if it is running with the OWNER restricion policy. "
                            "The administrator must be declared in the authentication file. "
                            "The very special --multi option allows to start a Hub in multi-instance mode. "
                            "Multi-instance mode is a non-standard Hub behavior that enables "
                            "multiple contemporaneous running Hubs. Multi-instance hubs place "
                            "their non-standard lock-files within the <home directory>/.samp-1 "
                            "directory naming them making use of the format: "
                            "samp-hub-<PID>-<ID>, where PID is the Hub process ID while ID is an "
                            "internal ID (integer).")

    adv_group.add_option("-l", "--label", dest="label", metavar="LABEL",
                         help="assign a LABEL to the Hub.")

    adv_group.add_option("-o", "--owner", dest="owner", metavar="OWNER",
                         help="assign a OWNER to the Hub.")

    adv_group.add_option("-g", "--group", dest="owner_group", metavar="GROUP",
                         help="assign a GROUP to the Hub.")

    if BDB_SUPPORT:
        adv_group.add_option("-a", "--auth-file", dest="auth_file", metavar="FILE",
                             help="filename for Basic Authentication.")

        adv_group.add_option("-r", "--restrict", dest="access_restrict", metavar="LEVEL",
                             help="set the access restriction to OWNER or GROUP.",
                             type="choice", choices=["OWNER", "GROUP"])

        adv_group.add_option("-A", "--admin", dest="admin", metavar="NAME",
                             help="set the name of the administrator user.")

    adv_group.add_option("-m", "--multi", dest="mode",
                         help="run the Hub in multi-instance mode generating a custom "
                         "lockfile with a random name.",
                         action="store_const", const=SAMP_HUB_MULTIPLE_INSTANCE)

    parser.add_option_group(adv_group)

    if SSL_SUPPORT:

        ssl_group = OptionGroup(parser, "SSL group", "Additional options to launch "
                                "the Hub instance using the Secure Sockets Layer (HTTPS). "
                                "The --key-file and --cert-file parameters specify optional "
                                "files which contain a certificate to be used to identify "
                                "the local side of the connection. "
                                "Often the private key is stored in the same file as the "
                                "certificate; in this case, only the --cert-file parameter need "
                                "be passed. If the private key is stored in a separate file, "
                                "both parameters must be used. If the private key is stored "
                                "in the certificate file, it should come before the first certificate "
                                "in the certificate chain.")

        ssl_group.add_option("-s", "--https", dest="https", action="store_true",
                             help="run the Hub using the Secure Sockets Layer.")

        ssl_group.add_option("-C", "--cert-file", dest="certfile", metavar="FILE",
                             help="set the certificate file.")

        ssl_group.add_option("-K", "--key-file", dest="keyfile", metavar="FILE",
                             help="set the key file. By default this option is ignored, "
                             "assuming that the private key is stored in the certificate file.")

        ssl_group.add_option("--cert-reqs", dest="cert_reqs", metavar="STRING",
                             help="this option specifies whether a certificate "
                             "is required from the client side of the connection, and whether "
                             "it will be validated if provided. It must be one of the three "
                             "values NONE (certificates ignored, default), OPTIONAL (not "
                             "required, but validated if provided), or REQUIRED (required "
                             "and validated). If the value of this option is not NONE, "
                             "then the --ca-certs option must point to a file of CA certificates.",
                             type="choice", choices=["NONE", "OPTIONAL", "REQUIRED"])

        ssl_group.add_option("--ca-certs", dest="ca_certs", metavar="FILE",
                             help="the --ca-certs file contains a set of concatenated "
                             "\"certification authority\" certificates, which are used to "
                             "validate certificates passed from the client end of the "
                             "connection.")

        ssl_group.add_option("--ssl-version", dest="ssl_version", metavar="STRING",
                             help="the --ssl-version option specifies which version of the "
                             "SSL protocol to use. Typically, the server chooses a particular "
                             "protocol version, and the client must adapt to the server's choice. "
                             "Most of the versions are not interoperable with the other versions. "
                             "If not specified the default SSL version is SSLv23. This version "
                             "provides the most compatibility with other versions client side. "
                             "Other SSL protocol versions are: SSLv2, SSLv3 and TLSv1.",
                             type="choice", choices=["SSLv23", "SSLv2", "SSLv3", "TLSv1"])

        parser.add_option_group(ssl_group)

        parser.set_defaults(https=False)
        parser.set_defaults(certfile=None)
        parser.set_defaults(keyfile=None)
        parser.set_defaults(cert_reqs="NONE")
        parser.set_defaults(ca_certs=None)
        parser.set_defaults(ssl_version="SSLv23")

    parser.set_defaults(web_profile=True)
    parser.set_defaults(pool_size=20)
    parser.set_defaults(timeout=0)
    parser.set_defaults(client_timeout=0)
    parser.set_defaults(loglevel="INFO")
    parser.set_defaults(logout="")
    parser.set_defaults(logerr="")
    parser.set_defaults(label="")
    parser.set_defaults(owner="")
    parser.set_defaults(owner_group="")
    parser.set_defaults(admin="admin")

    if BDB_SUPPORT:
        parser.set_defaults(auth_file=None)
        parser.set_defaults(access_restrict=None)

    parser.set_defaults(mode=SAMP_HUB_SINGLE_INSTANCE)

    (options, args) = parser.parse_args()

    try:

        if SSL_SUPPORT:

            # Set ssl options properly

            if options.cert_reqs == "OPTIONAL":
                options.cert_reqs = ssl.CERT_OPTIONAL
            elif options.cert_reqs == "REQUIRED":
                options.cert_reqs = ssl.CERT_REQUIRED
            else:
                options.cert_reqs = ssl.CERT_NONE

            if options.ssl_version == "SSLv2":
                options.ssl_version = ssl.PROTOCOL_SSLv2
            elif options.ssl_version == "SSLv3":
                options.ssl_version = ssl.PROTOCOL_SSLv3
            elif options.ssl_version == "TLSv1":
                options.ssl_version = ssl.PROTOCOL_TLSv1
            else:
                options.ssl_version = ssl.PROTOCOL_SSLv23

        level = SAMPLog.INFO
        stdout = sys.stdout
        stderr = sys.stderr

        if options.loglevel in ("OFF", "ERROR", "WARNING", "DEBUG"):
            if options.loglevel == "OFF":
                level = SAMPLog.OFF
            if options.loglevel == "ERROR":
                level = SAMPLog.ERROR
            if options.loglevel == "WARNING":
                level = SAMPLog.WARNING
            if options.loglevel == "DEBUG":
                level = SAMPLog.DEBUG

        if options.logout != "":
            stdout = open(options.logout, "a")
        if options.logerr != "":
            stderr = open(options.logerr, "a")

        log = SAMPLog(level=level, stdout=stdout, stderr=stderr)

        if BDB_SUPPORT:

            if options.access_restrict != None:
                if options.access_restrict == SAMP_RESTRICT_OWNER and options.owner == "":
                    log.warning("The access cannot be restricted to the owner if the owner "
                                "name is not specified!")
                    sys.exit()

                if options.access_restrict == SAMP_RESTRICT_GROUP and options.owner_group == "":
                    log.warning("The access cannot be restricted to the owner group if the owner "
                                "group name is not specified!")
                    sys.exit()

        args = copy.deepcopy(options.__dict__)
        del(args["loglevel"])
        del(args["logout"])
        del(args["logerr"])
        args["log"] = log

        hub = SAMPHubServer(**args)
        hub.start(False)

        if not timeout:
            while hub.isRunning():
                time.sleep(0.01)
        else:
            time.sleep(timeout)
            hub.stop()

    except KeyboardInterrupt:
        hub.stop()
    except IOError as e:
        print("[SAMP] Error: I/O error({0:%s}): {0:%s}".format(e.errno, e.strerror))
    except SAMPHubError:
        pass
    except SystemExit:
        pass
    except:
        err = StringIO()
        traceback.print_exc(file=err)
        txt = err.getvalue()
        print("[SAMP] Error: Unexpected error:\n" + txt)

import os
import re
import sys

from ...extern.six.moves.urllib.request import urlopen
from ...extern.six.moves import xmlrpc_client as xmlrpc

from .constants import SSL_SUPPORT

if SSL_SUPPORT:
    import ssl


def get_running_hubs():
    """
    Return a dictionary containing the lock-file contents of all the currently
    running hubs (single and/or multiple mode).

    The dictionary format is:

    ``{<lock-file>: {<token-name>: <token-string>, ...}, ...}``

    where ``{<lock-file>}`` is the lock-file name, ``{<token-name>}`` and
    ``{<token-string>}`` are the lock-file tokens (name and content).

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

    hub_is_running, lockfiledict = check_running_hub(lockfilename)

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
            if re.match('samp\\-hub\\-\d+\\-\d+', filename) is not None:
                lockfilename = os.path.join(lockfiledir, filename)
                hub_is_running, lockfiledict = check_running_hub(lockfilename)
                if hub_is_running:
                    hubs[lockfilename] = lockfiledict

    return hubs


def check_running_hub(lockfilename):
    """
    Test whether a hub identified by ``lockfilename`` is running or not.

    Parameters
    ----------
    lockfilename : str
        Lock-file name (path + file name) of the Hub to be tested.

    Returns
    -------
    is_running : bool
        Whether the hub is running
    hub_params : dict
        If the hub is running this contains the parameters from the lockfile
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
    except IOError:
        return is_running, lockfiledict

    for line in lockfile_content:
        if not line.startswith(b"#"):
            kw, val = line.split(b"=")
            lockfiledict[kw.decode().strip()] = val.decode().strip()

    if "samp.hub.xmlrpc.url" in lockfiledict:
        try:
            proxy = xmlrpc.ServerProxy(lockfiledict["samp.hub.xmlrpc.url"]
                                       .replace("\\", ""), allow_none=1)
            proxy.samp.hub.ping()
            is_running = True
        except xmlrpc.ProtocolError:
            # There is a protocol error (e.g. for authentication required),
            # but the server is alive
            is_running = True
        except:
            if SSL_SUPPORT:
                if sys.exc_info()[0] in [ssl.SSLError, ssl.SSLEOFError]:
                    # SSL connection refused for certifcate reasons...
                    # anyway the server is alive
                    is_running = True

    return is_running, lockfiledict


def remove_garbage_lock_files():

    lockfilename = ""

    # HUB SINGLE INSTANCE MODE

    if "HOME" in os.environ:
        # UNIX
        lockfilename = os.path.join(os.environ["HOME"], ".samp")
    else:
        # Windows
        lockfilename = os.path.join(os.environ["USERPROFILE"], ".samp")

    hub_is_running, lockfiledict = check_running_hub(lockfilename)

    if not hub_is_running:
        # If lockfilename belongs to a dead hub, then it is deleted
        if os.path.isfile(lockfilename):
            try:
                os.remove(lockfilename)
            except OSError:
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
            if re.match('samp\\-hub\\-\d+\\-\d+', filename) is not None:
                lockfilename = os.path.join(lockfiledir, filename)
                hub_is_running, lockfiledict = check_running_hub(lockfilename)
                if not hub_is_running:
                    # If lockfilename belongs to a dead hub, then it is deleted
                    if os.path.isfile(lockfilename):
                        try:
                            os.remove(lockfilename)
                        except OSError:
                            pass

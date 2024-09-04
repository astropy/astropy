# Licensed under a 3-clause BSD style license - see LICENSE.rst


# TODO: this file should be refactored to use a more thread-safe and
# race-condition-safe lockfile mechanism.

import datetime
import os
import stat
import warnings
import xmlrpc.client as xmlrpc
from contextlib import suppress
from pathlib import Path
from urllib.parse import urlparse

from astropy import log
from astropy.config.paths import _find_home
from astropy.utils.data import get_readable_fileobj

from .errors import SAMPHubError, SAMPWarning


def read_lockfile(lockfilename):
    """
    Read in the lockfile given by ``lockfilename`` into a dictionary.
    """
    # lockfilename may be a local file or a remote URL, but
    # get_readable_fileobj takes care of this.
    lockfiledict = {}
    with get_readable_fileobj(lockfilename) as f:
        for line in f:
            if not line.startswith("#"):
                kw, val = line.split("=")
                lockfiledict[kw.strip()] = val.strip()
    return lockfiledict


def write_lockfile(lockfilename, lockfiledict):
    lockfilename = Path(lockfilename)
    lockfilename.touch()
    lockfilename.chmod(stat.S_IREAD + stat.S_IWRITE)

    now_iso = datetime.datetime.now().isoformat()

    lockfilename.write_text(
        "".join(
            [
                f"# SAMP lockfile written on {now_iso}\n",
                "# Standard Profile required keys\n",
                *(f"{key}={value}\n" for key, value in lockfiledict.items()),
            ]
        )
    )


def create_lock_file(lockfilename=None, mode=None, hub_id=None, hub_params=None):
    # Remove lock-files of dead hubs
    remove_garbage_lock_files()

    lockfiledir = ""
    homedir = _find_home()

    # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
    if "SAMP_HUB" in os.environ:
        # For the time being I assume just the std profile supported.
        if os.environ["SAMP_HUB"].startswith("std-lockurl:"):
            lockfilename = os.environ["SAMP_HUB"][len("std-lockurl:") :]
            lockfile_parsed = urlparse(lockfilename)

            if lockfile_parsed[0] != "file":
                warnings.warn(
                    f"Unable to start a Hub with lockfile {lockfilename}. "
                    "Start-up process aborted.",
                    SAMPWarning,
                )
                return False
            else:
                lockfilename = lockfile_parsed[2]
    else:
        # If it is a fresh Hub instance
        if lockfilename is None:
            log.debug("Running mode: " + mode)

            if mode == "single":
                lockfilename = homedir / ".samp"
            else:
                lockfiledir = homedir / ".samp-1"

                # If missing create .samp-1 directory
                lockfiledir.mkdir(exist_ok=True)
                lockfiledir.chmod(stat.S_IREAD + stat.S_IWRITE + stat.S_IEXEC)

                lockfilename = lockfiledir / f"samp-hub-{hub_id}"

        else:
            log.debug("Running mode: multiple")

    hub_is_running, lockfiledict = check_running_hub(lockfilename)

    if hub_is_running:
        warnings.warn(
            "Another SAMP Hub is already running. Start-up process aborted.",
            SAMPWarning,
        )
        return False

    log.debug("Lock-file: %s", lockfilename)

    write_lockfile(lockfilename, hub_params)

    return lockfilename


def get_main_running_hub():
    """
    Get either the hub given by the environment variable SAMP_HUB, or the one
    given by the lockfile .samp in the user home directory.
    """
    hubs = get_running_hubs()

    if not hubs:
        raise SAMPHubError("Unable to find a running SAMP Hub.")

    # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
    if (SAMP_HUB := os.getenv("SAMP_HUB")) is not None:
        # For the time being I assume just the std profile supported.
        if not SAMP_HUB.startswith("std-lockurl:"):
            raise SAMPHubError("SAMP Hub profile not supported.")
        lockfilename = Path(SAMP_HUB.removeprefix("std-lockurl:"))
    else:
        lockfilename = _find_home() / ".samp"

    return hubs[lockfilename]


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
    homedir = _find_home()

    # HUB SINGLE INSTANCE MODE

    # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
    if (SAMP_HUB := os.getenv("SAMP_HUB")) is not None:
        # For the time being I assume just the std profile supported.
        if SAMP_HUB.startswith("std-lockurl:"):
            lockfilename = Path(SAMP_HUB.removeprefix("std-lockurl:"))
    else:
        lockfilename = homedir / ".samp"

    hub_is_running, lockfiledict = check_running_hub(lockfilename)

    if hub_is_running:
        hubs[lockfilename] = lockfiledict

    # HUB MULTIPLE INSTANCE MODE

    lockfiledir = homedir / ".samp-1"

    if lockfiledir.is_dir():
        for filename in lockfiledir.glob("*"):
            if not filename.startswith("samp-hub"):
                continue
            lockfilename = filename
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

    # Check whether a lockfile already exists
    try:
        lockfiledict = read_lockfile(lockfilename)
    except OSError:
        return is_running, lockfiledict

    if "samp.hub.xmlrpc.url" in lockfiledict:
        try:
            proxy = xmlrpc.ServerProxy(
                lockfiledict["samp.hub.xmlrpc.url"].replace("\\", ""), allow_none=1
            )
            proxy.samp.hub.ping()
            is_running = True
        except xmlrpc.ProtocolError:
            # There is a protocol error (e.g. for authentication required),
            # but the server is alive
            is_running = True
        except OSError:
            pass

    return is_running, lockfiledict


def remove_garbage_lock_files():
    lockfilename = ""

    # HUB SINGLE INSTANCE MODE

    homedir = _find_home()
    lockfilename = homedir / ".samp"

    hub_is_running, lockfiledict = check_running_hub(lockfilename)

    if not hub_is_running:
        # If lockfilename belongs to a dead hub, then it is deleted
        if lockfilename.is_file():
            with suppress(OSError):
                lockfilename.unlink()

    # HUB MULTIPLE INSTANCE MODE

    lockfiledir = homedir / ".samp-1"

    if not lockfiledir.is_dir():
        return

    for filename in lockfiledir.glob("*"):
        if not filename.startswith("samp-hub"):
            continue

        lockfilename = filename
        hub_is_running, lockfiledict = check_running_hub(lockfilename)
        if not hub_is_running:
            # If lockfilename belongs to a dead hub, then it is deleted
            if lockfilename.is_file():
                with suppress(OSError):
                    lockfilename.unlink()

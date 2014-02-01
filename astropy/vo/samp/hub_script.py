# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy
import time
import sys

from ...utils.compat import argparse

from ... import log, __version__

from .constants import SSL_SUPPORT
from .errors import SAMPHubError
from .hub import SAMPHubServer

if SSL_SUPPORT:
    import ssl

__all__ = ['main']


def hub_script(timeout=0):
    """
    This main function is executed by the ``samp_hub`` command line tool.
    """

    parser = argparse.ArgumentParser(prog="samp_hub " + __version__)

    parser.add_argument("-k", "--secret", dest="secret", metavar="CODE",
                        help="custom secret code.")

    parser.add_argument("-d", "--addr", dest="addr", metavar="ADDR",
                        help="listening address (or IP).")

    parser.add_argument("-p", "--port", dest="port", metavar="PORT", type=int,
                        help="listening port number.")

    parser.add_argument("-f", "--lockfile", dest="lockfile", metavar="FILE",
                        help="custom lockfile.")

    parser.add_argument("-w", "--no-web-profile", dest="web_profile", action="store_false",
                        help="run the Hub disabling the Web Profile.", default=True)

    parser.add_argument("-P", "--pool-size", dest="pool_size", metavar="SIZE", type=int,
                        help="the socket connections pool size.", default=20)

    timeout_group = parser.add_argument_group("Timeout group",
                                              "Special options to setup hub and client timeouts."
                                              "It contains a set of special options that allows to set up the Hub and "
                                              "clients inactivity timeouts, that is the Hub or client inactivity time "
                                              "interval after which the Hub shuts down or unregisters the client. "
                                              "Notification of samp.hub.disconnect MType is sent to the clients "
                                              "forcibly unregistered for timeout expiration.")

    timeout_group.add_argument("-t", "--timeout", dest="timeout", metavar="SECONDS",
                               help="set the Hub inactivity timeout in SECONDS. By default it "
                               "is set to 0, that is the Hub never expires.", type=int, default=0)

    timeout_group.add_argument("-c", "--client-timeout", dest="client_timeout", metavar="SECONDS",
                               help="set the client inactivity timeout in SECONDS. By default it "
                               "is set to 0, that is the client never expires.", type=int, default=0)

    parser.add_argument_group(timeout_group)

    log_group = parser.add_argument_group("Logging options",
                                          "Additional options which allow to customize the logging output. By "
                                          "default the SAMP Hub uses the standard output and standard error "
                                          "devices to print out INFO level logging messages. Using the options "
                                          "here below it is possible to modify the logging level and also "
                                          "specify the output files where redirect the logging messages.")

    log_group.add_argument("-L", "--log-level", dest="loglevel", metavar="LEVEL",
                           help="set the Hub instance log level (OFF, ERROR, WARNING, INFO, DEBUG).",
                           type=str, choices=["OFF", "ERROR", "WARNING", "INFO", "DEBUG"], default='INFO')

    log_group.add_argument("-O", "--log-output", dest="logout", metavar="FILE",
                           help="set the output file for the log messages.", default="")

    parser.add_argument_group(log_group)

    adv_group = parser.add_argument_group("Advanced group",
                                          "Advanced options addressed to facilitate administrative tasks and "
                                          "allow new non-standard Hub behaviors. In particular the --label "
                                          "options is used to assign a value to hub.label token and is used to "
                                          "assign a name to the Hub instance. "
                                          "The very special --multi option allows to start a Hub in multi-instance mode. "
                                          "Multi-instance mode is a non-standard Hub behavior that enables "
                                          "multiple contemporaneous running Hubs. Multi-instance hubs place "
                                          "their non-standard lock-files within the <home directory>/.samp-1 "
                                          "directory naming them making use of the format: "
                                          "samp-hub-<PID>-<ID>, where PID is the Hub process ID while ID is an "
                                          "internal ID (integer).")

    adv_group.add_argument("-l", "--label", dest="label", metavar="LABEL",
                           help="assign a LABEL to the Hub.", default="")

    adv_group.add_argument("-m", "--multi", dest="mode",
                           help="run the Hub in multi-instance mode generating a custom "
                           "lockfile with a random name.",
                           action="store_const", const='multiple', default='single')

    parser.add_argument_group(adv_group)

    if SSL_SUPPORT:

        ssl_group = parser.add_argument_group("SSL group", "Additional options to launch "
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

        ssl_group.add_argument("-s", "--https", dest="https", action="store_true",
                               help="run the Hub using the Secure Sockets Layer.", default=False)

        ssl_group.add_argument("-C", "--cert-file", dest="cert_file", metavar="FILE",
                               help="set the certificate file.", default=None)

        ssl_group.add_argument("-K", "--key-file", dest="key_file", metavar="FILE",
                               help="set the key file. By default this option is ignored, "
                               "assuming that the private key is stored in the certificate file.", default=None)

        ssl_group.add_argument("--cert-reqs", dest="cert_reqs", metavar="STRING",
                               help="this option specifies whether a certificate "
                               "is required from the client side of the connection, and whether "
                               "it will be validated if provided. It must be one of the three "
                               "values NONE (certificates ignored, default), OPTIONAL (not "
                               "required, but validated if provided), or REQUIRED (required "
                               "and validated). If the value of this option is not NONE, "
                               "then the --ca-certs option must point to a file of CA certificates.",
                               type=str, choices=["NONE", "OPTIONAL", "REQUIRED"], default="NONE")

        ssl_group.add_argument("--ca-certs", dest="ca_certs", metavar="FILE",
                               help="the --ca-certs file contains a set of concatenated "
                               "\"certification authority\" certificates, which are used to "
                               "validate certificates passed from the client end of the "
                               "connection.", default=None)

        ssl_group.add_argument("--ssl-version", dest="ssl_version", metavar="STRING",
                               help="the --ssl-version option specifies which version of the "
                               "SSL protocol to use. Typically, the server chooses a particular "
                               "protocol version, and the client must adapt to the server's choice. "
                               "Most of the versions are not interoperable with the other versions. "
                               "If not specified the default SSL version is SSLv23. This version "
                               "provides the most compatibility with other versions client side. "
                               "Other SSL protocol versions are: SSLv2, SSLv3 and TLSv1.",
                               type=str, choices=["SSLv23", "SSLv2", "SSLv3", "TLSv1"], default="SSLv23")

        parser.add_argument_group(ssl_group)

    options = parser.parse_args()

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

        if options.loglevel in ("OFF", "ERROR", "WARNING", "DEBUG", "INFO"):
            log.setLevel(options.loglevel)

        if options.logout != "":
            context = log.log_to_file(options.logout)
        else:
            class dummy_context(object):

                def __enter__(self):
                    pass

                def __exit__(self, exc_type, exc_value, traceback):
                    pass
            context = dummy_context()

        with context:

            args = copy.deepcopy(options.__dict__)
            del(args["loglevel"])
            del(args["logout"])

            hub = SAMPHubServer(**args)
            hub.start(False)

            if not timeout:
                while hub.is_running:
                    time.sleep(0.01)
            else:
                time.sleep(timeout)
                hub.stop()

    except KeyboardInterrupt:
        try:
            hub.stop()
        except NameError:
            pass
    except IOError as e:
        print("[SAMP] Error: I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit(1)
    except SystemExit:
        pass

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import time
import copy
import traceback
import warnings
from optparse import OptionParser, OptionGroup

from ...extern.six import StringIO
from ... import log, __version__

from .constants import (SAMP_HUB_SINGLE_INSTANCE, SAMP_RESTRICT_GROUP,
                        SAMP_RESTRICT_OWNER, SAMP_HUB_MULTIPLE_INSTANCE)
from .errors import SAMPWarning, SAMPHubError
from .hub import SAMPHubServer
from .utils import BDB_SUPPORT, SSL_SUPPORT

if SSL_SUPPORT:
    import ssl

__all__ = ['main']


def main(timeout=0):
    """This main function is executed by the `samp_hub` command line tool."""

    parser = OptionParser(version="%prog " + __version__)

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
                            "default the SAMP Hub uses the standard output and standard error "
                            "devices to print out INFO level logging messages. Using the options "
                            "here below it is possible to modify the logging level and also "
                            "specify the output files where redirect the logging messages.")

    log_group.add_option("-L", "--log-level", dest="loglevel", metavar="LEVEL",
                         help="set the Hub instance log level (OFF, ERROR, WARNING, INFO, DEBUG).",
                         type="choice", choices=["OFF", "ERROR", "WARNING", "INFO", "DEBUG"])

    log_group.add_option("-O", "--log-output", dest="logout", metavar="FILE",
                         help="set the output file for the log messages.")

    # The following was implemented in SAMPy but is not easy to implement here
    # without modifying the Astropy logger.
    # log_group.add_option("-E", "--log-error", dest="logerr", metavar="FILE",
    #                      help="set the output file for WARNING and ERROR messages.")

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

            if BDB_SUPPORT:

                if options.access_restrict != None:
                    if options.access_restrict == SAMP_RESTRICT_OWNER and options.owner == "":
                        warnings.warn("The access cannot be restricted to the owner if the owner "
                                      "name is not specified!", SAMPWarning)
                        sys.exit()

                    if options.access_restrict == SAMP_RESTRICT_GROUP and options.owner_group == "":
                        warnings.warn("The access cannot be restricted to the owner group if the owner "
                                      "group name is not specified!", SAMPWarning)
                        sys.exit()

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

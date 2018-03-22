# Licensed under a 3-clause BSD style license - see LICENSE.rst


import copy
import time
import sys
import argparse

from .. import log, __version__

from .hub import SAMPHubServer

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

    options = parser.parse_args()

    try:

        if options.loglevel in ("OFF", "ERROR", "WARNING", "DEBUG", "INFO"):
            log.setLevel(options.loglevel)

        if options.logout != "":
            context = log.log_to_file(options.logout)
        else:
            class dummy_context:

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
    except OSError as e:
        print("[SAMP] Error: I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit(1)
    except SystemExit:
        pass

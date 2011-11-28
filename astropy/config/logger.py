# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module defines a logging class based on the built-in logging module"""

from __future__ import division

import logging

from . import ConfigurationItem

LEVEL = ConfigurationItem('log_level', logging.INFO,
                          'Threshold for the logging messages. Logging messages which are less severe than this level will be ignored. The levels are 10 (DEBUG), 20 (INFO), 30 (WARNING), 40 (ERROR), and 50 (CRITICAL)')

COLOR = ConfigurationItem('log_color', True,
                          'Whether to color-code the logging messages according to severity.')


def add_coloring_to_emit_ansi(fn):
    def new(*args):
        levelno = args[1].levelno
        if(levelno >= 50):
            color = '\x1b[31m'  # red
        elif(levelno >= 40):
            color = '\x1b[31m'  # red
        elif(levelno >= 30):
            color = '\x1b[33m'  # yellow
        elif(levelno >= 20):
            color = '\x1b[32m'  # green
        elif(levelno >= 10):
            color = '\x1b[35m'  # pink
        else:
            color = '\x1b[0m'  # normal
        args[1].levelname = color + args[1].levelname + '\x1b[0m'  # normal
        return fn(*args)
    return new

# Initialize logger
logging.basicConfig(format="%(levelname)s: %(message)s", level=LEVEL())
logger = logging.getLogger()

if COLOR():
    f = logging.Formatter("%(levelname)s: %(message)s")
    if len(logger.handlers) > 0:
        logger.handlers[0].setFormatter(f)
    logging.StreamHandler.emit = add_coloring_to_emit_ansi(logging.StreamHandler.emit)

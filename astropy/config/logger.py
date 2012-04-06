# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module defines a logging class based on the built-in logging module"""

import sys
import logging
import warnings

from . import ConfigurationItem
from ..utils.console import color_print

# Save original showwarning and excepthook functions so they can be restored
_showwarning = warnings.showwarning
_excepthook = sys.excepthook


def set_catch_warnings(catch):
    if catch:
        def logging_showwarning(*args, **kwargs):
            logger.warn(args[0].message)
        warnings.showwarning = logging_showwarning
    else:
        warnings.showwarning = _showwarning


def set_catch_exceptions(catch):
    if catch:
        def handle_exceptions(type, value, exception):
            print type(value)
            logger.error(value.message)
            _excepthook(type, value, exception)
        sys.excepthook = handle_exceptions
    else:
        sys.excepthook = _excepthook

# Set default levels
DEFAULT_LEVELS = {'debug': 10, 'info': 20, 'warn': 30, 'error': 40}


def numerical_level(log_level):
    if log_level in DEFAULT_LEVELS:
        return DEFAULT_LEVELS[log_level]
    else:
        raise ValueError("{:s} is not a valid log level (should be one of "
                         "debug/info/warning/error)".format(log_level))

# Read in configuration

log_level = ConfigurationItem('log_level', DEFAULT_LEVELS['info'],
                              "Threshold for the logging messages. Logging "
                              "messages that are less severe than this level "
                              "will be ignored. The levels are 'debug', "
                              "'info', 'warning', 'error'")

log_level = numerical_level(log_level)

use_color = ConfigurationItem('use_color', True,
                              "Whether to use color for the level names")

catch_warnings = ConfigurationItem('catch_warnings', False,
                                   "Whether to catch warnings.warn calls and "
                                   "output them via the logger")

set_catch_warnings(catch_warnings)
del catch_warnings  # future changes should be done directly by the function

catch_exceptions = ConfigurationItem('catch_exceptions', False,
                                     "Whether to output an entry for "
                                     "exceptions in the logger")

set_catch_exceptions(catch_exceptions)
del catch_exceptions  # future changes should be done directly by the function


# Initialize logger
logger = logging.getLogger('astropy')

# Set up the stdout handler
sh = logging.StreamHandler()


def stream_formatter(record):
    if record.levelno < 10 or not use_color:
        print(record.levelname),
    elif(record.levelno < 20):
        color_print(record.levelname, 'pink', end='')
    elif(record.levelno < 30):
        color_print(record.levelname, 'green', end='')
    elif(record.levelno < 40):
        color_print(record.levelname, 'yellow', end='')
    else:
        color_print(record.levelname, 'red', end='')
    print(": " + record.msg)

sh.emit = stream_formatter
logger.addHandler(sh)


class FilterOrigin(object):
    '''A filter for the record origin'''
    def __init__(self, origin):
        self.origin = origin

    def filter(self, record):
        return record.module.startswith(self.origin)


class log_to_file(object):
    '''A context manager to log to a file'''

    def __init__(self, filename, filter_level='info', filter_origin=None):
        self.filename = filename
        self.filter_level = numerical_level(filter_level)
        self.filter_origin = filter_origin

    def __enter__(self):
        self.fh = logging.FileHandler(self.filename)
        self.fh.setLevel(self.filter_level)
        self.fh.addFilter(FilterOrigin(self.filter_origin))
        logger.addHandler(self.fh)

    def __exit__(self, *exc_info):
        logger.removeHandler(self.fh)


class ListHandler(logging.Handler):
    '''A handler that can be used to capture the records in a list'''

    def __init__(self, filter_level=None, filter_origin=None):
        logging.Handler.__init__(self)
        self.log_list = []

    def emit(self, record):
        self.log_list.append(record)


class log_to_list(object):
    '''A context manager to log to a list'''

    def __init__(self, filter_level='info', filter_origin=None):
        self.filter_level = numerical_level(filter_level)
        self.filter_origin = filter_origin

    def __enter__(self):
        self.lh = ListHandler()
        self.lh.setLevel(self.filter_level)
        self.lh.addFilter(FilterOrigin(self.filter_origin))
        logger.addHandler(self.lh)
        return self.lh.log_list

    def __exit__(self, *exc_info):
        logger.removeHandler(self.lh)

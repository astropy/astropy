# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module defines a logging class based on the built-in logging module"""

from __future__ import print_function

import os
import sys
import inspect
import logging
import warnings
from contextlib import contextmanager

from . import ConfigurationItem
from ..utils.console import color_print
from ..utils.misc import find_current_module

# Save original showwarning and excepthook functions so they can be restored
_showwarning = warnings.showwarning
_excepthook = sys.excepthook

FILE_FORMAT = "%(asctime)s - %(origin)s - %(levelname)s - %(message)s"

class FilterOrigin(object):
    '''A filter for the record origin'''
    def __init__(self, origin):
        self.origin = origin

    def filter(self, record):
        return record.module.startswith(self.origin)


class ListHandler(logging.Handler):
    '''A handler that can be used to capture the records in a list'''

    def __init__(self, filter_level=None, filter_origin=None):
        logging.Handler.__init__(self)
        self.log_list = []

    def emit(self, record):
        self.log_list.append(record)

Logger = logging.getLoggerClass()


class AstropyLogger(Logger):

    def debug(self, *args, **kwargs):
        if 'extra' not in kwargs:
            kwargs['extra'] = {}
        if 'origin' not in kwargs['extra']:
            kwargs['extra']['origin'] = find_current_module(2).__name__
        Logger.debug(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        if 'extra' not in kwargs:
            kwargs['extra'] = {}
        if 'origin' not in kwargs['extra']:
            kwargs['extra']['origin'] = find_current_module(2).__name__
        Logger.info(self, *args, **kwargs)

    def warning(self, *args, **kwargs):
        if 'extra' not in kwargs:
            kwargs['extra'] = {}
        if 'origin' not in kwargs['extra']:
            kwargs['extra']['origin'] = find_current_module(2).__name__
        Logger.warning(self, *args, **kwargs)

    warn = warning

    def error(self, *args, **kwargs):
        if 'extra' not in kwargs:
            kwargs['extra'] = {}
        if 'origin' not in kwargs['extra']:
            kwargs['extra']['origin'] = find_current_module(2).__name__
        Logger.error(self, *args, **kwargs)

    def exception(self, *args, **kwargs):
        if 'extra' not in kwargs:
            kwargs['extra'] = {}
        if 'origin' not in kwargs['extra']:
            kwargs['extra']['origin'] = find_current_module(2).__name__
        Logger.exception(self, *args, **kwargs)

    def critical(self, *args, **kwargs):
        if 'extra' not in kwargs:
            kwargs['extra'] = {}
        if 'origin' not in kwargs['extra']:
            kwargs['extra']['origin'] = find_current_module(2).__name__
        Logger.critical(self, *args, **kwargs)

    fatal = critical

    def log(self, *args, **kwargs):
        if 'extra' not in kwargs:
            kwargs['extra'] = {}
        kwargs['extra']['origin'] = find_current_module(2).__name__
        Logger.log(self, *args, **kwargs)

    def set_catch_warnings(self, catch):
        if catch:
            def logging_showwarning(*args, **kwargs):
                self.warn(args[0].message)
            warnings.showwarning = logging_showwarning
        else:
            warnings.showwarning = _showwarning

    def set_catch_exceptions(self, catch):
        if catch:
            def handle_exceptions(type, value, exception):
                try:
                    origin = inspect.getmodule(exception.tb_next).__name__
                except:
                    origin = inspect.getmodule(exception).__name__
                self.error(value.message, extra={'origin': origin})
                _excepthook(type, value, exception)
            sys.excepthook = handle_exceptions
        else:
            sys.excepthook = _excepthook

    def setColor(self, use_color):
        self._use_color = use_color

    def stream_formatter(self, record):
        if record.levelno < logging.DEBUG or not self._use_color:
            print(record.levelname, end='')
        elif(record.levelno < logging.INFO):
            color_print(record.levelname, 'magenta', end='')
        elif(record.levelno < logging.WARN):
            color_print(record.levelname, 'green', end='')
        elif(record.levelno < logging.ERROR):
            color_print(record.levelname, 'brown', end='')
        else:
            color_print(record.levelname, 'red', end='')
        print(": " + record.msg + " [{:s}]".format(record.origin))

    @contextmanager
    def log_to_file(self, filename, filter_level='INFO', filter_origin=None):
        fh = logging.FileHandler(filename)
        fh.setLevel(filter_level)
        if filter_origin is not None:
            fh.addFilter(FilterOrigin(filter_origin))
        f = logging.Formatter(FILE_FORMAT)
        fh.setFormatter(f)
        self.addHandler(fh)
        yield
        self.removeHandler(fh)

    @contextmanager
    def log_to_list(self, filter_level='INFO', filter_origin=None):
        lh = ListHandler()
        lh.setLevel(filter_level)
        if filter_origin is not None:
            lh.addFilter(FilterOrigin(filter_origin))
        self.addHandler(lh)
        yield lh.log_list
        self.removeHandler(lh)

logging.setLoggerClass(AstropyLogger)

# Read in configuration

LOG_LEVEL = ConfigurationItem('log_level', 'WARN',
                              "Threshold for the logging messages. Logging "
                              "messages that are less severe than this level "
                              "will be ignored. The levels are 'DEBUG', "
                              "'INFO', 'WARNING', 'ERROR'")

USE_COLOR = ConfigurationItem('use_color', True,
                              "Whether to use color for the level names")

CATCH_WARNINGS = ConfigurationItem('catch_warnings', False,
                                   "Whether to catch warnings.warn calls and "
                                   "output them via the logger")

CATCH_EXCEPTIONS = ConfigurationItem('catch_exceptions', False,
                                     "Whether to output an entry for "
                                     "exceptions in the logger")

LOG_TO_FILE = ConfigurationItem('log_to_file', True,
                                "Whether to always log messages to a log "
                                "file")

LOG_FILE_PATH = ConfigurationItem('log_file_path', '~/.astropy/astropy.log',
                                  "The file to log messages to")

LOG_FILE_LEVEL = ConfigurationItem('log_file_level', 'WARN',
                                   "Threshold for logging messages to "
                                   "log_file_path")


# Initialize logger
logger = logging.getLogger('astropy')
logger.setLevel(LOG_LEVEL())
logger.setColor(USE_COLOR())

# Set up the stdout handler
sh = logging.StreamHandler()
sh.emit = logger.stream_formatter
logger.addHandler(sh)

# Set up the main log file handler if requested (but this might fail if
# configuration directory or log file is not writeable).
if LOG_TO_FILE():
    try:
        fh = logging.FileHandler(os.path.expanduser(LOG_FILE_PATH()))
    except IOError:
        pass
    else:
        formatter = logging.Formatter(FILE_FORMAT)
        fh.setFormatter(formatter)
        fh.setLevel(LOG_FILE_LEVEL())
        logger.addHandler(fh)

logger.set_catch_warnings(CATCH_WARNINGS())
logger.set_catch_exceptions(CATCH_EXCEPTIONS())
